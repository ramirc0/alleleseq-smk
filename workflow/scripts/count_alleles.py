"""Count reads per allele at each het SNV from haplotype pileups.

Port of AlleleSeq2's ``pileup2counts.py``. Each het SNV has one position on
each haplotype. Counts from both haplotypes are summed per base, the
haplotype bases in the pileups decide which allele sits on which haplotype,
and a binomial test compares the reference allele to 0.5.
"""

import argparse
import sys

from binomial import binomial_test
from pileup import BASES, SiteCounts, pileup_counts

HEADER = [
    "#chr",
    "ref_coord",
    "hap1_coord",
    "hap2_coord",
    "ref_allele",
    "hap1_allele",
    "hap2_allele",
    "cA",
    "cC",
    "cG",
    "cT",
    "cN",
    "ref_allele_ratio",
    "sum_ref_n_alt_cnts",
    "p_binom",
    "warn_hap1",
    "warn_hap2",
]
NOT_LIFTED = "not_lifted?"


def build_parser():
    """Build the count_alleles.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the count_alleles.py flags.
    """
    p = argparse.ArgumentParser(description="Count alleles at het SNVs.")
    p.add_argument(
        "--ref-bed", required=True, help="Het SNVs in reference coordinates."
    )
    p.add_argument("--hap1-bed", required=True, help="Het SNVs lifted to hap1.")
    p.add_argument("--hap2-bed", required=True, help="Het SNVs lifted to hap2.")
    p.add_argument("--discarded", required=True, help="Log of discarded SNVs to write.")
    p.add_argument("--output", required=True, help="Counts TSV to write.")
    p.add_argument(
        "--pileups", nargs="+", required=True, help="Unique-read mpileups, any order."
    )
    return p


def main(argv=None):
    """Write the raw allele counts table.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with open(args.discarded, "w") as discarded, open(args.output, "w") as out:
        snvs = read_lifted_snvs(args.hap1_bed, args.hap2_bed, discarded)
        sites = pileup_counts(args.pileups)
        out.write("\t".join(HEADER) + "\n")
        with open(args.ref_bed) as ref_bed:
            for line in ref_bed:
                fields = line.split("\t")
                key = f"{fields[0]}_{fields[2]}"
                if key in snvs:
                    row = count_row(key, snvs[key], sites, discarded)
                    if row is not None:
                        out.write("\t".join(row) + "\n")


def read_lifted_snvs(hap1_bed, hap2_bed, discarded):
    """Collect each het SNV's haplotype positions and alleles.

    Parameters
    ----------
    hap1_bed : str
        Het SNVs lifted to hap1. The name column is
        ``<chr>_<pos>_<ref>_<allele1>_<allele2>``.
    hap2_bed : str
        Het SNVs lifted to hap2.
    discarded : file
        Log for SNVs with non-ACGTN alleles or shared positions.

    Returns
    -------
    dict of str to dict
        Keyed by ``<chr>_<pos>``. Values hold ``hap1_pos``, ``hap2_pos`` (when
        lifted), ``r_a``, ``hap1_a`` and ``hap2_a``.
    """
    snvs = {}
    invalid = set()
    duplicates = {}
    for hap, path in (("hap1", hap1_bed), ("hap2", hap2_bed)):
        with open(path) as bed:
            for line in bed:
                h_chr, _, h_crd, ref_info = line.strip().split("\t")
                r_chr, r_crd, r_a, a1, a2 = ref_info.rsplit("_", 4)
                key = f"{r_chr}_{r_crd}"
                if a1 not in BASES or a2 not in BASES:
                    if hap == "hap1" or key not in invalid:
                        discarded.write(f"{r_chr}\t{r_crd}\t{ref_info}\n")
                    if hap == "hap1":
                        invalid.add(key)
                    continue
                snv = snvs.setdefault(key, {})
                if f"{hap}_pos" in snv or key in duplicates:
                    duplicates[key] = None
                    continue
                snv[f"{hap}_pos"] = f"{h_chr}_{h_crd}"
                snv.update(r_a=r_a, hap1_a=a1, hap2_a=a2)
    for key in duplicates:
        del snvs[key]
        chrom, pos = key.rsplit("_", 1)
        print(
            f"WARN: same coords of two SNVs: {chrom} {pos} ? Removing.", file=sys.stderr
        )
        discarded.write(f"{chrom}\t{pos}\tduplicate_pos\n")
    return snvs


def count_row(key, snv, sites, discarded):
    """Build the counts row of one het SNV.

    Parameters
    ----------
    key : str
        ``<chr>_<pos>`` of the SNV.
    snv : dict
        Haplotype positions and alleles, from `read_lifted_snvs`.
    sites : dict of str to SiteCounts
        Pooled pileup counts by haplotype position.
    discarded : file
        Log for SNVs whose alleles disagree with the haplotypes.

    Returns
    -------
    list of str or None
        The row, or None when the SNV has no reads or is discarded.
    """
    absent = SiteCounts(dict.fromkeys(BASES, 0), ".", ".")
    hap1 = sites.get(snv.get("hap1_pos"), absent)
    hap2 = sites.get(snv.get("hap2_pos"), absent)
    counts = {base: hap1.counts[base] + hap2.counts[base] for base in BASES}
    a1, a2 = snv["hap1_a"], snv["hap2_a"]
    total = counts[a1] + counts[a2]
    if total < 1:
        return None
    chrom, pos = key.rsplit("_", 1)
    expected = [(a1, a2), (a1, "."), (".", a2)]
    if (hap1.allele, hap2.allele) in expected:
        hap1_allele, hap2_allele = a1, a2
    elif (hap2.allele, hap1.allele) in expected:
        hap1_allele, hap2_allele = a2, a1
    else:
        print(
            f"WARN: confused with phasing in .vcf vs .mpileup, skipping: {key}",
            file=sys.stderr,
        )
        discarded.write(f"{chrom}\t{pos}\tvcf_vs_mpileup_alleles\n")
        return None
    ref = counts[snv["r_a"]]
    if ref > total:
        discarded.write(
            f"{chrom}\t{pos}\tref_allele:{ref}__hap1:{counts[a1]}__hap2:{counts[a2]}\n"
        )
        return None
    return [
        chrom,
        pos,
        snv.get("hap1_pos", NOT_LIFTED),
        snv.get("hap2_pos", NOT_LIFTED),
        snv["r_a"],
        hap1_allele,
        hap2_allele,
        *(str(counts[base]) for base in BASES),
        repr(ref / total),
        str(total),
        repr(binomial_test(ref, total)),
        hap1.warning,
        hap2.warning,
    ]


if __name__ == "__main__":
    main()
