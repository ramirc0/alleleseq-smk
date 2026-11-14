"""Filter het SNV counts for copy number, phase warnings and multi-mapping reads.

Port of AlleleSeq2's ``filter_cnv_sites.py``, ``filter_phase_warnings.py``
and ``filter_sites_w_mmaps.py adjust``, applied to each row in that order.

1. Keep SNVs whose WGS read depth ratio lies in [0.5, 1.5].
2. Drop SNVs with any pileup warning or a missing haplotype position.
3. Add multi-mapping reads that carry the weaker allele to its count, up to
   the stronger allele's count, so that unaccounted multi-mapping reads
   cannot create the imbalance.
"""

import argparse

from binomial import binomial_test
from pileup import BASES, SiteCounts, pileup_counts

RD_MIN, RD_MAX = 0.5, 1.5


def build_parser():
    """Build the filter_counts.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the filter_counts.py flags.
    """
    p = argparse.ArgumentParser(description="Filter raw het SNV counts.")
    p.add_argument(
        "--counts", required=True, help="Raw counts TSV from count_alleles.py."
    )
    p.add_argument("--read-depth", required=True, help="WGS read depth table.")
    p.add_argument(
        "--mmap-pileups", nargs="+", required=True, help="Multi-mapping mpileups."
    )
    p.add_argument("--output", required=True, help="Filtered counts TSV to write.")
    p.add_argument(
        "--cnv-log", required=True, help="Log of SNVs failing the depth filter."
    )
    p.add_argument(
        "--warning-log", required=True, help="Log of SNVs with phase warnings."
    )
    p.add_argument(
        "--mmap-log", required=True, help="Log of multi-mapping adjustments."
    )
    return p


def main(argv=None):
    """Write the filtered counts table and the three logs.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    depth = read_depth_ratios(args.read_depth)
    mmap = pileup_counts(args.mmap_pileups)
    with (
        open(args.counts) as counts,
        open(args.output, "w") as out,
        open(args.cnv_log, "w") as cnv_log,
        open(args.warning_log, "w") as warning_log,
        open(args.mmap_log, "w") as mmap_log,
    ):
        cnv_log.write("#chr\tref_coord\tcA_cC_cG_cT__rd\n")
        warning_log.write("#chr\tref_coord\tcA_cC_cG_cT__warning_hap1__warning_hap2\n")
        mmap_log.write(
            "#chr\tpos\t(max_mm):cA_cC_cG_cT\tmm_hap1_warn;mm_hap2_warn;mm_log\n"
        )
        header = counts.readline().rstrip("\n").split("\t")
        out.write("\t".join(header[:15] + ["cnv", "mmap_log"]) + "\n")
        for line in counts:
            row = line.rstrip("\n").split("\t")
            rd = passes_depth(row, depth, cnv_log)
            if rd is None or not passes_warnings(row, warning_log):
                continue
            out.write(
                "\t".join(adjust_for_mmaps(row[:15] + [repr(rd)], mmap, mmap_log))
                + "\n"
            )


def read_depth_ratios(path):
    """Read WGS read depth ratios by SNV.

    Parameters
    ----------
    path : str
        Table with ``#chr``, ``pos`` and ``rd`` columns.

    Returns
    -------
    dict of str to float
        Ratio keyed by ``<chr>_<pos>``.
    """
    ratios = {}
    with open(path) as handle:
        for line in handle:
            if not line.startswith("#"):
                chrom, pos, rd = line.rstrip("\n").split("\t")
                ratios[f"{chrom}_{pos}"] = float(rd)
    return ratios


def passes_depth(row, depth, log):
    """Apply the copy number filter to one row.

    Parameters
    ----------
    row : list of str
        Raw counts row.
    depth : dict of str to float
        WGS read depth ratio by SNV.
    log : file
        Log for SNVs that fail. A SNV absent from `depth` is logged as ``NA``.

    Returns
    -------
    float or None
        The depth ratio when the SNV passes, else None.
    """
    chrom, coord = row[0], row[1]
    rd = depth.get(f"{chrom}_{coord}")
    if rd is not None and RD_MIN <= rd <= RD_MAX:
        return rd
    shown = "NA" if rd is None else repr(rd)
    log.write(f"{chrom}\t{coord}\t{'_'.join(row[7:11])}__{shown}\n")
    return None


def passes_warnings(row, log):
    """Apply the phase warning filter to one row.

    Parameters
    ----------
    row : list of str
        Raw counts row.
    log : file
        Log for SNVs that fail.

    Returns
    -------
    bool
        True when both haplotypes are clean and lifted.
    """
    hap1_coord, hap2_coord, warn1, warn2 = row[2], row[3], row[15], row[16]
    clean = all(w == "." for w in warn1.split(",")) and all(
        w == "." for w in warn2.split(",")
    )
    if clean and hap1_coord != "not_lifted?" and hap2_coord != "not_lifted?":
        return True
    log.write(f"{row[0]}\t{row[1]}\t{'_'.join(row[7:11])}__{warn1}__{warn2}\n")
    return False


def adjust_for_mmaps(row, mmap, log):
    """Add weaker-allele multi-mapping reads to one row.

    Parameters
    ----------
    row : list of str
        Counts row ending with the ``cnv`` column.
    mmap : dict of str to SiteCounts
        Multi-mapping pileup counts by haplotype position.
    log : file
        Log of the multi-mapping counts and the action taken.

    Returns
    -------
    list of str
        The row with a trailing ``mmap_log`` column.
    """
    absent = SiteCounts(dict.fromkeys(BASES, 0), ".", ".")
    hap1 = mmap.get(row[2], absent)
    hap2 = mmap.get(row[3], absent)
    mm = {base: max(hap1.counts[base], hap2.counts[base]) for base in BASES}
    ref_allele, hap1_allele, hap2_allele = row[4], row[5], row[6]
    unique = {base: int(row[7 + i]) for i, base in enumerate(BASES)}
    mm_text = "_".join(str(mm[base]) for base in "ACGT")
    warnings = f"{hap1.warning};{hap2.warning}"

    if unique[hap1_allele] == unique[hap2_allele]:
        log.write(f"{row[0]}\t{row[1]}\t{mm_text}\t{warnings};eq_uniq_cnts\n")
        return row + ["."]
    weaker = hap2_allele if unique[hap1_allele] > unique[hap2_allele] else hap1_allele
    if mm[weaker] == 0:
        log.write(f"{row[0]}\t{row[1]}\t{mm_text}\t{warnings};{weaker}:0\n")
        return row + [f"{weaker}:0"]

    stronger = max(unique[hap1_allele], unique[hap2_allele])
    adjusted = min(unique[weaker] + mm[weaker], stronger)
    added = adjusted - unique[weaker]
    unique[weaker] = adjusted
    total = unique[hap1_allele] + unique[hap2_allele]
    ref = unique[ref_allele]
    log.write(f"{row[0]}\t{row[1]}\t{mm_text}\t{warnings};{weaker}:+{added}\n")
    return [
        *row[:7],
        *(str(unique[base]) for base in BASES),
        repr(ref / total),
        str(total),
        repr(binomial_test(ref, total)),
        row[15],
        f"{weaker}:+{added}",
    ]


if __name__ == "__main__":
    main()
