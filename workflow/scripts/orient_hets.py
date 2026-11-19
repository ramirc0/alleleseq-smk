"""Phase every genotype so that ``bcftools consensus -H 1/-H 2`` can split them.

Phased genotypes without a phase set (``PS``) keep their order. A phase set,
such as GATK's physical phasing, only fixes phase within its block, so each
block gets one random orientation. Each unphased het gets a random
orientation from a seeded generator, as the AlleleSeq paper does for hets the
parents cannot phase. The orientations and the phased fraction are written
out, because gene-level aggregation is only meaningful for phased hets.
"""

import argparse
import gzip
import io
import random


def build_parser():
    """Build the orient_hets.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the orient_hets.py flags.
    """
    p = argparse.ArgumentParser(description="Orient unphased hets at random.")
    p.add_argument("--vcf", required=True, help="Single-sample VCF, plain or gzipped.")
    p.add_argument("--seed", type=int, default=0, help="Seed for the orientations.")
    p.add_argument(
        "--output", required=True, help="Oriented VCF to write (plain text)."
    )
    p.add_argument(
        "--orientations", required=True, help="Gzipped TSV of oriented hets."
    )
    p.add_argument("--summary", required=True, help="TSV with the phased fraction.")
    return p


def main(argv=None):
    """Write the oriented VCF, the orientation table and the summary.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    rng = random.Random(args.seed)
    flips = {}
    hets = phased = 0
    opener = gzip.open if args.vcf.endswith(".gz") else open
    with (
        opener(args.vcf, "rt") as vcf,
        open(args.output, "w") as out,
        open_deterministic(args.orientations) as table,
    ):
        table.write("chrom\tpos\tref\talt\tgenotype\toriented\n")
        for line in vcf:
            if line.startswith("#"):
                out.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            genotype, oriented, chromosome_wide = orient(fields, rng, flips)
            alleles = oriented.split("|")
            if len(alleles) == 2 and alleles[0] != alleles[1]:
                hets += 1
                if chromosome_wide:
                    phased += 1
                else:
                    row = [*fields[:2], fields[3], fields[4], genotype, oriented]
                    table.write("\t".join(row) + "\n")
            out.write("\t".join(fields) + "\n")
    with open(args.summary, "w") as out:
        fraction = phased / hets if hets else float("nan")
        out.write("het_sites\tphased_in_vcf\tfraction_phased\n")
        out.write(f"{hets}\t{phased}\t{fraction:.4f}\n")


def open_deterministic(path):
    """Open a gzipped text file without a timestamp, so reruns are byte-identical.

    Parameters
    ----------
    path : str
        Output path.

    Returns
    -------
    io.TextIOWrapper
        Text stream over a gzip file with mtime 0.
    """
    return io.TextIOWrapper(gzip.GzipFile(path, "wb", mtime=0))


def orient(fields, rng, flips):
    """Phase the genotype of one VCF record in place.

    Parameters
    ----------
    fields : list of str
        VCF record columns. The sample is the tenth column.
    rng : random.Random
        Source of orientations for unphased hets and phase sets.
    flips : dict
        Orientation drawn for each ``(chrom, PS)`` phase set so far. Updated
        in place.

    Returns
    -------
    genotype : str
        The genotype as read.
    oriented : str
        The phased genotype written back into `fields`.
    chromosome_wide : bool
        Whether the genotype was phased along the whole chromosome.
    """
    keys = fields[8].split(":")
    values = fields[9].split(":")
    calls = dict(zip(keys, values, strict=False))
    genotype = calls["GT"]
    alleles = genotype.replace("/", "|").split("|")
    het = len(alleles) == 2 and alleles[0] != alleles[1]
    phase_set = calls.get("PS", ".")
    chromosome_wide = "|" in genotype and phase_set == "."
    if het and not chromosome_wide:
        if "|" in genotype:
            block = (fields[0], phase_set)
            if block not in flips:
                flips[block] = rng.random() < 0.5
            flip = flips[block]
        else:
            flip = rng.random() < 0.5
        if flip:
            alleles.reverse()
    oriented = "|".join(alleles)
    values[keys.index("GT")] = oriented
    fields[9] = ":".join(values)
    return genotype, oriented, chromosome_wide


if __name__ == "__main__":
    main()
