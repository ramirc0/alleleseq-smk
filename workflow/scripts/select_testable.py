"""Keep the rows that enter the allele-specific tests.

Port of AlleleSeq2's ``filter_non-autosomal_chr.py`` and
``filter_by_counts.py`` and their region twins. A row is testable when it
lies on an autosome or a kept chromosome and its two alleles (or haplotypes)
reach the total and per-allele read thresholds.
"""

import argparse

from chromosomes import is_testable


def build_parser():
    """Build the select_testable.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the select_testable.py flags.
    """
    p = argparse.ArgumentParser(description="Select testable SNVs or regions.")
    p.add_argument("--counts", required=True, help="Filtered SNV or region counts TSV.")
    p.add_argument("--output", required=True, help="Testable counts TSV to write.")
    p.add_argument(
        "--min-total", type=int, default=6, help="Minimum reads on both alleles."
    )
    p.add_argument(
        "--min-per-allele", type=int, default=0, help="Minimum reads per allele."
    )
    p.add_argument(
        "--keep-chromosomes", nargs="*", default=[], help="Extra chromosomes, e.g. X."
    )
    p.add_argument(
        "--regions", help="BED4 of the regions. Marks the counts as region-level."
    )
    return p


def main(argv=None):
    """Write the testable rows.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    region_chromosome = read_region_chromosomes(args.regions) if args.regions else None
    with open(args.counts) as counts, open(args.output, "w") as out:
        out.write(counts.readline())
        for line in counts:
            row = line.rstrip("\n").split("\t")
            if region_chromosome is None:
                chrom = row[0]
                bases = dict(zip("ACGTN", map(int, row[7:12])))
                allele1, allele2 = bases[row[5]], bases[row[6]]
            else:
                chrom = region_chromosome[row[0]]
                allele1, allele2 = int(row[1]), int(row[2])
            if (
                is_testable(chrom, args.keep_chromosomes)
                and allele1 + allele2 >= args.min_total
                and allele1 >= args.min_per_allele
                and allele2 >= args.min_per_allele
            ):
                out.write(line)


def read_region_chromosomes(path):
    """Map each region name to its chromosome.

    Parameters
    ----------
    path : str
        BED4 file of regions.

    Returns
    -------
    dict of str to str
        Chromosome by region name. A repeated name keeps its last line.
    """
    with open(path) as bed:
        return {
            fields[3]: fields[0]
            for fields in (line.rstrip("\n").split("\t") for line in bed)
        }


if __name__ == "__main__":
    main()
