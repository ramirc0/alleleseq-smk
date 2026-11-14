"""Normalize WGS read counts around het SNVs by their median.

Port of AlleleSeq2's ``calculate_rd.py``. The input is ``bedtools intersect
-c`` output: each window followed by its read count. The median uses only
windows with reads, as upstream did. Windows without reads get a ratio of 0.
"""

import argparse
import statistics
import sys


def build_parser():
    """Build the read_depth.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the read_depth.py flags.
    """
    p = argparse.ArgumentParser(description="Read depth ratio per het SNV window.")
    p.add_argument("--output", required=True, help="Read depth table to write.")
    return p


def main(argv=None):
    """Read window counts from stdin and write the depth ratios.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    counts = {}
    for line in sys.stdin:
        fields = line.rstrip("\n").split("\t")
        counts[fields[3]] = counts.get(fields[3], 0) + int(fields[4])
    median = statistics.median(c for c in counts.values() if c > 0)
    print(f"median rd:\t{median}", file=sys.stderr)
    with open(args.output, "w") as out:
        out.write("#chr\tpos\trd\n")
        for name, count in counts.items():
            chrom, pos = name.split("_")[:2]
            out.write(f"{chrom}\t{pos}\t{count / median!r}\n")


if __name__ == "__main__":
    main()
