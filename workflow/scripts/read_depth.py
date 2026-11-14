"""Normalize WGS read counts around het SNVs by their median.

Port of AlleleSeq2's ``calculate_rd.py``. The input is upstream's ``bedtools
intersect -a <reads> -b <windows> -sorted -wb -bed`` output: one line per
read and window overlap, with the window name last. Streaming the reads keeps
memory bounded by the windows, even over a pileup of millions of reads. The
median uses only windows with reads, as upstream did. Windows without reads
get a ratio of 0.
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
    p.add_argument(
        "--windows", required=True, help="Window BED with SNV names in column 4."
    )
    p.add_argument("--output", required=True, help="Read depth table to write.")
    return p


def main(argv=None):
    """Count overlaps from stdin per window and write the depth ratios.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with open(args.windows) as windows:
        counts = {line.split("\t")[3].rstrip("\n"): 0 for line in windows}
    for line in sys.stdin:
        counts[line.rstrip("\n").rsplit("\t", 1)[1]] += 1
    median = statistics.median(c for c in counts.values() if c > 0)
    print(f"median rd:\t{median}", file=sys.stderr)
    with open(args.output, "w") as out:
        out.write("#chr\tpos\trd\n")
        for name, count in counts.items():
            chrom, pos = name.split("_")[:2]
            out.write(f"{chrom}\t{pos}\t{count / median!r}\n")


if __name__ == "__main__":
    main()
