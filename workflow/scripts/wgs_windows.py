"""Make 2 kb WGS windows around het SNVs in the alignment's contig order.

Port of the window step of AlleleSeq2's ``makePersonalGenome.mk`` and of
``feed_sorted_check_names.py``. Each window is ``[pos - 1000, pos + 1000)``
around the 1-based SNV position. Sites with a negative start are dropped, as
upstream's ``grep -v "-"`` did.
"""

import argparse
from collections import defaultdict

from chromosomes import harmonize, primary_names

HALF_WIDTH = 1000


def build_parser():
    """Build the wgs_windows.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the wgs_windows.py flags.
    """
    p = argparse.ArgumentParser(description="Windows around het SNVs for read depth.")
    p.add_argument(
        "--hets", required=True, help="Het SNV BED in reference coordinates."
    )
    p.add_argument("--header", required=True, help="SAM header of the WGS alignment.")
    p.add_argument("--windows", required=True, help="Window BED to write.")
    p.add_argument("--genome", required=True, help="bedtools genome file to write.")
    return p


def main(argv=None):
    """Write the genome file and the sorted windows.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    contigs = []
    with open(args.header) as header, open(args.genome, "w") as genome:
        for line in header:
            if line.startswith("@SQ"):
                tags = dict(
                    field.split(":", 1) for field in line.rstrip("\n").split("\t")[1:]
                )
                contigs.append(tags["SN"])
                genome.write(f"{tags['SN']}\t{tags['LN']}\n")
    primary = primary_names(contigs)
    windows = defaultdict(list)
    with open(args.hets) as hets:
        for line in hets:
            chrom, _, pos, name = line.rstrip("\n").split("\t")
            start = int(pos) - HALF_WIDTH
            if start < 0:
                continue
            contig = harmonize(chrom, primary)
            if contig is None:
                raise SystemExit(f"{chrom} has no match in the WGS alignment header")
            windows[contig].append(
                f"{contig}\t{start}\t{int(pos) + HALF_WIDTH}\t{name}\n"
            )
    with open(args.windows, "w") as out:
        for contig in primary:
            out.writelines(windows[contig])


if __name__ == "__main__":
    main()
