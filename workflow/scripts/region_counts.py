"""Count distinct reads per haplotype over each region's het SNVs.

Port of AlleleSeq2's ``intersect2counts.py``. The input is the output of
``bedtools intersect -a <region SNVs> -b <BAMs> -split -wb -bed`` on stdin. A read
counts once per region and haplotype, however many of the region's SNVs it
covers. Mate suffixes ``/1`` and ``/2`` are dropped, so a pair counts once.
"""

import argparse
import sys
from collections import Counter

from binomial import binomial_test


def build_parser():
    """Build the region_counts.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the region_counts.py flags.
    """
    p = argparse.ArgumentParser(description="Count region reads per haplotype.")
    p.add_argument(
        "--kind", choices=["uniq", "mmap"], required=True, help="Read class."
    )
    p.add_argument(
        "--hap1-bed", required=True, help="Region SNVs on hap1, for SNV counts."
    )
    p.add_argument("--output", required=True, help="Region counts TSV to write.")
    return p


def main(argv=None):
    """Write the region counts table.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    regions = collect_reads(sys.stdin)
    with open(args.hap1_bed) as bed:
        snv_counts = Counter(line.rstrip("\n").split("\t")[3] for line in bed)
    with open(args.output, "w") as out:
        if args.kind == "mmap":
            out.write("#region\thap1_count\thap2_count\n")
            out.writelines(
                f"{name}\t{len(region['hap1'])}\t{len(region['hap2'])}\n"
                for name, region in regions.items()
            )
            return
        out.write(
            "#region\thap1_count\thap2_count\thap1_allele_ratio\tp_binom\t"
            "snv_count\tsnv_hap1_hap2_coords\n"
        )
        for name, region in regions.items():
            hap1, hap2 = len(region["hap1"]), len(region["hap2"])
            coords = ":" + "".join(f"{snv}:" for snv in region["snvs"])
            fields = [
                name,
                str(hap1),
                str(hap2),
                repr(hap1 / (hap1 + hap2)),
                repr(binomial_test(hap1, hap1 + hap2)),
                str(snv_counts[name]),
                coords,
            ]
            out.write("\t".join(fields) + "\n")


def collect_reads(lines):
    """Group read names by region and haplotype.

    Parameters
    ----------
    lines : iterable of str
        ``bedtools intersect -wb -bed`` lines. Ten columns for one BAM, eleven
        when a file index column precedes the read.

    Returns
    -------
    dict of str to dict
        Per region, in order of first appearance: read name sets under
        ``hap1`` and ``hap2``, and the covered SNV positions under ``snvs``.
    """
    regions = {}
    for line in lines:
        fields = line.rstrip("\n").split("\t")
        if len(fields) == 10:
            read = fields[7]
        elif len(fields) == 11:
            read = fields[8]
        else:
            raise ValueError(f"unexpected intersectBed line:\n{line}")
        read = read.split("/")[0]
        name = fields[3]
        hap = fields[0].rsplit("_", 1)[1]
        region = regions.setdefault(name, {"hap1": set(), "hap2": set(), "snvs": {}})
        region[hap].add(read)
        region["snvs"].setdefault(f"{fields[0]}_{fields[2]}")
    return regions


if __name__ == "__main__":
    main()
