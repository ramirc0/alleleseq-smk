"""Place each region's het SNVs on both haplotypes.

Port of AlleleSeq2's ``refbed2hapcoords.py`` together with the
``intersectBed`` call before it. Upstream kept only the last region of a SNV
that falls in overlapping regions. Here a SNV is placed in every region that
contains it.
"""

import argparse
import heapq
from collections import defaultdict


def build_parser():
    """Build the region_hets.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the region_hets.py flags.
    """
    p = argparse.ArgumentParser(description="Assign het SNVs to regions per haplotype.")
    p.add_argument("--counts", required=True, help="Filtered het SNV counts TSV.")
    p.add_argument("--regions", required=True, help="BED4 of regions.")
    p.add_argument("--hap1", required=True, help="BED of region SNVs on hap1 to write.")
    p.add_argument("--hap2", required=True, help="BED of region SNVs on hap2 to write.")
    return p


def main(argv=None):
    """Write the hap1 and hap2 region SNV BEDs.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with open(args.counts) as counts:
        counts.readline()
        rows = [line.rstrip("\n").split("\t") for line in counts]
    regions = regions_by_position(args.regions, [(row[0], int(row[1])) for row in rows])
    with open(args.hap1, "w") as hap1, open(args.hap2, "w") as hap2:
        for row in rows:
            chrom, pos = row[0], int(row[1])
            hap1_pos = int(row[2].rsplit("_", 1)[-1])
            hap2_pos = int(row[3].rsplit("_", 1)[-1])
            for name in regions.get((chrom, pos), []):
                hap1.write(f"{chrom}_hap1\t{hap1_pos - 1}\t{hap1_pos}\t{name}\n")
                hap2.write(f"{chrom}_hap2\t{hap2_pos - 1}\t{hap2_pos}\t{name}\n")


def regions_by_position(bed_path, positions):
    """Find the regions that contain each 1-based position.

    Parameters
    ----------
    bed_path : str
        BED4 of regions (0-based, half-open).
    positions : list of tuple of (str, int)
        Chromosome and 1-based position of each SNV.

    Returns
    -------
    dict of tuple to list of str
        Region names per position, in BED file order.
    """
    by_chrom = defaultdict(list)
    with open(bed_path) as bed:
        for order, line in enumerate(bed):
            chrom, start, end, name = line.rstrip("\n").split("\t")[:4]
            by_chrom[chrom].append((int(start), int(end), order, name))
    wanted = defaultdict(set)
    for chrom, pos in positions:
        wanted[chrom].add(pos)

    hits = {}
    for chrom, points in wanted.items():
        intervals = sorted(by_chrom.get(chrom, []))
        active = []
        next_interval = 0
        for pos in sorted(points):
            base = pos - 1
            while (
                next_interval < len(intervals) and intervals[next_interval][0] <= base
            ):
                start, end, order, name = intervals[next_interval]
                heapq.heappush(active, (end, order, name))
                next_interval += 1
            while active and active[0][0] <= base:
                heapq.heappop(active)
            if active:
                hits[(chrom, pos)] = [
                    name for _, _, name in sorted(active, key=lambda a: a[1])
                ]
    return hits


if __name__ == "__main__":
    main()
