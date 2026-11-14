"""Add weaker-haplotype multi-mapping reads to region counts.

Port of AlleleSeq2's ``filter_regions_w_mmaps.py adjust``. The weaker
haplotype gains its multi-mapping reads, up to the stronger haplotype's
count, and the binomial test is redone.
"""

import argparse

from binomial import binomial_test


def build_parser():
    """Build the filter_region_counts.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the filter_region_counts.py flags.
    """
    p = argparse.ArgumentParser(
        description="Adjust region counts for multi-mapping reads."
    )
    p.add_argument("--counts", required=True, help="Unique-read region counts TSV.")
    p.add_argument("--mmap", required=True, help="Multi-mapping region counts TSV.")
    p.add_argument(
        "--output", required=True, help="Filtered region counts TSV to write."
    )
    p.add_argument("--log", required=True, help="Log of adjustments to write.")
    return p


def main(argv=None):
    """Write the filtered region counts and the log.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with open(args.mmap) as mmap_file:
        mmap_file.readline()
        mmap = {}
        for line in mmap_file:
            region, hap1, hap2 = line.rstrip("\n").split("\t")
            mmap[region] = (hap1, hap2)
    with (
        open(args.counts) as counts,
        open(args.output, "w") as out,
        open(args.log, "w") as log,
    ):
        log.write("#region\tmm_count:hap1_hap2\tmm_log\n")
        out.write(counts.readline().rstrip("\n") + "\tmmap_log\n")
        for line in counts:
            row = line.rstrip("\n").split("\t")
            out.write("\t".join(adjust(row, mmap, log)) + "\n")


def adjust(row, mmap, log):
    """Adjust one region row.

    Parameters
    ----------
    row : list of str
        Unique-read region row.
    mmap : dict of str to tuple of str
        Multi-mapping hap1 and hap2 counts by region.
    log : file
        Log of the multi-mapping counts and the action taken.

    Returns
    -------
    list of str
        The row with a trailing ``mmap_log`` column.
    """
    region = row[0]
    if region not in mmap:
        log.write(f"{region}\t0_0\tno_mmap_reads\n")
        return row + ["no_mmap_reads"]
    mm1, mm2 = mmap[region]
    hap1, hap2 = int(row[1]), int(row[2])
    if hap1 == hap2:
        log.write(f"{region}\t{mm1}_{mm2}\teq_uniq_cnts\n")
        return row + ["."]
    weaker = "hap2" if hap1 > hap2 else "hap1"
    weaker_mm = int(mm2 if weaker == "hap2" else mm1)
    if weaker_mm == 0:
        log.write(f"{region}\t{mm1}_{mm2}\t{weaker}:0\n")
        return row + [f"{weaker}:0"]
    weaker_unique = min(hap1, hap2)
    adjusted = min(weaker_unique + weaker_mm, max(hap1, hap2))
    added = adjusted - weaker_unique
    if weaker == "hap1":
        hap1 = adjusted
    else:
        hap2 = adjusted
    log.write(f"{region}\t{mm1}_{mm2}\t{weaker}:+{added}\n")
    return [
        region,
        str(hap1),
        str(hap2),
        repr(hap1 / (hap1 + hap2)),
        repr(binomial_test(hap1, hap1 + hap2)),
        row[5],
        row[6],
        f"{weaker}:+{added}",
    ]


if __name__ == "__main__":
    main()
