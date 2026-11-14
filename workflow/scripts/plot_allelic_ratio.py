"""Plot the allele ratio distribution of a counts table.

Port of AlleleSeq2's ``plot_AllelicRatio_distribution.R``: a 41-bin
histogram of the reference (SNV) or hap1 (region) allele ratio.
"""

import argparse

import matplotlib.pyplot as plt
import poikilos as pk

from count_tables import read_count_table, totals_and_ratios


def build_parser():
    """Build the plot_allelic_ratio.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the plot_allelic_ratio.py flags.
    """
    p = argparse.ArgumentParser(description="Histogram of allele ratios.")
    p.add_argument("--counts", required=True, help="SNV or region counts TSV.")
    p.add_argument("--title", required=True, help="Plot title.")
    p.add_argument("--plot", required=True, help="SVG to write (PNG beside it).")
    return p


def main(argv=None):
    """Write the histogram.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    table = read_count_table(args.counts)
    _, ratios = totals_and_ratios(table)
    label = "hap1 allele ratio" if table.is_region else "reference allele ratio"
    pk.use("plain")
    fig, ax = plt.subplots(figsize=(4, 3.5))
    ax.hist(ratios, bins=41, range=(0, 1), color="0.4")
    ax.set(xlabel=label, ylabel="count", title=f"{args.title} (n = {len(ratios):,})")
    pk.despine(ax)
    pk.save_figure(fig, args.plot)
    plt.close(fig)


if __name__ == "__main__":
    main()
