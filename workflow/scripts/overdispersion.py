"""Estimate beta-binomial overdispersion from the allele ratio histogram.

Port of AlleleDB's ``alleledb_calcOverdispersion.R`` (Chen et al. 2016) as
AlleleSeq2 runs it. The empirical histogram of allele ratios is compared to
null histograms built from binomial or beta-binomial masses, weighted by how
many rows have each total. A coarse grid over the overdispersion ``b`` stops
at the first increase in weighted squared error. Finer grids around the best
``b`` repeat until the error no longer changes in three significant digits.
The control flow follows the R script line by line, including its quirks.
"""

import argparse
import sys
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import poikilos as pk

from count_tables import read_count_table, totals_and_ratios
from r_numerics import (
    HISTOGRAM_EDGES,
    dbetabinom,
    dbinom,
    r_hist_counts,
    r_seq_by,
    r_signif,
    r_sum,
)

MAX_TOTAL = 2500
SEARCH_ROWS = 50


def build_parser():
    """Build the overdispersion.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the overdispersion.py flags.
    """
    p = argparse.ArgumentParser(description="Fit the beta-binomial overdispersion.")
    p.add_argument("--counts", required=True, help="Testable SNV or region counts TSV.")
    p.add_argument(
        "--min-total", type=int, default=6, help="Smallest total in the fit."
    )
    p.add_argument("--summary", required=True, help="TSV of the chosen b to write.")
    p.add_argument("--search", required=True, help="TSV of every b tried to write.")
    p.add_argument(
        "--plot", required=True, help="Fit plot SVG to write (PNG beside it)."
    )
    return p


def main(argv=None):
    """Fit b and write the summary, the search table and the plot.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    totals, ratios = totals_and_ratios(read_count_table(args.counts))
    if totals.size == 0:
        fit = Fit(0.0, float("nan"), np.zeros((1, 2)), 0, float("nan"), None, None)
    else:
        fit = fit_overdispersion(totals, ratios, args.min_total)
    with open(args.summary, "w") as out:
        out.write("b\tsse\tmax_total\tfraction_within_max_total\n")
        out.write(
            f"{fit.b:.15g}\t{fit.sse:.15g}\t{fit.max_total}\t{fit.fraction:.15g}\n"
        )
    with open(args.search, "w") as out:
        out.write("b\tsse\n")
        out.writelines(f"{b:.15g}\t{sse:.15g}\n" for b, sse in fit.search)
    plot_fit(fit, args.plot)


@dataclass
class Fit:
    """Result of the overdispersion search.

    Attributes
    ----------
    b : float
        Chosen overdispersion.
    sse : float
        Weighted squared error at the last accepted step.
    search : numpy.ndarray of float
        Rows of ``(b, sse)`` as the R script keeps them.
    max_total : int
        Largest total used in the fit.
    fraction : float
        Fraction of rows with a total of at most 2500.
    empirical : numpy.ndarray of float or None
        Empirical ratio density per bin.
    nulls : dict of str to numpy.ndarray or None
        Binomial and chosen beta-binomial null densities per bin.
    """

    b: float
    sse: float
    search: np.ndarray
    max_total: int
    fraction: float
    empirical: np.ndarray | None
    nulls: dict | None


def null_density(min_total, max_total, weights, rho=None):
    """Bin the weighted null masses of every total into the ratio histogram.

    Parameters
    ----------
    min_total, max_total : int
        Range of totals, inclusive.
    weights : numpy.ndarray of int
        Rows per total, indexed by total.
    rho : float, optional
        Beta-binomial overdispersion. None uses the binomial.

    Returns
    -------
    numpy.ndarray of float
        Density per histogram bin.
    """
    sizes = np.arange(min_total, max_total + 1)
    size = np.repeat(sizes, sizes + 1)
    successes = np.concatenate([np.arange(n + 1) for n in sizes])
    ratio = successes / size
    if rho is None:
        mass = dbinom(successes, size, 0.5)
    else:
        mass = dbetabinom(successes, size, 0.5, rho)
    weighted = mass * weights[size]
    order = np.lexsort((weighted, ratio))
    ratio, weighted = ratio[order], weighted[order]
    edges = HISTOGRAM_EDGES
    binned = np.empty(len(edges) - 1)
    for i in range(len(edges) - 1):
        side = "left" if i == 0 else "right"
        lo = np.searchsorted(ratio, edges[i], side)
        hi = np.searchsorted(ratio, edges[i + 1], "right")
        binned[i] = r_sum(weighted[lo:hi])
    return binned / r_sum(binned)


def gradient_weights():
    """Return the bin weights of the squared error, low at the edges.

    Returns
    -------
    numpy.ndarray of float
        One weight per histogram bin.
    """
    n_bins = len(HISTOGRAM_EDGES) - 1
    steps = r_seq_by(0, 1, 1 / (n_bins / 2))[1:]
    tail = steps if n_bins % 2 == 0 else steps[:-1]
    return np.concatenate((steps, np.sort(tail)[::-1]))


def fit_overdispersion(totals, ratios, min_total):
    """Search the overdispersion as ``alleledb_calcOverdispersion.R`` does.

    Parameters
    ----------
    totals : numpy.ndarray of int
        Reads on the two alleles per row.
    ratios : numpy.ndarray of float
        Tested allele ratio per row.
    min_total : int
        Smallest total in the fit.

    Returns
    -------
    Fit
        The chosen b and the search trace.
    """
    max_total = int(min(totals.max(), MAX_TOTAL))
    fraction = np.count_nonzero(totals <= MAX_TOTAL) / totals.size
    w_grad = gradient_weights()
    match = (totals <= max_total) & (totals >= min_total)
    counts = r_hist_counts(ratios[match])
    empirical = counts / counts.sum()
    weights = np.bincount(totals)

    binomial = null_density(min_total, max_total, weights)
    sse = r_sum((empirical - binomial) ** 2)
    search = np.zeros((SEARCH_ROWS, 2))
    b_choice = 0.0
    ctr = 1
    by = 0.1
    for k in r_seq_by(0, 0.99, by):
        density = null_density(min_total, max_total, weights, k)
        if b_choice == 0:
            search[0] = (b_choice, sse)
        sse_bbin = r_sum(w_grad * (empirical - density) ** 2)
        _store(search, ctr, k, sse_bbin)
        if sse_bbin < sse:
            sse, b_choice = sse_bbin, k
        elif sse_bbin > sse:
            break
        ctr += 1

    searching = b_choice < 0.9
    while searching:
        start = max(0, b_choice - by / 2)
        stop = b_choice + by / 2
        by = by / 4
        sse = search[0, 1]
        b_choice = 0.0
        for k in r_seq_by(start, stop, by):
            density = null_density(min_total, max_total, weights, k)
            sse_bbin = r_sum(w_grad * (empirical - density) ** 2)
            _store(search, ctr + 1, k, sse_bbin)
            if sse_bbin < sse:
                sse, b_choice = sse_bbin, k
            elif sse_bbin > sse:
                break
            ctr += 1
        print(f"b.chosen= {b_choice:.15g} , SSE.chosen= {sse:.15g}", file=sys.stderr)
        if ctr + 1 >= SEARCH_ROWS:
            raise IndexError(
                "subscript out of bounds in the b search (R's 50-row matrix)"
            )
        if r_signif(search[ctr + 1, 1], 3) == r_signif(search[ctr, 1], 3):
            searching = False

    nulls = {
        "binomial": binomial,
        "beta-binomial": null_density(min_total, max_total, weights, b_choice),
    }
    return Fit(b_choice, sse, search[: ctr + 2], max_total, fraction, empirical, nulls)


def _store(search, row, b, sse):
    """Write one search row, failing where R's matrix would overflow.

    Parameters
    ----------
    search : numpy.ndarray of float
        Search matrix.
    row : int
        Zero-based row.
    b, sse : float
        Values to store.
    """
    if row >= SEARCH_ROWS:
        raise IndexError("subscript out of bounds in the b search (R's 50-row matrix)")
    search[row] = (b, sse)


def plot_fit(fit, path):
    """Plot the empirical histogram against the binomial and chosen nulls.

    Parameters
    ----------
    fit : Fit
        Search result.
    path : str
        SVG path. The PNG goes beside it.
    """
    pk.use("plain")
    fig, (left, right) = plt.subplots(1, 2, figsize=(9, 3.5))
    if fit.empirical is not None:
        mids = (HISTOGRAM_EDGES[:-1] + HISTOGRAM_EDGES[1:]) / 2
        left.bar(mids, fit.empirical, width=0.018, color="0.8", label="empirical")
        left.plot(mids, fit.nulls["binomial"], marker="o", ms=3, label="binomial")
        left.plot(
            mids,
            fit.nulls["beta-binomial"],
            marker="o",
            ms=3,
            label=f"beta-binomial, b={fit.b:.3g}",
        )
        left.legend(frameon=False, fontsize=7)
        trace = fit.search[np.argsort(fit.search[:, 0], kind="stable")]
        right.plot(trace[:, 0], trace[:, 1], marker="o", ms=3)
        right.plot([fit.b], [fit.sse], marker="*", ms=10, color="C3")
    left.set(
        xlabel="allele ratio", ylabel="density", title=f"totals up to {fit.max_total}"
    )
    right.set(xlabel="b", ylabel="weighted SSE", title="overdispersion search")
    pk.despine(left)
    pk.despine(right)
    pk.save_figure(fig, path)
    plt.close(fig)


if __name__ == "__main__":
    main()
