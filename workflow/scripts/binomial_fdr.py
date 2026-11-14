"""Pick a binomial p-value cutoff by simulation and keep the rows below it.

Port of AlleleSeq2's ``FalsePos.py`` and ``filter_by_pval.py`` (and their
region twins). Each simulation flips a fair coin once per read at every
testable site and records the binomial p-values. At a cutoff ``pv``, the FDR
is the mean number of simulated p-values at or below ``pv`` divided by one
plus the index of ``pv`` among the sorted observed p-values. The walk starts
at the largest observed p-value and stops at the first one whose FDR is
below the target.

The simulation reproduces upstream's random stream exactly: Python 2's
``random.randint(0, 1)`` after ``random.seed(0)`` equals ``floor(2 * u)``
for successive ``random.random()`` draws, and NumPy's legacy MT19937 yields
the same draws once it holds the same state.
"""

import argparse
import random
from bisect import bisect_right

import numpy as np

from binomial import NORMAL_APPROXIMATION_MIN_TRIALS, binomial_pvalues
from count_tables import observed_pvalues, read_count_table, two_largest

GRID = (
    [e * 0.001 for e in range(10)]
    + [e * 0.01 for e in range(1, 10)]
    + [e * 0.1 for e in range(1, 10)]
)
FLIPS_PER_DRAW = 1 << 24
SNV_COLUMNS = [0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14]
REGION_COLUMNS = [0, 1, 2, 3, 4, 5]


def build_parser():
    """Build the binomial_fdr.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the binomial_fdr.py flags.
    """
    p = argparse.ArgumentParser(description="Binomial FDR cutoff by simulation.")
    p.add_argument("--counts", required=True, help="Testable SNV or region counts TSV.")
    p.add_argument("--report", required=True, help="FDR report to write.")
    p.add_argument("--output", required=True, help="Significant rows TSV to write.")
    p.add_argument("--fdr", type=float, default=0.05, help="Target FDR.")
    p.add_argument(
        "--simulations", type=int, default=500, help="Number of simulations."
    )
    return p


def main(argv=None):
    """Write the FDR report and the significant rows.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    table = read_count_table(args.counts)
    higher, lower = two_largest(table)
    pvalues = observed_pvalues(table)
    report, cutoff = fdr_report(higher + lower, pvalues, args.simulations, args.fdr)
    with open(args.report, "w") as out:
        out.write(f"# binomial_fdr.py {args.counts} {args.simulations} {args.fdr}\n")
        out.write("\n".join(report) + "\n")
    columns = REGION_COLUMNS if table.is_region else SNV_COLUMNS
    with open(args.output, "w") as out:
        out.write("\t".join(table.header[i] for i in columns) + "\n")
        for row, pvalue in zip(table.rows, pvalues):
            if pvalue <= cutoff:
                out.write("\t".join(row[i] for i in columns) + "\n")


def upstream_random_state():
    """Return a NumPy RandomState holding Python's state after ``seed(0)``.

    Returns
    -------
    numpy.random.RandomState
        Generator whose ``random_sample`` repeats ``random.random()``.
    """
    state = random.Random(0).getstate()[1]
    generator = np.random.RandomState()
    generator.set_state(("MT19937", np.array(state[:-1], dtype=np.uint32), state[-1]))
    return generator


def simulated_heads(trials, generator):
    """Flip a fair coin `trials[i]` times for each site, in site order.

    Parameters
    ----------
    trials : numpy.ndarray of int
        Flips per site.
    generator : numpy.random.RandomState
        Source of uniform draws, consumed in order.

    Returns
    -------
    numpy.ndarray of int
        Heads per site.
    """
    heads = np.empty(len(trials), dtype=np.int64)
    ends = np.cumsum(trials)
    start_site = 0
    while start_site < len(trials):
        first_flip = ends[start_site] - trials[start_site]
        stop_site = max(
            int(np.searchsorted(ends, first_flip + FLIPS_PER_DRAW, "right")),
            start_site + 1,
        )
        n_flips = int(ends[stop_site - 1] - first_flip)
        flips = (generator.random_sample(n_flips) * 2).astype(np.int64)
        running = np.concatenate(([0], np.cumsum(flips)))
        block_ends = ends[start_site:stop_site] - first_flip
        heads[start_site:stop_site] = (
            running[block_ends] - running[block_ends - trials[start_site:stop_site]]
        )
        start_site = stop_site
    return heads


def pvalue_lookup(trials):
    """Precompute exact-branch p-values for every (heads, trials) pair.

    Parameters
    ----------
    trials : numpy.ndarray of int
        Reads per site.

    Returns
    -------
    numpy.ndarray of float
        Square table indexed by ``[trials, heads]`` for trials below
        ``NORMAL_APPROXIMATION_MIN_TRIALS``.
    """
    size = NORMAL_APPROXIMATION_MIN_TRIALS
    n, k = np.meshgrid(np.arange(size), np.arange(size), indexing="ij")
    table = np.zeros((size, size))
    valid = k <= n
    table[valid] = binomial_pvalues(k[valid], n[valid])
    return table


def fdr_report(trials, observed, simulations, target):
    """Simulate the null and walk to the p-value cutoff.

    Parameters
    ----------
    trials : numpy.ndarray of int
        Reads per site, in file order.
    observed : numpy.ndarray of float
        Observed p-values, in file order.
    simulations : int
        Number of simulations.
    target : float
        Target FDR.

    Returns
    -------
    lines : list of str
        Report lines, as ``FalsePos.py`` prints them.
    cutoff : float
        P-value cutoff as printed with six decimals, or 0 when none exists.
    """
    act = np.sort(observed)
    grid = np.array(GRID)
    at_act = np.zeros(len(act), dtype=np.int64)
    at_grid = np.zeros(len(grid), dtype=np.int64)
    generator = upstream_random_state()
    table = pvalue_lookup(trials)
    exact = trials < NORMAL_APPROXIMATION_MIN_TRIALS
    for _ in range(simulations):
        heads = simulated_heads(trials, generator)
        simulated = np.empty(len(trials))
        simulated[exact] = table[trials[exact], heads[exact]]
        simulated[~exact] = binomial_pvalues(heads[~exact], trials[~exact])
        at_act += np.cumsum(
            np.bincount(np.searchsorted(act, simulated, "left"), minlength=len(act) + 1)
        )[:-1]
        at_grid += np.cumsum(
            np.bincount(
                np.searchsorted(grid, simulated, "left"), minlength=len(grid) + 1
            )
        )[:-1]

    lines = ["pval\tP\tFP\tFDR"]
    for pv, n_sim in zip(GRID, at_grid):
        n_act = bisect_right(act, pv)
        mean_sim = n_sim / simulations
        lines.append(f"{pv:f}\t{n_act}\t{mean_sim:f}\t{mean_sim / (n_act + 1):f}")

    last_fdr = last_pv = 0.0
    for index in range(len(act) - 1, -1, -1):
        pv = act[index]
        fdr = at_act[index] / simulations / (index + 1)
        if fdr < target:
            lines.append(f"target {target:f}")
            lines.append(f"before {last_fdr:f} {last_pv:f}")
            lines.append(f"after  {fdr:f} {pv:f}")
            lines.append(f"Target {target:f} FDR {fdr:f} pv {pv:f}")
            return lines, float(f"{pv:f}")
        last_fdr, last_pv = fdr, pv
    lines.append(
        "# couldn't calculate pvalue threshold: no snvs in input or none with small "
        "enough p-value (?)\n# will set p value threshold to 0: \n"
        f"Target {target} FDR NA pv 0"
    )
    return lines, 0.0


if __name__ == "__main__":
    main()
