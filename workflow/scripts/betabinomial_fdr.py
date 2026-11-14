"""Test each row against a beta-binomial null and pick an FDR cutoff.

Port of AlleleDB's ``alleledb_alleleseqBetabinomial.R`` as AlleleSeq2 runs
it, including the user's patch for degenerate bisection intervals. Each row's
p-value is twice the lower tail of the null at its second-largest allele
count, capped at 1. For a cutoff, the expected false positives sum, over every
total, the largest null p-value at or below the cutoff times the number of
rows with that total. A grid picks a starting cutoff and a bisection refines
it. The control flow follows the R script line by line.
"""

import argparse
import sys
from dataclasses import dataclass

import numpy as np

from count_tables import read_count_table, two_largest
from r_numerics import dbetabinom, pbinom, r_prefix_sums, r_seq_by, r_signif, r_sum

STEP = 0.0001
PSEUDO_ZERO_RHO = 1e-10
BISECTION_ROWS = 100
RD_MIN, RD_MAX = 0.5, 1.5


def build_parser():
    """Build the betabinomial_fdr.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the betabinomial_fdr.py flags.
    """
    p = argparse.ArgumentParser(description="Beta-binomial test with an FDR cutoff.")
    p.add_argument("--counts", required=True, help="Testable SNV or region counts TSV.")
    p.add_argument(
        "--overdispersion", required=True, help="Summary from overdispersion.py."
    )
    p.add_argument("--fdr", type=float, default=0.05, help="Target FDR.")
    p.add_argument(
        "--min-total", type=int, default=6, help="Smallest total in the FDR."
    )
    p.add_argument("--all", required=True, help="Counts with p.betabin to write.")
    p.add_argument("--output", required=True, help="Significant rows TSV to write.")
    p.add_argument("--report", required=True, help="FDR report to write.")
    return p


def main(argv=None):
    """Write the tested table, the significant rows and the FDR report.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    table = read_count_table(args.counts)
    with open(args.overdispersion) as handle:
        handle.readline()
        b = float(handle.readline().split("\t")[0])
    higher, lower = two_largest(table)
    if table.rows:
        result = beta_binomial_fdr(higher + lower, lower, b, args.fdr, args.min_total)
    else:
        result = FdrResult(np.array([]), 0.0, ["# no testable rows"])
    keep = result.p_betabin <= result.cutoff_betabin
    if "cnv" in table.header:
        cnv = np.array(table.column("cnv"), dtype=np.float64)
        keep &= (cnv >= RD_MIN) & (cnv <= RD_MAX)
    header = "\t".join([*table.header, "p.betabin"]) + "\n"
    with open(args.all, "w") as tested, open(args.output, "w") as significant:
        tested.write(header)
        significant.write(header)
        for row, pvalue, kept in zip(table.rows, result.p_betabin, keep):
            line = "\t".join([*row, repr(float(pvalue))]) + "\n"
            tested.write(line)
            if kept:
                significant.write(line)
    with open(args.report, "w") as out:
        out.write("\n".join(result.report) + "\n")


@dataclass
class FdrResult:
    """Outcome of the beta-binomial FDR search.

    Attributes
    ----------
    p_betabin : numpy.ndarray of float
        Beta-binomial p-value per row.
    cutoff_betabin : float
        Chosen beta-binomial p-value cutoff.
    report : list of str
        Report lines in the layout of AlleleDB's FDR file.
    """

    p_betabin: np.ndarray
    cutoff_betabin: float
    report: list


@dataclass
class NullPvalues:
    """Null p-values of every count, for every observed total.

    Attributes
    ----------
    totals : numpy.ndarray of int
        Distinct totals, ascending.
    rows_per_total : numpy.ndarray of int
        Rows with each total.
    pvalues : list of numpy.ndarray
        For each total ``t``, the p-values of counts ``0..t``.
    """

    totals: np.ndarray
    rows_per_total: np.ndarray
    pvalues: list

    def pvalue(self, total, count):
        """Return the null p-value of one observed count.

        Parameters
        ----------
        total : int
            Trials.
        count : int
            Second-largest allele count.

        Returns
        -------
        float
            Twice the lower tail, capped at 1.
        """
        index = np.searchsorted(self.totals, total)
        if index == len(self.totals) or self.totals[index] != total:
            raise ValueError(f"total {total} is below the minimum total of the FDR")
        return self.pvalues[index][count]

    def false_positives(self, cutoff):
        """Return AlleleDB's expected false positives at a cutoff.

        Parameters
        ----------
        cutoff : float
            P-value cutoff.

        Returns
        -------
        float
            Sum over totals of the rows times the largest null p-value at or
            below the cutoff.
        """
        largest = [
            np.max(np.where(p <= cutoff, p, 0.0) * n)
            for p, n in zip(self.pvalues, self.rows_per_total)
        ]
        return r_sum(largest)


def thresholds():
    """Return AlleleDB's grid of p-value cutoffs.

    Returns
    -------
    numpy.ndarray of float
        0 to 0.01 by 0.001, 0.02 to 0.1 by 0.01, then 0.2 to 1 by 0.1.
    """
    return np.concatenate(
        (
            r_seq_by(0, 0.01, 0.001),
            r_seq_by(0.01, 0.1, 0.01)[1:],
            r_seq_by(0.1, 1, 0.1)[1:],
        )
    )


def null_pvalues(totals, min_total, rho):
    """Tabulate twice the lower-tail null probability for each observed total.

    Parameters
    ----------
    totals : numpy.ndarray of int
        Total per row.
    min_total : int
        Smallest total kept.
    rho : float or None
        Beta-binomial overdispersion. None uses the binomial.

    Returns
    -------
    NullPvalues
        P-values of every count for every distinct total of at least
        `min_total`.
    """
    distinct, rows = np.unique(totals, return_counts=True)
    kept = distinct >= min_total
    distinct, rows = distinct[kept], rows[kept]
    pvalues = [np.minimum(2 * lower_tail(t, rho), 1.0) for t in distinct]
    return NullPvalues(distinct, rows, pvalues)


def lower_tail(total, rho):
    """Return the null CDF at counts ``0..total``.

    Parameters
    ----------
    total : int
        Trials.
    rho : float or None
        Beta-binomial overdispersion. None uses the binomial.

    Returns
    -------
    numpy.ndarray of float
        R's ``pbinom`` or VGAM's ``pbetabinom`` (a running sum of masses).
    """
    counts = np.arange(total + 1)
    if rho is None:
        return pbinom(counts, total, 0.5)
    return r_prefix_sums(dbetabinom(counts, total, 0.5, rho))


def beta_binomial_fdr(totals, lower, b, target, min_total):
    """Run AlleleDB's binomial and beta-binomial FDR searches.

    Parameters
    ----------
    totals : numpy.ndarray of int
        Sum of the two largest allele counts per row.
    lower : numpy.ndarray of int
        Second-largest allele count per row.
    b : float
        Overdispersion from the fit. 0 becomes 1e-10, as in AlleleDB.
    target : float
        Target FDR.
    min_total : int
        Smallest total in the false-positive sums.

    Returns
    -------
    FdrResult
        P-values, the beta-binomial cutoff and the report.
    """
    rho = PSEUDO_ZERO_RHO if b == 0 else b
    grid = thresholds()
    nulls = {
        "bin": null_pvalues(totals, min_total, None),
        "betabin": null_pvalues(totals, min_total, rho),
    }
    observed = {
        "bin": np.minimum(2 * pbinom(lower, totals, 0.5), 1.0),
        "betabin": np.array(
            [nulls["betabin"].pvalue(t, k) for t, k in zip(totals, lower)],
            dtype=np.float64,
        ),
    }
    columns, old, new = {}, {}, {}
    for kind in ("bin", "betabin"):
        fp = np.array([nulls[kind].false_positives(c) for c in grid])
        tp = np.array([np.count_nonzero(observed[kind] <= c) for c in grid]) + 1
        fdr = fp / tp
        passing = fdr <= target
        old[kind] = grid[passing].max()
        new[kind] = bisect(
            observed[kind], nulls[kind], old[kind], fdr[passing].max(), target
        )
        columns[kind] = (tp, fp, fdr)

    report = ["pval\tP.bin\tFP.bin\tFDR.bin\tP.betabin\tFP.betabin\tFDR.betabin"]
    for i, cutoff in enumerate(grid[:-1]):
        cells = [cutoff]
        for kind in ("bin", "betabin"):
            tp, fp, fdr = columns[kind]
            cells += [tp[i], fp[i], fdr[i] if np.isfinite(fdr[i]) else 0]
        report.append("\t".join(f"{c:.15g}" for c in cells))
    report += [
        f"FDR.threshold = {target:.15g}",
        f"p.choice.bin.old = {old['bin']:.15g}",
        f"p.choice.betabin.old = {old['betabin']:.15g}",
        f"p.choice.bin = {new['bin']:.15g}",
        f"p.choice.betabin = {new['betabin']:.15g}",
    ]
    return FdrResult(observed["betabin"], new["betabin"], report)


def bisect(observed, nulls, choice, fdr, target):
    """Refine the cutoff as AlleleDB's ``bisect()`` does, with the user's patch.

    Parameters
    ----------
    observed : numpy.ndarray of float
        Observed p-values.
    nulls : NullPvalues
        Null p-values for the false-positive sums.
    choice : float
        Starting cutoff from the grid.
    fdr : float
        FDR at the starting cutoff.
    target : float
        Target FDR.

    Returns
    -------
    float
        Cutoff of the last evaluated row whose FDR stayed below the target.
    """
    trace = np.zeros((BISECTION_ROWS, 3))
    trace[0] = (choice, fdr, 10)
    ctr = 0
    by = STEP
    searching = True
    while searching:
        start = max(0, choice - by / 2)
        stop = choice + by / 2
        by = by / 4
        if start == 0:
            start = 5e-4
        cutoffs = r_seq_by(start, stop, by) if start < stop and by > 0 else [start]
        for i in cutoffs:
            tp = np.count_nonzero(observed <= i)
            fp = nulls.false_positives(i)
            with np.errstate(divide="ignore", invalid="ignore"):
                fdr_i = np.float64(fp) / np.float64(tp)
            if np.isnan(fdr_i):
                raise ValueError("missing value where TRUE/FALSE needed (FP = TP = 0)")
            error = target - fdr_i
            ctr += 1
            if ctr >= BISECTION_ROWS:
                raise IndexError(
                    "subscript out of bounds in bisect (R's 100-row matrix)"
                )
            trace[ctr] = (i, fdr_i, error)
            choice = i
            print(
                f"i= {i:.15g} fdr.ind= {fdr_i:.15g} fp= {fp:.15g} tp= {tp}",
                file=sys.stderr,
            )
            if error < 0:
                break
        if r_signif(trace[ctr - 1, 2], 3) == r_signif(trace[ctr, 2], 3):
            searching = False
    return float(trace[trace[:, 2] > 0][-1, 0])


if __name__ == "__main__":
    main()
