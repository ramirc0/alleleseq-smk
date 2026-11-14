"""R arithmetic that the AlleleDB beta-binomial code depends on.

AlleleDB's overdispersion fit and FDR search stop on exact comparisons, so
the port reproduces R's numerics rather than approximating them. The
distribution functions are R's own, called through rpy2. ``seq``,
``signif``, ``sum`` and ``hist`` are re-implemented and unit-tested against R.
"""

import math

import numpy as np
from rpy2.robjects import default_converter, numpy2ri
from rpy2.robjects.packages import importr

# R's pretty(0:1, 40): 51 edges, 0.02 apart.
HISTOGRAM_EDGES = np.arange(51) * 0.02

_CONVERTER = default_converter + numpy2ri.converter
_STATS = importr("stats")
_VGAM = importr("VGAM")

if np.finfo(np.longdouble).nmant != 63:
    raise RuntimeError(
        "R's sum() accumulates in 80-bit long double; this platform lacks it."
    )


def _call(function, *args):
    """Call an R function on NumPy arrays and return a float64 array.

    Parameters
    ----------
    function : rpy2 function
        The R function.
    *args : array_like or float
        Arguments, converted to R vectors.

    Returns
    -------
    numpy.ndarray of float
        The R result.
    """
    with _CONVERTER.context():
        converted = [np.atleast_1d(np.asarray(a, dtype=np.float64)) for a in args]
        return np.asarray(function(*converted), dtype=np.float64)


def dbinom(x, size, prob):
    """Return R's ``dbinom(x, size, prob)``.

    Parameters
    ----------
    x : array_like of int
        Successes.
    size : array_like of int
        Trials, recycled against `x`.
    prob : float
        Success probability.

    Returns
    -------
    numpy.ndarray of float
        Probability masses.
    """
    return _call(_STATS.dbinom, x, size, prob)


def pbinom(q, size, prob):
    """Return R's ``pbinom(q, size, prob)``.

    Parameters
    ----------
    q : array_like of int
        Quantiles.
    size : array_like of int
        Trials, recycled against `q`.
    prob : float
        Success probability.

    Returns
    -------
    numpy.ndarray of float
        Lower-tail probabilities.
    """
    return _call(_STATS.pbinom, q, size, prob)


def dbetabinom(x, size, prob, rho):
    """Return VGAM's ``dbetabinom(x, size, prob, rho)``.

    Parameters
    ----------
    x : array_like of int
        Successes.
    size : array_like of int
        Trials, recycled against `x`.
    prob : float
        Mean success probability.
    rho : float
        Correlation (overdispersion) parameter. 0 gives the binomial.

    Returns
    -------
    numpy.ndarray of float
        Probability masses.
    """
    return _call(_VGAM.dbetabinom, x, size, prob, rho)


def r_sum(values):
    """Sum like R's ``sum()``: in order, accumulated in long double.

    Parameters
    ----------
    values : array_like of float
        Values to add.

    Returns
    -------
    float
        The sum rounded to double.
    """
    values = np.asarray(values, dtype=np.float64)
    if values.size == 0:
        return 0.0
    return float(np.cumsum(values, dtype=np.longdouble)[-1])


def r_prefix_sums(values):
    """Return R's ``sum(values[1:i])`` for every i at once.

    Parameters
    ----------
    values : array_like of float
        Values to add.

    Returns
    -------
    numpy.ndarray of float
        Running sums, each rounded to double.
    """
    return np.cumsum(np.asarray(values, dtype=np.float64), dtype=np.longdouble).astype(
        np.float64
    )


def r_seq_by(start, stop, by):
    """Return R's ``seq(start, stop, by = by)`` for doubles.

    Parameters
    ----------
    start, stop : float
        Range ends.
    by : float
        Step.

    Returns
    -------
    numpy.ndarray of float
        ``start + (0:n) * by``, capped at `stop`.

    Raises
    ------
    ValueError
        If `by` points away from `stop`, as R's "wrong sign in 'by'".
    """
    delta = stop - start
    if delta == 0 and stop == 0:
        return np.array([stop])
    n = delta / by
    if n < 0:
        raise ValueError("wrong sign in 'by' argument")
    if abs(delta) / max(abs(stop), abs(start)) < 100 * np.finfo(np.float64).eps:
        return np.array([start])
    steps = int(n + 1e-10)
    values = start + np.arange(steps + 1) * by
    return np.minimum(values, stop) if by > 0 else np.maximum(values, stop)


def r_signif(x, digits):
    """Round to significant digits like R's ``signif()`` (``fprec``).

    Parameters
    ----------
    x : float
        Value to round.
    digits : int
        Significant digits.

    Returns
    -------
    float
        The rounded value.
    """
    if math.isnan(x) or not math.isfinite(x) or x == 0:
        return x
    sign = -1.0 if x < 0 else 1.0
    x = abs(x)
    l10 = math.log10(x)
    e10 = int(digits - 1 - math.floor(l10))
    max10e = 308
    if abs(l10) >= max10e - 2:
        raise ValueError("signif() of extreme magnitudes is not ported")
    p10 = 1.0
    if e10 > max10e:
        p10 = 10.0 ** (e10 - max10e)
        e10 = max10e
    if e10 > 0:
        pow10 = 10.0**e10
        return sign * (round((x * pow10) * p10) / pow10) / p10
    pow10 = 10.0**-e10
    return sign * (round(x / pow10) * pow10)


def r_hist_counts(values, edges=HISTOGRAM_EDGES):
    """Count values per bin like ``hist(values, breaks=edges, right=TRUE)``.

    Parameters
    ----------
    values : array_like of float
        Values inside the edge range.
    edges : numpy.ndarray of float
        Bin edges.

    Returns
    -------
    numpy.ndarray of int
        Counts per right-closed bin, the first bin closed on both sides.
    """
    diddle = 1e-7 * np.median(np.diff(edges))
    fuzzy = edges + np.concatenate(([-diddle], np.full(len(edges) - 1, diddle)))
    bins = np.searchsorted(fuzzy, np.asarray(values, dtype=np.float64), "left") - 1
    return np.bincount(bins, minlength=len(edges) - 1)
