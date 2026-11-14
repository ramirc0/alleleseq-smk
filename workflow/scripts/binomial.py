"""Two-sided binomial test at p = 0.5, as AlleleSeq2 computes it.

AlleleSeq2 (from AlleleSeq v1.2a) uses a continuity-corrected normal
approximation when ``n * 0.5 > 50`` and SciPy's exact ``binom_test``
otherwise. The exact branch is written out with ``bdtr`` and ``bdtrc`` the way
SciPy 1.2 evaluated it, so ties between p-values behave as upstream.
"""

import numpy as np
from scipy.special import bdtr, bdtrc, ndtr

NORMAL_APPROXIMATION_MIN_TRIALS = 101


def binomial_pvalues(successes, trials):
    """Return two-sided binomial p-values for a fair coin.

    Parameters
    ----------
    successes : array_like of int
        Successes per test. Floats are truncated, as SciPy 1.2 did.
    trials : array_like of int
        Trials per test, broadcast against `successes`.

    Returns
    -------
    numpy.ndarray of float
        One p-value per test, capped at 1.
    """
    k = np.asarray(successes).astype(np.int64)
    n = np.asarray(trials).astype(np.int64)
    k, n = np.broadcast_arrays(k, n)
    pvalues = np.empty(k.shape, dtype=np.float64)
    normal = n >= NORMAL_APPROXIMATION_MIN_TRIALS
    pvalues[normal] = _normal_approximation(k[normal], n[normal])
    pvalues[~normal] = _exact(k[~normal], n[~normal])
    return pvalues


def binomial_test(successes, trials):
    """Return the two-sided binomial p-value of one test.

    Parameters
    ----------
    successes : int
        Successes.
    trials : int
        Trials.

    Returns
    -------
    float
        The p-value, capped at 1.
    """
    return float(binomial_pvalues([successes], [trials])[0])


def _exact(k, n):
    """Evaluate SciPy 1.2's two-sided ``binom_test`` at p = 0.5.

    Parameters
    ----------
    k : numpy.ndarray of int
        Successes.
    n : numpy.ndarray of int
        Trials.

    Returns
    -------
    numpy.ndarray of float
        P-values.
    """
    expected = n * 0.5
    below = k < expected
    above = k > expected
    pvalues = np.ones(k.shape, dtype=np.float64)
    kb, nb = k[below], n[below]
    pvalues[below] = bdtr(kb, nb, 0.5) + bdtrc(nb - kb - 1, nb, 0.5)
    ka, na = k[above], n[above]
    pvalues[above] = bdtr(na - ka, na, 0.5) + bdtrc(ka - 1, na, 0.5)
    return np.minimum(1.0, pvalues)


def _normal_approximation(k, n):
    """Evaluate AlleleSeq's continuity-corrected normal approximation.

    Parameters
    ----------
    k : numpy.ndarray of int
        Successes.
    n : numpy.ndarray of int
        Trials.

    Returns
    -------
    numpy.ndarray of float
        P-values.
    """
    mean = n * 0.5
    sd = np.sqrt(n * 0.5 * 0.5)
    lower = 2 * ndtr((k + 0.5 - mean) / sd)
    upper = 2 * (1 - ndtr((k - 0.5 - mean) / sd))
    pvalues = np.where(k < mean, lower, upper)
    return np.where(np.abs(k - mean) < 1e-5, 1.0, pvalues)
