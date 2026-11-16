import numpy as np

from binomial import binomial_pvalues, binomial_test


def test_balanced_counts_give_one():
    assert binomial_test(5, 10) == 1.0
    assert binomial_test(60, 120) == 1.0


def test_exact_tail():
    assert np.isclose(binomial_test(0, 10), 2 / 1024, rtol=1e-12)


def test_normal_approximation_above_100_trials():
    p = binomial_test(40, 101)
    assert 0 < p < 0.05
    assert np.isclose(p, binomial_test(61, 101), rtol=1e-12)


def test_vectorized_matches_scalar():
    k = np.array([0, 3, 7, 50, 80, 101])
    n = np.array([10, 10, 10, 100, 150, 200])
    expected = [binomial_test(a, b) for a, b in zip(k, n)]
    assert binomial_pvalues(k, n).tolist() == expected


def test_float_successes_are_truncated():
    assert binomial_test(3.0, 10) == binomial_test(3, 10)
