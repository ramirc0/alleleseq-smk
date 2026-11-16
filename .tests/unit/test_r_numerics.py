import numpy as np
import pytest
from rpy2 import robjects

from r_numerics import (
    HISTOGRAM_EDGES,
    r_hist_counts,
    r_prefix_sums,
    r_seq_by,
    r_signif,
    r_sum,
)


def r(expr):
    return np.asarray(robjects.r(expr), dtype=np.float64)


def test_edges_equal_r_pretty():
    assert HISTOGRAM_EDGES.tolist() == r("pretty(0:1, 40)").tolist()


@pytest.mark.parametrize(
    "start,stop,by",
    [
        (0, 0.01, 0.001),
        (0.01, 0.1, 0.01),
        (0.1, 1, 0.1),
        (0, 0.99, 0.1),
        (0, 1, 0.04),
        (0.0005, 0.0525, 0.0125),
        (0.0125, 0.0375, 0.003125),
        (5e-4, 5e-4, 2.5e-5),
    ],
)
def test_seq_by_equals_r(start, stop, by):
    assert (
        r_seq_by(start, stop, by).tolist()
        == r(f"seq({start!r}, {stop!r}, by={by!r})").tolist()
    )


def test_seq_wrong_sign():
    with pytest.raises(ValueError):
        r_seq_by(5e-4, 5e-5, 2.5e-5)


@pytest.mark.parametrize(
    "x",
    [
        0.000632054121995964,
        0.0278757249376457,
        0.00345385508339163,
        -0.0231808032331031,
        0.0006325,
        1.2345,
        0.1235,
        2.5e-7,
        123456.0,
    ],
)
def test_signif_equals_r(x):
    assert r_signif(x, 3) == r(f"signif({x!r}, 3)")[0]


def test_sum_equals_r():
    rng = np.random.default_rng(0)
    values = rng.random(10001) * 1e-3 + 1e3
    robjects.globalenv["v"] = robjects.FloatVector(values)
    assert r_sum(values) == r("sum(v)")[0]
    assert (
        r_prefix_sums(values)[[0, 99, 10000]].tolist()
        == r("c(sum(v[1]), sum(v[1:100]), sum(v))").tolist()
    )


def test_hist_counts_equal_r():
    values = np.array([0.0, 0.02, 0.025, 3 / 40, 0.5, 37 / 40, 1.0, 0.48, 1 / 3])
    robjects.globalenv["v"] = robjects.FloatVector(values)
    expected = r("hist(v, breaks=pretty(0:1, 40), right=TRUE, plot=FALSE)$counts")
    assert r_hist_counts(values).tolist() == expected.tolist()
