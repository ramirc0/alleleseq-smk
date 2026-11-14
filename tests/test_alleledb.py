import numpy as np
from rpy2 import robjects

from betabinomial_fdr import NullPvalues, bisect, lower_tail, thresholds
from overdispersion import gradient_weights


def r(expr):
    return np.asarray(robjects.r(expr), dtype=np.float64)


def test_thresholds_equal_r_grid():
    expected = r(
        "c(seq(0,0.01,by=0.001), seq(0.01,0.1,by=0.01)[-1], seq(0.1,1,by=0.1)[-1])"
    )
    assert thresholds().tolist() == expected.tolist()


def test_gradient_weights_equal_r():
    robjects.r("library(VGAM)")
    expected = r(
        "{bins=pretty(0:1,40); r=seq(0,1,(1-0)/((length(bins)-1)/2)); r=r[2:length(r)];"
        " c(r, sort(r[1:length(r)], decreasing=TRUE))}"
    )
    assert gradient_weights().tolist() == expected.tolist()


def test_lower_tail_equals_vgam_pbetabinom():
    robjects.r("library(VGAM)")
    expected = r("pbetabinom(0:40, 40, 0.5, 0.0199)")
    assert lower_tail(40, 0.0199).tolist() == expected.tolist()
    assert lower_tail(40, None).tolist() == r("pbinom(0:40, 40, 0.5)").tolist()


def test_bisection_evaluates_start_on_a_degenerate_interval():
    # Starting at 0, start becomes 5e-4 > stop = 5e-5. Unpatched R stopped here
    # with "wrong sign in 'by'". The patch evaluates start alone.
    nulls = NullPvalues(np.array([10]), np.array([1]), [np.full(11, 1e-3)])
    observed = np.array([1e-4, 1e-4, 0.5])
    assert bisect(observed, nulls, 0.0, 0.0, 0.05) == 0.0005125
