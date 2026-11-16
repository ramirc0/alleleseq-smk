import random

import numpy as np

from binomial import binomial_test
from binomial_fdr import fdr_report, simulated_heads, upstream_random_state


def test_random_state_repeats_python_random():
    reference = random.Random(0)
    generator = upstream_random_state()
    assert generator.random_sample(1000).tolist() == [
        reference.random() for _ in range(1000)
    ]


def test_heads_follow_upstream_coin_flips():
    trials = np.array([3, 0, 7, 1, 20])
    reference = random.Random(0)
    expected = [sum(int(reference.random() * 2) for _ in range(n)) for n in trials]
    assert simulated_heads(trials, upstream_random_state()).tolist() == expected


def test_heads_across_draw_blocks(monkeypatch):
    monkeypatch.setattr("binomial_fdr.FLIPS_PER_DRAW", 4)
    trials = np.array([3, 0, 7, 1, 20])
    reference = random.Random(0)
    expected = [sum(int(reference.random() * 2) for _ in range(n)) for n in trials]
    assert simulated_heads(trials, upstream_random_state()).tolist() == expected


def falsepos_reference(trials, observed, sims, target):
    """Line-for-line Python 3 transcription of FalsePos.py's arithmetic."""
    import bisect

    rng = random.Random(0)
    act = sorted(observed)
    sim_pvals = []
    for _ in range(sims):
        row = []
        for n in trials:
            a = sum(int(rng.random() * 2) for _ in range(n))
            row.append(binomial_test(a, n))
        sim_pvals.append(sorted(row))
    for nact, pv in sorted(enumerate(act), reverse=True):
        mean = np.mean([bisect.bisect(sim_pvals[i], pv) for i in range(sims)])
        fdr = mean / (nact + 1)
        if fdr < target:
            return f"Target {target:f} FDR {fdr:f} pv {pv:f}"
    return None


def test_fdr_report_matches_reference_walk():
    rng = np.random.default_rng(1)
    trials = rng.integers(6, 140, size=60)
    heads = np.array(
        [rng.binomial(n, 0.8 if i < 15 else 0.5) for i, n in enumerate(trials)]
    )
    observed = np.array([binomial_test(h, n) for h, n in zip(heads, trials)])
    lines, cutoff = fdr_report(trials, observed, 20, 0.2)
    assert lines[-1] == falsepos_reference(trials, observed, 20, 0.2)
    assert cutoff == float(lines[-1].split()[-1])


def test_fdr_report_without_cutoff():
    lines, cutoff = fdr_report(np.array([10]), np.array([1.0]), 5, 0.05)
    assert cutoff == 0.0
    assert lines[-1].endswith("Target 0.05 FDR NA pv 0")
