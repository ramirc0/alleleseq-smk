"""Check the integration run against the planted truth.

Run after ``pixi run test``. Skipped when the integration results are absent.
"""

from pathlib import Path

import pytest

RUN = Path(".tests/integration")
RESULTS = RUN / "results/default"
UNITS = ["RNA1", "RNA2", "POOL"]

pytestmark = pytest.mark.skipif(
    not RESULTS.exists(), reason="run `pixi run test` first"
)


def column(path, index=0):
    lines = path.read_text().splitlines()[1:]
    return {line.split("\t")[index] for line in lines}


@pytest.mark.parametrize("unit", UNITS)
@pytest.mark.parametrize("test", ["binomial", "betabinomial"])
def test_exactly_the_planted_genes_are_called(unit, test):
    planted = column(RUN / "resources/planted.tsv")
    called = column(RESULTS / f"{unit}/genes/{unit}.{test}.significant.tsv")
    assert called == planted


@pytest.mark.parametrize("unit", UNITS)
def test_cnv_filter_removes_the_doubled_stretch(unit):
    rows = [
        line.split("\t")
        for line in (RESULTS / f"{unit}/snvs/{unit}.discarded.cnv.log")
        .read_text()
        .splitlines()[1:]
    ]
    depth = [
        (chrom, int(pos)) for chrom, pos, cell in rows if not cell.endswith("__NA")
    ]
    missing = [int(pos) for _, pos, cell in rows if cell.endswith("__NA")]
    assert depth
    assert all(chrom == "chr2" and 30_000 <= pos <= 40_100 for chrom, pos in depth)
    # Windows that would start before the contig are dropped, as upstream did.
    assert all(pos < 1000 for pos in missing)


def test_named_tables_carry_gene_names():
    named = (
        (RESULTS / "POOL/genes/POOL.betabinomial.significant.named.tsv")
        .read_text()
        .splitlines()
    )
    assert named[0].split("\t")[:3] == ["#region", "gene_name", "gene_type"]
    assert all(line.split("\t")[1].startswith("TOY") for line in named[1:])


def test_phasing_summary():
    hets, phased, _ = (
        (RESULTS.parent / "refs/IND1/snps-indels/phasing.tsv")
        .read_text()
        .splitlines()[1]
        .split("\t")
    )
    assert int(hets) > int(phased) > 0
