"""Check the integration run against the planted truth.

Run after ``pixi run test``. Skipped when the integration results are absent.
"""

import gzip
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


def calls(path):
    """Return each record's last-column FORMAT values, keyed by site and alleles."""
    with gzip.open(path, "rt") as vcf:
        rows = [line.rstrip("\n").split("\t") for line in vcf if line[0] != "#"]
    return {
        (r[0], r[1], r[3], r[4]): dict(zip(r[8].split(":"), r[-1].split(":")))
        for r in rows
    }


def test_hets_take_confident_phase_from_the_phased_vcf():
    joint = calls(RUN / "resources/joint.vcf.gz")
    phased = calls(RUN / "resources/phased.vcf.gz")
    oriented = calls(RESULTS / "refs/IND1/oriented.vcf.gz")
    confident = {
        key: call["GT"]
        for key, call in phased.items()
        if call["PP"] == "." or float(call["PP"]) >= 0.99
    }
    transferred = [k for k in confident if k in joint and joint[k]["GT"] == "0/1"]
    assert transferred
    assert all(oriented[k]["GT"] == confident[k] for k in transferred)
    assert set(oriented) == set(joint)
    already_phased = sum(call["GT"] in ("0|1", "1|0") for call in joint.values())
    hets, n_phased, _ = (
        (RESULTS / "refs/IND1/phasing.tsv").read_text().splitlines()[1].split("\t")
    )
    assert int(n_phased) == already_phased + len(transferred)
    assert int(hets) > int(n_phased)


@pytest.mark.parametrize("sample", ["RNA1", "RNA2"])
def test_star_temp_dir_sits_beside_its_outputs(sample):
    # STAR execs its readFilesCommand script from the temp dir, and SLURM job
    # /tmp is mounted noexec on zlurm.
    log = (RESULTS / f"{sample}/{sample}.Log.out").read_text()
    assert "--outTmpDir" not in log


def test_personal_genome_keeps_bgzipped_haplotypes_and_drops_the_index():
    refs = RESULTS / "refs/IND1"
    assert not (refs / "star_index").exists()
    for hap in (1, 2):
        fasta = refs / f"hap{hap}.fa.gz"
        assert fasta.read_bytes()[:4] == b"\x1f\x8b\x08\x04"
        assert (refs / f"hap{hap}.fa.gz.fai").exists()
        assert (refs / f"hap{hap}.fa.gz.gzi").exists()
