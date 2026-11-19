import gzip
import io
import random
import sys

from genes_bed import main as genes_bed_main
from hetsnv_bed import het_snv
from orient_hets import orient
from read_depth import main as read_depth_main
from rename_contigs import rename_chain, rename_fasta
from wgs_windows import main as wgs_windows_main

REFERENCE = [f"chr{c}" for c in [*range(1, 23), "X", "Y", "M"]]


def record(chrom, pos, ref, alt, gt, status="PASS"):
    return [chrom, str(pos), ".", ref, alt, "50", status, ".", "GT:DP", f"{gt}:10"]


def test_het_snv_keeps_pass_biallelic_snv_hets():
    assert het_snv(record("chr1", 100, "A", "G", "0/1"), 9, REFERENCE) == [
        "chr1",
        "99",
        "100",
        "chr1_100_A_A_G",
    ]
    assert (
        het_snv(record("chr1", 100, "A", "G", "1|0"), 9, REFERENCE)[3]
        == "chr1_100_A_G_A"
    )
    assert (
        het_snv(record("chr1", 100, "A", "G,T", "1/2"), 9, REFERENCE)[3]
        == "chr1_100_A_G_T"
    )


def test_het_snv_rejects_the_rest():
    assert het_snv(record("chr1", 100, "A", "G", "1/1"), 9, REFERENCE) is None
    assert (
        het_snv(record("chr1", 100, "A", "G", "0/1", "LowQual"), 9, REFERENCE) is None
    )
    assert het_snv(record("chr1", 100, "AT", "A", "0/1"), 9, REFERENCE) is None
    assert het_snv(record("chr1", 100, "A", "G", "./."), 9, REFERENCE) is None
    assert het_snv(record("chrUn_x", 100, "A", "G", "0/1"), 9, REFERENCE) is None


def test_het_snv_renames_to_reference_convention():
    assert het_snv(record("1", 5, "C", "T", "0/1"), 9, REFERENCE)[0] == "chr1"
    assert het_snv(record("MT", 5, "C", "T", "0/1"), 9, REFERENCE)[0] == "chrM"


def with_phase_set(chrom, pos, gt, phase_set):
    return [
        chrom,
        str(pos),
        ".",
        "A",
        "G",
        "50",
        "PASS",
        ".",
        "GT:PS",
        f"{gt}:{phase_set}",
    ]


def test_orient_keeps_phased_and_swaps_unphased_at_random():
    rng = random.Random(0)
    phased = record("chr1", 1, "A", "G", "1|0")
    assert orient(phased, rng, {}) == ("1|0", "1|0", True)
    reference = random.Random(0)
    unphased = record("chr1", 1, "A", "G", "0/1")
    expected = "1|0" if reference.random() < 0.5 else "0|1"
    assert orient(unphased, random.Random(0), {}) == ("0/1", expected, False)
    assert unphased[9] == f"{expected}:10"


def test_orient_leaves_homozygous_unswapped():
    rng = random.Random(1)
    assert orient(record("chr1", 1, "A", "G", "1/1"), rng, {})[1] == "1|1"
    assert rng.random() == random.Random(1).random()


def test_orient_treats_a_missing_phase_set_as_chromosome_wide():
    assert orient(with_phase_set("chr1", 1, "1|0", "."), random.Random(0), {}) == (
        "1|0",
        "1|0",
        True,
    )


def test_orient_flips_each_phase_set_as_one_block():
    rng, flips = random.Random(3), {}
    block = [
        with_phase_set("chr1", p, gt, 100) for p, gt in ((100, "0|1"), (150, "1|0"))
    ]
    other = with_phase_set("chr1", 900, "0|1", 900)
    first, second, third = (orient(r, rng, flips) for r in (*block, other))
    assert first[1] != second[1]
    assert not (first[2] or second[2] or third[2])
    reference = random.Random(3)
    flip_block, flip_other = reference.random() < 0.5, reference.random() < 0.5
    assert first[1] == ("1|0" if flip_block else "0|1")
    assert third[1] == ("1|0" if flip_other else "0|1")
    assert rng.random() == reference.random()


def test_rename_fasta_drops_contigs_and_descriptions():
    lines = [">chr1 desc\n", "ACGT\n", ">chrM\n", "AC\n", ">chr2\n", "GG\n"]
    assert list(rename_fasta(lines, "_hap1", {"chrM"})) == [
        ">chr1_hap1\n",
        "ACGT\n",
        ">chr2_hap1\n",
        "GG\n",
    ]


def test_rename_chain_changes_only_the_query():
    lines = ["chain 19 chr1 20 + 0 20 chr1 19 + 0 19 1\n", "15 1 0\n"]
    assert (
        next(rename_chain(lines, "_hap2", set()))
        == "chain 19 chr1 20 + 0 20 chr1_hap2 19 + 0 19 1\n"
    )


def test_wgs_windows_order_by_header_and_drop_negative_starts(write, tmp_path):
    hets = write(
        "hets.bed",
        "chr2\t4999\t5000\tchr2_5000_A_A_G\nchr1\t499\t500\tchr1_500_A_A_G\nchr1\t1999\t2000\tchr1_2000_C_C_T\n",
    )
    header = write("h.sam", "@HD\tVN:1.6\n@SQ\tSN:1\tLN:9000\n@SQ\tSN:2\tLN:9000\n")
    windows, genome = tmp_path / "w.bed", tmp_path / "g.txt"
    wgs_windows_main(
        [
            "--hets",
            hets,
            "--header",
            header,
            "--windows",
            str(windows),
            "--genome",
            str(genome),
        ]
    )
    assert windows.read_text().splitlines() == [
        "1\t1000\t3000\tchr1_2000_C_C_T",
        "2\t4000\t6000\tchr2_5000_A_A_G",
    ]
    assert genome.read_text() == "1\t9000\n2\t9000\n"


def test_read_depth_counts_overlaps_per_window(tmp_path, monkeypatch):
    names = [f"chr1_{p}_A_A_G" for p in (10, 20, 30, 40)]
    windows = tmp_path / "w.bed"
    windows.write_text("".join(f"chr1\t0\t1\t{n}\n" for n in names))
    overlaps = "".join(
        f"chr1\t5\t105\tread{i}\t255\t+\tchr1\t0\t1\t{n}\n"
        for n, c in zip(names, (10, 20, 0, 40), strict=True)
        for i in range(c)
    )
    monkeypatch.setattr(sys, "stdin", io.StringIO(overlaps))
    out = tmp_path / "rd.tsv"
    read_depth_main(["--windows", str(windows), "--output", str(out)])
    assert out.read_text().splitlines() == [
        "#chr\tpos\trd",
        "chr1\t10\t0.5",
        "chr1\t20\t1.0",
        "chr1\t30\t0.0",
        "chr1\t40\t2.0",
    ]


def test_genes_bed_uses_unique_gene_ids(tmp_path):
    gtf = tmp_path / "a.gtf.gz"
    with gzip.open(gtf, "wt") as out:
        out.write(
            'chr1\tx\tgene\t11\t20\t.\t+\t.\tgene_id "G1.1"; gene_type "lncRNA"; gene_name "DUP";\n'
        )
        out.write('chr1\tx\texon\t11\t20\t.\t+\t.\tgene_id "G1.1"; gene_name "DUP";\n')
        out.write(
            'chr2\tx\tgene\t5\t9\t.\t-\t.\tgene_id "G2.1"; gene_type "protein_coding"; gene_name "DUP";\n'
        )
    bed, names = tmp_path / "g.bed", tmp_path / "g.tsv"
    genes_bed_main(["--gtf", str(gtf), "--bed", str(bed), "--names", str(names)])
    assert bed.read_text().splitlines() == ["chr1\t10\t20\tG1.1", "chr2\t4\t9\tG2.1"]
    assert names.read_text().splitlines()[1:] == [
        "G1.1\tDUP\tlncRNA",
        "G2.1\tDUP\tprotein_coding",
    ]
