from count_alleles import main


def run(write, tmp_path, hap1_pileup, hap2_pileup, extra_bed=""):
    ref = write(
        "ref.bed", "chr1\t99\t100\tchr1_100_A_A_G\n" + extra_bed.replace("_hapX", "")
    )
    hap1 = write(
        "hap1.bed",
        "chr1_hap1\t99\t100\tchr1_100_A_A_G\n" + extra_bed.replace("_hapX", "_hap1"),
    )
    hap2 = write(
        "hap2.bed",
        "chr1_hap2\t99\t100\tchr1_100_A_A_G\n" + extra_bed.replace("_hapX", "_hap2"),
    )
    p1 = write("h1.mpileup", hap1_pileup)
    p2 = write("h2.mpileup", hap2_pileup)
    out, log = tmp_path / "counts.tsv", tmp_path / "discarded.log"
    main(
        [
            "--ref-bed",
            ref,
            "--hap1-bed",
            hap1,
            "--hap2-bed",
            hap2,
            "--discarded",
            str(log),
            "--output",
            str(out),
            "--pileups",
            p1,
            p2,
        ]
    )
    rows = [line.split("\t") for line in out.read_text().splitlines()]
    return rows, log.read_text()


def test_counts_both_haplotypes(write, tmp_path):
    rows, log = run(
        write,
        tmp_path,
        "chr1_hap1\t100\tA\t3\t...\tIII\n",
        "chr1_hap2\t100\tG\t1\t.\tI\n",
    )
    assert rows[1][:15] == [
        "chr1",
        "100",
        "chr1_hap1_100",
        "chr1_hap2_100",
        "A",
        "A",
        "G",
        "3",
        "0",
        "1",
        "0",
        "0",
        "0.75",
        "4",
        "0.625",
    ]
    assert log == ""


def test_haplotype_bases_decide_allele_order(write, tmp_path):
    rows, _ = run(
        write,
        tmp_path,
        "chr1_hap1\t100\tG\t2\t..\tII\n",
        "chr1_hap2\t100\tA\t1\t.\tI\n",
    )
    assert rows[1][5:7] == ["G", "A"]


def test_unexpected_haplotype_base_is_discarded(write, tmp_path):
    rows, log = run(
        write,
        tmp_path,
        "chr1_hap1\t100\tT\t2\t..\tII\n",
        "chr1_hap2\t100\tA\t1\t.\tI\n",
    )
    assert len(rows) == 1
    assert log == "chr1\t100\tvcf_vs_mpileup_alleles\n"


def test_shared_position_removes_both(write, tmp_path):
    rows, log = run(
        write,
        tmp_path,
        "chr1_hap1\t100\tA\t1\t.\tI\n",
        "",
        "chr1_hapX\t99\t100\tchr1_100_A_A_C\n",
    )
    assert len(rows) == 1
    assert "duplicate_pos" in log
