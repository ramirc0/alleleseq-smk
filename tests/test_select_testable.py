from select_testable import main


def test_snv_thresholds_and_chromosomes(write, tmp_path):
    header = "#chr\tref_coord\th1\th2\tref\tA1\tA2\tcA\tcC\tcG\tcT\tcN\n"
    rows = [
        "chr1\t1\tx\tx\tA\tA\tG\t3\t0\t3\t0\t0",
        "chr1\t2\tx\tx\tA\tA\tG\t3\t0\t2\t0\t0",
        "chrX\t3\tx\tx\tA\tA\tG\t9\t0\t9\t0\t0",
    ]
    counts = write("c.tsv", header + "\n".join(rows) + "\n")
    out = tmp_path / "t.tsv"
    main(["--counts", counts, "--output", str(out)])
    assert [r.split("\t")[1] for r in out.read_text().splitlines()[1:]] == ["1"]
    main(
        [
            "--counts",
            counts,
            "--output",
            str(out),
            "--keep-chromosomes",
            "X",
            "--min-per-allele",
            "3",
        ]
    )
    assert [r.split("\t")[1] for r in out.read_text().splitlines()[1:]] == ["1", "3"]


def test_region_chromosome_from_bed(write, tmp_path):
    counts = write("r.tsv", "#region\thap1_count\thap2_count\nG1\t4\t4\nG2\t9\t9\n")
    bed = write("g.bed", "chr2\t0\t10\tG1\nchrY\t0\t10\tG2\n")
    out = tmp_path / "t.tsv"
    main(["--counts", counts, "--output", str(out), "--regions", bed])
    assert out.read_text().splitlines()[1:] == ["G1\t4\t4"]
