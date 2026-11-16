from filter_counts import main

HEADER = (
    "#chr\tref_coord\thap1_coord\thap2_coord\tref_allele\thap1_allele\thap2_allele\t"
    "cA\tcC\tcG\tcT\tcN\tref_allele_ratio\tsum_ref_n_alt_cnts\tp_binom\twarn_hap1\twarn_hap2"
)


def row(pos, a, g, warn1=".", warn2="."):
    return "\t".join(
        [
            "chr1",
            str(pos),
            f"chr1_hap1_{pos}",
            f"chr1_hap2_{pos}",
            "A",
            "A",
            "G",
            str(a),
            "0",
            str(g),
            "0",
            "0",
            "x",
            str(a + g),
            "x",
            warn1,
            warn2,
        ]
    )


def run(write, tmp_path, rows, depth, mmap=""):
    counts = write("raw.tsv", HEADER + "\n" + "\n".join(rows) + "\n")
    rd = write("rd.tab", "#chr\tpos\trd\n" + depth)
    mm = write("mm.mpileup", mmap)
    paths = {name: tmp_path / f"{name}.txt" for name in ("out", "cnv", "warn", "mmap")}
    main(
        [
            "--counts",
            counts,
            "--read-depth",
            rd,
            "--mmap-pileups",
            mm,
            "--output",
            str(paths["out"]),
            "--cnv-log",
            str(paths["cnv"]),
            "--warning-log",
            str(paths["warn"]),
            "--mmap-log",
            str(paths["mmap"]),
        ]
    )
    return {name: path.read_text().splitlines() for name, path in paths.items()}


def test_depth_filter_and_missing_depth(write, tmp_path):
    out = run(
        write,
        tmp_path,
        [row(10, 5, 5), row(20, 5, 5), row(30, 5, 5)],
        "chr1\t10\t1.0\nchr1\t20\t2.0\n",
    )
    assert [line.split("\t")[1] for line in out["out"][1:]] == ["10"]
    assert out["cnv"][1:] == ["chr1\t20\t5_0_5_0__2.0", "chr1\t30\t5_0_5_0__NA"]


def test_warnings_are_dropped(write, tmp_path):
    out = run(
        write, tmp_path, [row(10, 5, 5, warn1="1_other_alleles")], "chr1\t10\t1.0\n"
    )
    assert len(out["out"]) == 1
    assert out["warn"][1:] == ["chr1\t10\t5_0_5_0__1_other_alleles__."]


def test_weaker_allele_gets_multimapping_reads_up_to_the_stronger(write, tmp_path):
    out = run(
        write,
        tmp_path,
        [row(10, 8, 2)],
        "chr1\t10\t1.0\n",
        "chr1_hap2\t10\tG\t9\t.........\tIIIIIIIII\n",
    )
    fields = out["out"][1].split("\t")
    assert fields[7:12] == ["8", "0", "8", "0", "0"]
    assert fields[12:14] == ["0.5", "16"]
    assert fields[-2:] == ["1.0", "G:+6"]
