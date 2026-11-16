from filter_region_counts import adjust
from region_counts import collect_reads
from region_hets import main as region_hets_main


def test_snv_in_overlapping_regions_counts_for_each(write, tmp_path):
    counts = write(
        "c.tsv",
        "#chr\tref_coord\thap1_coord\thap2_coord\n"
        "chr1\t50\tchr1_hap1_50\tchr1_hap2_51\n"
        "chr1\t500\tchr1_hap1_500\tchr1_hap2_501\n",
    )
    # Upstream kept only GENE_B for position 50.
    bed = write(
        "g.bed", "chr1\t10\t100\tGENE_A\nchr1\t40\t60\tGENE_B\nchr1\t499\t500\tGENE_C\n"
    )
    hap1, hap2 = tmp_path / "h1.bed", tmp_path / "h2.bed"
    region_hets_main(
        ["--counts", counts, "--regions", bed, "--hap1", str(hap1), "--hap2", str(hap2)]
    )
    assert hap1.read_text().splitlines() == [
        "chr1_hap1\t49\t50\tGENE_A",
        "chr1_hap1\t49\t50\tGENE_B",
        "chr1_hap1\t499\t500\tGENE_C",
    ]
    assert hap2.read_text().splitlines()[0] == "chr1_hap2\t50\t51\tGENE_A"


def intersect_line(contig, end, region, read):
    return f"{contig}\t{end - 1}\t{end}\t{region}\t{contig}\t0\t100\t{read}\t255\t+\n"


def test_reads_count_once_per_region_and_haplotype():
    lines = [
        intersect_line("chr1_hap1", 10, "G", "r1/1"),
        intersect_line("chr1_hap1", 20, "G", "r1/2"),
        intersect_line("chr1_hap1", 20, "G", "r2"),
        intersect_line("chr1_hap2", 10, "G", "r3"),
    ]
    region = collect_reads(lines)["G"]
    assert (len(region["hap1"]), len(region["hap2"])) == (2, 1)


def test_snv_coordinates_are_not_lost_to_substring_matches():
    # Upstream skipped chr1_hap1_12 because ":chr1_hap1_123:" contains it.
    lines = [
        intersect_line("chr1_hap1", 123, "G", "r1"),
        intersect_line("chr1_hap1", 12, "G", "r2"),
    ]
    assert list(collect_reads(lines)["G"]["snvs"]) == ["chr1_hap1_123", "chr1_hap1_12"]


class Log(list):
    def write(self, text):
        self.append(text)


def test_region_weaker_haplotype_gets_multimapping_reads():
    row = ["G", "10", "2", "0.83", "x", "3", ":c:"]
    out = adjust(row, {"G": ("1", "20")}, Log())
    assert out[:3] == ["G", "10", "10"]
    assert out[-1] == "hap2:+8"


def test_region_without_multimapping_reads():
    log = Log()
    assert (
        adjust(["G", "3", "1", "0.75", "x", "1", ":c:"], {}, log)[-1] == "no_mmap_reads"
    )
    assert log == ["G\t0_0\tno_mmap_reads\n"]
