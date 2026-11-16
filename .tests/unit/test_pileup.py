import pytest

from pileup import parse_line, pileup_counts, read_bases


def test_read_bases_strips_starts_ends_and_indels():
    assert read_bases("^~.,$A+2CG-10NNNNNNNNNNt*") == ".,At*"


def test_read_start_mapq_character_is_not_a_base():
    # MAPQ 32 prints as "A" after "^". Upstream's parser counted it as a base.
    assert read_bases("^A.") == "."


def test_parse_line_counts_and_warnings():
    key, site = parse_line("chr1_hap1\t100\tc\t6\t..,,Ag\tIIIIII")
    assert key == "chr1_hap1_100"
    assert site.allele == "C"
    assert site.counts == {"A": 1, "C": 4, "G": 1, "T": 0, "N": 0}
    assert site.warning == "2_other_alleles"


def test_hap_allele_not_largest():
    _, site = parse_line("chr1_hap1\t100\tC\t3\t.TT\tIII")
    assert site.warning == "hap_allele_not_largest_cnt"


def test_depth_mismatch_raises():
    with pytest.raises(ValueError):
        parse_line("chr1_hap1\t100\tC\t5\t..\tII")


def test_deletions_and_skips_count_toward_depth():
    _, site = parse_line("chr1_hap1\t100\tC\t4\t.*><\tIIII")
    assert site.counts["C"] == 1


def test_pooled_files_sum_every_base_including_n(write):
    a = write("a.mpileup", "chr1_hap1\t5\tA\t2\t.N\tII\n")
    b = write("b.mpileup", "chr1_hap1\t5\tA\t3\t.nN\tIII\n")
    site = pileup_counts([a, b])["chr1_hap1_5"]
    assert site.counts == {"A": 2, "C": 0, "G": 0, "T": 0, "N": 3}
    assert site.warning == "1_other_alleles,hap_allele_not_largest_cnt"
