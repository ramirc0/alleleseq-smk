import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    "haplotypes", Path(__file__).parents[1] / "golden/haplotypes.py"
)
haplotypes = importlib.util.module_from_spec(spec)
spec.loader.exec_module(haplotypes)


def test_golden_orientation_clears_stale_phase_sets(write, tmp_path):
    vcf = write(
        "in.vcf",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        "chr1\t10\t.\tA\tG\t50\tPASS\t.\tGT:PS\t0|1:10\n"
        "chr1\t20\t.\tC\tT\t50\tPASS\t.\tGT:DP\t0/1:9\n",
    )
    out = tmp_path / "out.vcf"
    haplotypes.write_phased_vcf(vcf, out, {("chr1", 10): "G", ("chr1", 20): "C"})
    assert [line.split("\t")[9] for line in out.read_text().splitlines()[1:]] == [
        "1|0:.",
        "0|1:9",
    ]
