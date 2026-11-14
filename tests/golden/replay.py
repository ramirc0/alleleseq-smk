"""Replay the ported steps on the upstream AlleleSeq2 run and compare outputs.

The upstream run of IND1_RNA_rep1 (individual IND1) is controlled
data, so it is read in place and never copied into the repository. Each step
runs on upstream's own inputs, so a mismatch points at one step.

Usage::

    pixi run golden [<upstream run dir>] [<scratch dir>]

The upstream directory defaults to ``$ALLELESEQ_GOLDEN``. The scratch
directory defaults to a fresh temporary directory.
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

from compare import compare_tables

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow/scripts"
PREFIX = "IND1_RNA_rep1"
POOL = "IND1"
SAMPLE = "IND1"
SNV_TAG = "chrs1-22.6-tot_0-min"


def run(script, *args, stdin=None):
    """Run a ported script with the current interpreter.

    Parameters
    ----------
    script : str
        Script file name in ``workflow/scripts``.
    *args : str
        Arguments.
    stdin : file, optional
        Standard input.
    """
    subprocess.run(
        [sys.executable, str(SCRIPTS / script), *map(str, args)],
        stdin=stdin,
        check=True,
        stderr=subprocess.DEVNULL,
    )


def main(golden=None, scratch=None):
    """Replay every step and print one line per comparison.

    Parameters
    ----------
    golden : str, optional
        Upstream run directory.
    scratch : str, optional
        Directory for the ported outputs.
    """
    g = Path(golden or os.environ["ALLELESEQ_GOLDEN"])
    s = Path(scratch or tempfile.mkdtemp(prefix="alleleseq-golden-"))
    s.mkdir(parents=True, exist_ok=True)
    pg = g / "pgenome_IND1"
    results = []

    def check(name, ours, theirs, **kwargs):
        problems = compare_tables(str(ours), str(theirs), **kwargs)
        results.append(not problems)
        print(
            f"{'ok  ' if not problems else 'FAIL'} {name}"
            + (f": {problems[:3]}" if problems else ""),
            flush=True,
        )

    run(
        "hetsnv_bed.py",
        "--vcf",
        g / f"resources/{SAMPLE}.snvs.vcf",
        "--sample",
        SAMPLE,
        "--fai",
        g / "resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.fai",
        "--output",
        s / "ref.bed",
    )
    check("het SNV BED", s / "ref.bed", pg / f"{SAMPLE}_hetSNVs_ref.bed")

    run(
        "count_alleles.py",
        "--ref-bed",
        pg / f"{SAMPLE}_hetSNVs_ref.bed",
        "--hap1-bed",
        pg / f"{SAMPLE}_hetSNVs_hap1.bed",
        "--hap2-bed",
        pg / f"{SAMPLE}_hetSNVs_hap2.bed",
        "--discarded",
        s / "discarded.log",
        "--output",
        s / "raw.tsv",
        "--pileups",
        g / f"{PREFIX}_hap1_uniqreads.mpileup",
        g / f"{PREFIX}_hap2_uniqreads.mpileup",
    )
    check("raw counts", s / "raw.tsv", g / f"{PREFIX}_raw_counts.tsv")

    run(
        "filter_counts.py",
        "--counts",
        g / f"{PREFIX}_raw_counts.tsv",
        "--read-depth",
        pg / f"{SAMPLE}_hetSNVs_rd.tab",
        "--mmap-pileups",
        g / f"{PREFIX}_hap1_mmapreads.mpileup",
        g / f"{PREFIX}_hap2_mmapreads.mpileup",
        "--output",
        s / "filtered.tsv",
        "--cnv-log",
        s / "cnv.log",
        "--warning-log",
        s / "warn.log",
        "--mmap-log",
        s / "mmap.log",
    )
    check("filtered counts", s / "filtered.tsv", g / f"{PREFIX}_filtered_counts.tsv")
    check("CNV log", s / "cnv.log", g / f"{PREFIX}_discarded_HetSNVs_potential-CNV.log")
    check(
        "phase warning log",
        s / "warn.log",
        g / f"{PREFIX}_discarded_HetSNVs_warn-haplotype.log",
    )
    check(
        "multi-mapping log", s / "mmap.log", g / f"{PREFIX}_mmap_reads_over_hetSNVs.log"
    )

    run(
        "select_testable.py",
        "--counts",
        g / f"{PREFIX}_filtered_counts.tsv",
        "--output",
        s / "testable.tsv",
    )
    check(
        "testable SNVs",
        s / "testable.tsv",
        g / f"{PREFIX}_filtered_counts.{SNV_TAG}_cnt.tsv",
    )

    testable = g / f"{PREFIX}_filtered_counts.{SNV_TAG}_cnt.tsv"
    run(
        "binomial_fdr.py",
        "--counts",
        testable,
        "--report",
        s / "binom.txt",
        "--output",
        s / "binom.tsv",
    )
    check(
        "binomial FDR report",
        s / "binom.txt",
        g / f"{PREFIX}_FDR-0.05.binom.{SNV_TAG}_cnt.txt",
        rtol=0,
        skip_prefix="#",
    )
    check(
        "binomial significant SNVs",
        s / "binom.tsv",
        g / f"{PREFIX}_interestingHets.FDR-0.05.binom.{SNV_TAG}_cnt.tsv",
    )

    run(
        "overdispersion.py",
        "--counts",
        testable,
        "--summary",
        s / "od.tsv",
        "--search",
        s / "od_search.tsv",
        "--plot",
        s / "od.svg",
    )
    check(
        "overdispersion search",
        s / "od_search.tsv",
        g / f"{PREFIX}_FDR-0.05.betabinomial.{SNV_TAG}/b_and_sse.grad.txt",
        rtol=1e-12,
        columns=[0, 1],
    )
    run(
        "betabinomial_fdr.py",
        "--counts",
        testable,
        "--overdispersion",
        s / "od.tsv",
        "--all",
        s / "bb_all.tsv",
        "--output",
        s / "bb.tsv",
        "--report",
        s / "bb.txt",
    )
    check(
        "beta-binomial FDR report",
        s / "bb.txt",
        g / f"{PREFIX}_FDR-0.05.betabinom.{SNV_TAG}_cnt.txt",
        rtol=1e-12,
    )
    check(
        "beta-binomial p-values",
        s / "bb_all.tsv",
        g / f"{PREFIX}_counts.FDR-0.05.betabinom.{SNV_TAG}_cnt.tsv",
    )
    check(
        "beta-binomial significant SNVs",
        s / "bb.tsv",
        g / f"{PREFIX}_interestingHets.FDR-0.05.betabinom.{SNV_TAG}_cnt.tsv",
    )

    genes = g / "resources/gencode.v49_genes.bed"
    for kind, bam in (("uniq", "uniqreads"), ("mmap", "mmapreads")):
        regions = subprocess.Popen(
            ["cat", g / "hets_regions_hap1.bed", g / "hets_regions_hap2.bed"],
            stdout=subprocess.PIPE,
        )
        intersect = subprocess.Popen(
            [
                "bedtools",
                "intersect",
                "-a",
                "stdin",
                "-b",
                g / f"{PREFIX}_ASE-params_crdsorted_{bam}_over_hetSNVs.bam",
                "-split",
                "-wb",
                "-bed",
            ],
            stdin=regions.stdout,
            stdout=subprocess.PIPE,
        )
        run(
            "region_counts.py",
            "--kind",
            kind,
            "--hap1-bed",
            g / "hets_regions_hap1.bed",
            "--output",
            s / f"region_{kind}.tsv",
            stdin=intersect.stdout,
        )
        intersect.wait()
        check(
            f"region {kind} counts",
            s / f"region_{kind}.tsv",
            g / f"{POOL}_region_raw_counts_{kind}.tsv",
        )
    run(
        "filter_region_counts.py",
        "--counts",
        s / "region_uniq.tsv",
        "--mmap",
        s / "region_mmap.tsv",
        "--output",
        s / "region_filtered.tsv",
        "--log",
        s / "region_mmap.log",
    )
    check(
        "filtered region counts",
        s / "region_filtered.tsv",
        g / f"{POOL}_region_filtered_counts.tsv",
    )
    run(
        "select_testable.py",
        "--counts",
        s / "region_filtered.tsv",
        "--output",
        s / "region_testable.tsv",
        "--regions",
        genes,
    )
    check(
        "testable regions",
        s / "region_testable.tsv",
        g / f"{POOL}_region_filtered_counts.{SNV_TAG}_cnt.tsv",
    )
    run(
        "binomial_fdr.py",
        "--counts",
        s / "region_testable.tsv",
        "--report",
        s / "region_binom.txt",
        "--output",
        s / "region_binom.tsv",
    )
    check(
        "region binomial FDR report",
        s / "region_binom.txt",
        g / f"{POOL}_region_FDR-0.05.binom.{SNV_TAG}_cnt.txt",
        rtol=0,
        skip_prefix="#",
    )
    check(
        "region binomial significant",
        s / "region_binom.tsv",
        g / f"{POOL}_interesting_regions.FDR-0.05.binom.{SNV_TAG}_cnt.tsv",
    )
    run(
        "overdispersion.py",
        "--counts",
        s / "region_testable.tsv",
        "--summary",
        s / "region_od.tsv",
        "--search",
        s / "region_od_search.tsv",
        "--plot",
        s / "region_od.svg",
    )
    check(
        "region overdispersion search",
        s / "region_od_search.tsv",
        g / f"{POOL}_region_FDR-0.05.betabinomial.{SNV_TAG}/b_and_sse.grad.txt",
        rtol=1e-12,
    )
    run(
        "betabinomial_fdr.py",
        "--counts",
        s / "region_testable.tsv",
        "--overdispersion",
        s / "region_od.tsv",
        "--all",
        s / "region_bb_all.tsv",
        "--output",
        s / "region_bb.tsv",
        "--report",
        s / "region_bb.txt",
    )
    check(
        "region beta-binomial FDR report",
        s / "region_bb.txt",
        g / f"{POOL}_region_FDR-0.05.betabinom.{SNV_TAG}_cnt.txt",
        rtol=1e-12,
    )
    # Upstream wrote region p-values over the pooled SNV table (same file name).
    upstream = s / "region_bb_upstream.tsv"
    with (
        open(g / f"{POOL}_counts.FDR-0.05.betabinom.{SNV_TAG}_cnt.tsv") as src,
        open(upstream, "w") as dst,
    ):
        dst.writelines(
            "\t".join([r[0], r[-1]]) + "\n"
            for r in (line.rstrip("\n").split("\t") for line in src)
        )
    ours = s / "region_bb_pvalues.tsv"
    with open(s / "region_bb_all.tsv") as src, open(ours, "w") as dst:
        dst.writelines(
            "\t".join([r[0], r[-1]]) + "\n"
            for r in (line.rstrip("\n").split("\t") for line in src)
        )
    check("region beta-binomial p-values", ours, upstream, columns=[0, 1])

    print(f"{sum(results)}/{len(results)} comparisons agree. Outputs in {s}")
    sys.exit(0 if all(results) else 1)


if __name__ == "__main__":
    main(*sys.argv[1:])
