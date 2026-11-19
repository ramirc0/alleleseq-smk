# alleleseq-smk

AlleleSeq2 allele-specific expression (ASE) analysis as a Snakemake workflow, rewritten in Python 3.

## Description

The workflow builds a diploid personal genome for each individual, aligns RNA-seq reads to both haplotypes with STAR, and counts reads per allele at heterozygous SNVs.
It then calls allele-specific SNVs and genes with the binomial and beta-binomial tests of AlleleSeq2 (Rozowsky et al. 2011, Chen et al. 2016).
Every step of the upstream Makefiles, Python 2 scripts and AlleleDB R scripts is ported to Python.
Only R's `dbinom`, `pbinom` and VGAM's `dbetabinom` are called, through rpy2.
On a full-size upstream AlleleSeq2 run, the ported steps reproduce every count table, FDR cutoff and significant set (`pixi run golden`).

### Layout

```
.tests/
  unit/                  pytest tests of the scripts
  integration/           synthetic data with planted imbalance, made by .tests/make_toy_data.py
  golden/                replays an upstream AlleleSeq2 run
config/                  config.yaml.template, individuals and samples examples
workflow/rules/          common, genome, align, hetsnvs, regions, stats
workflow/scripts/        one script per step, plus shared modules
profiles/                local and SLURM execution
docs/                    workflow diagrams: workflow.gv (every rule), workflow-condensed.gv
```

### Steps

1. Personal genome, once per individual and run (`results/<run_id>/refs/<individual>/`):
   - Subset the VCF to one sample's PASS biallelic SNVs (and indels).
   - Take each het's phase from `phased_vcf` where it calls the same het with `PP` ≥ `min_phase_confidence` or no `PP`. Genotypes stay those of the VCF.
   - Orient the remaining hets with a seeded random draw, one per `PS` phase set and one per unphased het.
   - Build both haplotypes and their chains with `bcftools consensus`. This replaces vcf2diploid.
   - Lift het SNVs and the GENCODE annotation onto each haplotype.
   - Build a STAR index of both haplotypes plus the unplaced contigs.
   - Compute each het's WGS read depth ratio in a 2 kb window.
2. Per sample: align with STAR in two-pass mode, remove duplicates with Picard, split the reads over hets into unique and multi-mapping, and pile them up on each haplotype.
3. Per unit (a sample, or a `group` of samples pooled):
   - Count alleles, then drop SNVs with a WGS depth ratio outside [0.5, 1.5] or with phase warnings.
   - Add weaker-allele multi-mapping reads to the weaker allele's count.
   - Run the simulated binomial FDR and the AlleleDB beta-binomial FDR.
4. Per unit, at the gene level: sum distinct reads per haplotype over each gene's hets, then run the same tests.

### Outputs

Per unit under `results/<run_id>/<unit>/`, `snvs/` and `genes/` each hold:

| File | Content |
|---|---|
| `<unit>.counts.raw.tsv` | allele counts before filtering (`raw_uniq` for genes) |
| `<unit>.counts.filtered.tsv` | after the depth, phase-warning and multi-mapping filters |
| `<unit>.counts.testable.tsv` | autosomes (plus `keep_chromosomes`) with enough reads |
| `<unit>.binomial.{fdr.txt,significant.tsv}` | simulated binomial FDR report and calls |
| `<unit>.betabinomial.{tsv,fdr.txt,significant.tsv}` | beta-binomial p-values, FDR report and calls |
| `<unit>.overdispersion.{tsv,search.tsv,svg,png}` | fitted overdispersion and its search |
| `<unit>.allelic_ratio.{raw,filtered,testable}.{svg,png}` | allele ratio histograms |

Gene tables use `gene_id`. The `*.significant.named.tsv` copies add gene names and types.

### Differences from upstream AlleleSeq2

- A het SNV inside overlapping genes counts for each of them. Upstream kept only the last gene.
- Pooled pileups sum `cN` across files. A het without a WGS window is dropped and logged instead of crashing.
- The region tables keep their own names and correct column labels. Upstream overwrote the pooled SNV beta-binomial table.
- AlleleDB's minimum depth follows `min_total`. The default of 6 leaves results unchanged.

## Getting Started

### Dependencies

- Linux x86-64 and [pixi](https://pixi.sh).
- Network access to GitHub. Snakemake, its plugins and poikilos install from git.
- Per individual:
  - a VCF with the sample's genotypes (a joint VCF works)
  - a WGS alignment (CRAM or BAM)
- Per sample: RNA-seq FASTQs.
- An indexed reference FASTA (`samtools faidx`) and its GENCODE GTF.
- About 200 GB of RAM for each STAR index, and about 100 GB for each alignment.

### Installing

```bash
git clone <repo> alleleseq-smk && cd alleleseq-smk
pixi install
cp config/config.yaml.template config/config.yaml
cp config/individuals.example.tsv config/individuals.tsv
cp config/samples.example.tsv config/samples.tsv
```

Edit the three files:

- `individuals.tsv` maps each individual to a VCF, the VCF sample column and a WGS alignment. An optional `phased_vcf` column names a statistically phased VCF, such as SHAPEIT5's, with the same sample name.
- `samples.tsv` maps each RNA sample to an individual and its FASTQs. Leave `r2` empty for single-end data. Samples that share a `group` are also pooled.

### Executing program

```bash
pixi run snakemake -n -p --profile profiles/local   # dry run
pixi run slurm <run_id>                             # SLURM, sets run_id
pixi run unit                                       # unit tests
pixi run test                                       # integration test on toy data
pixi run golden --run <dir> --rna <PREFIX> --pool <PREFIX> --vcf-sample <id> --pgenome <dir>  # replay an upstream run
pixi run lint
```

## Help

- Gene-level calls need phased hets. Each individual's `phasing.tsv` reports the phased fraction. Unphased hets are oriented at random, so hap1 and hap2 sums across a gene mix both haplotypes. Phase the cohort first, for example with SHAPEIT5, and set `phased_vcf`.
- `variant_types: [snps]` writes only SNVs into the haplotypes. `[snps, indels]` also writes indels, which reduces mapping bias near them.
- Personal genomes belong to their run. Use a new `run_id` for a new VCF or `variant_types`, so earlier runs keep their STAR indexes.
- The profiles ignore software-environment changes when deciding what to rerun, so an env edit does not rebuild every STAR index. After a tool upgrade, use `--forcerun <rule>`.
- STAR refuses to start when its temp directory exists. The `align` rule removes it first.

## Authors

Christian Ramirez ([@ramirc0](https://github.com/ramirc0))

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE).

## Acknowledgments

- AlleleSeq2 (gersteinlab), AlleleSeq (Rozowsky et al. 2011, doi:10.1038/msb.2011.54) and AlleleDB (Chen et al. 2016, doi:10.1038/ncomms11101).
- [snakemake-template](https://github.com/ramirc0/snakemake-template).
