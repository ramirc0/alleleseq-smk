# Per-unit gene-level counts: reads over each gene's het SNVs, per haplotype.
# Summing hap1 against hap2 across a gene needs phased hets; see phasing.tsv.


rule region_hets:
    input:
        counts=f"{OUTDIR}/{{unit}}/snvs/{{unit}}.counts.filtered.tsv",
        regions=f"{REFDIR}/genome/genes.bed",
    output:
        hap1=f"{OUTDIR}/{{unit}}/genes/{{unit}}.hets.hap1.bed",
        hap2=f"{OUTDIR}/{{unit}}/genes/{{unit}}.hets.hap2.bed",
    log:
        f"{LOGDIR}/region_hets/{{unit}}.txt",
    benchmark:
        f"{BENCHDIR}/region_hets/{{unit}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/region_hets.py \
            --counts {input.counts:q} \
            --regions {input.regions:q} \
            --hap1 {output.hap1:q} \
            --hap2 {output.hap2:q}

        sync {output:q}
        """


rule region_counts:
    input:
        hap1=f"{OUTDIR}/{{unit}}/genes/{{unit}}.hets.hap1.bed",
        hap2=f"{OUTDIR}/{{unit}}/genes/{{unit}}.hets.hap2.bed",
        bams=lambda w: [f"{sample_prefix(s)}.hets.{w.kind}.bam" for s in samples_of(w.unit)],
    output:
        f"{OUTDIR}/{{unit}}/genes/{{unit}}.counts.raw_{{kind}}.tsv",
    resources:
        mem_mb=16000,
    log:
        f"{LOGDIR}/region_counts/{{unit}}.{{kind}}.txt",
    benchmark:
        f"{BENCHDIR}/region_counts/{{unit}}.{{kind}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        # -split: a read spliced across a het does not count for it.
        cat {input.hap1:q} {input.hap2:q} \
        | bedtools intersect -a stdin -b {input.bams:q} -split -wb -bed \
        | python {workflow.basedir:q}/scripts/region_counts.py \
            --kind {wildcards.kind} \
            --hap1-bed {input.hap1:q} \
            --output {output:q}

        sync {output:q}
        """


rule filter_region_counts:
    input:
        counts=f"{OUTDIR}/{{unit}}/genes/{{unit}}.counts.raw_uniq.tsv",
        mmap=f"{OUTDIR}/{{unit}}/genes/{{unit}}.counts.raw_mmap.tsv",
    output:
        counts=f"{OUTDIR}/{{unit}}/genes/{{unit}}.counts.filtered.tsv",
        log=f"{OUTDIR}/{{unit}}/genes/{{unit}}.mmap.log",
    log:
        f"{LOGDIR}/filter_region_counts/{{unit}}.txt",
    benchmark:
        f"{BENCHDIR}/filter_region_counts/{{unit}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/filter_region_counts.py \
            --counts {input.counts:q} \
            --mmap {input.mmap:q} \
            --output {output.counts:q} \
            --log {output.log:q}

        sync {output:q}
        """


rule name_genes:
    input:
        table=f"{OUTDIR}/{{unit}}/genes/{{unit}}.{{test}}.significant.tsv",
        names=f"{REFDIR}/genome/genes.tsv",
    output:
        f"{OUTDIR}/{{unit}}/genes/{{unit}}.{{test}}.significant.named.tsv",
    log:
        f"{LOGDIR}/name_genes/{{unit}}.{{test}}.txt",
    benchmark:
        f"{BENCHDIR}/name_genes/{{unit}}.{{test}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/name_genes.py \
            --table {input.table:q} \
            --names {input.names:q} \
            --output {output:q}

        sync {output:q}
        """
