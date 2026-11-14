# Per-sample alignment to the diploid personal genome and pileups at het SNVs.


rule fastqc:
    input:
        reads=lambda w: [SAMPLES[w.sample]["r1"]]
        + ([SAMPLES[w.sample]["r2"]] if SAMPLES[w.sample].get("r2") else []),
    output:
        directory(f"{OUTDIR}/{{sample}}/fastqc"),
    threads: 2
    log:
        f"{LOGDIR}/fastqc/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/fastqc/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        mkdir -p {output:q}
        fastqc --threads {threads} --outdir {output:q} {input.reads:q}

        sync {output:q}
        """


rule align:
    input:
        reads=lambda w: [SAMPLES[w.sample]["r1"]]
        + ([SAMPLES[w.sample]["r2"]] if SAMPLES[w.sample].get("r2") else []),
        index=lambda w: f"{pgenome(SAMPLES[w.sample]['individual_id'])}/star_index",
        gtf=lambda w: f"{pgenome(SAMPLES[w.sample]['individual_id'])}/diploid.gtf",
    output:
        bam=f"{OUTDIR}/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam"
        if not config["remove_duplicates"]
        else temp(f"{OUTDIR}/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam"),
        bai=f"{OUTDIR}/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam.bai",
        flagstat=f"{OUTDIR}/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.flagstat",
        final=f"{OUTDIR}/{{sample}}/{{sample}}.Log.final.out",
        junctions=f"{OUTDIR}/{{sample}}/{{sample}}.SJ.out.tab",
    params:
        prefix=lambda w: f"{sample_prefix(w.sample)}.",
        reads=lambda w: star_read_files(w.sample),
        read_group=lambda w: [
            f"ID:{w.sample}",
            f"SM:{INDIVIDUALS[SAMPLES[w.sample]['individual_id']]['vcf_sample']}",
            f"LB:{w.sample}",
            "PL:ILLUMINA",
        ],
        star=config["star"],
    threads: 16
    log:
        f"{LOGDIR}/align/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/align/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        # A leftover STAR temp directory makes STAR refuse to start.
        rm -rf {params.prefix:q}_STARtmp
        STAR \
            --runThreadN {threads} \
            --genomeDir {input.index:q} \
            --twopassMode Basic \
            --readFilesIn {params.reads:q} \
            --outFileNamePrefix {params.prefix:q} \
            --outSAMattrRGline {params.read_group:q} \
            --outSAMattributes All \
            --outFilterMismatchNoverReadLmax {params.star[mismatch_n_over_read_lmax]} \
            --outFilterMatchNminOverLread {params.star[match_n_min_over_lread]} \
            --outFilterMultimapNmax 999999 \
            --scoreGenomicLengthLog2scale 0.0 \
            --sjdbOverhang {params.star[sjdb_overhang]} \
            --sjdbGTFfile {input.gtf:q} \
            --limitSjdbInsertNsj {params.star[limit_sjdb_insert_nsj]} \
            --limitBAMsortRAM {params.star[limit_bam_sort_ram]} \
            --outSAMtype BAM SortedByCoordinate
        rm -rf {params.prefix:q}_STARgenome {params.prefix:q}_STARpass1
        samtools index {output.bam:q}
        samtools flagstat {output.bam:q} > {output.flagstat:q}

        sync {output:q}
        """


rule remove_duplicates:
    input:
        bam=f"{OUTDIR}/{{sample}}/{{sample}}.Aligned.sortedByCoord.out.bam",
    output:
        bam=f"{OUTDIR}/{{sample}}/{{sample}}.dedup.bam",
        bai=f"{OUTDIR}/{{sample}}/{{sample}}.dedup.bam.bai",
        metrics=f"{OUTDIR}/{{sample}}/{{sample}}.dedup.metrics",
        flagstat=f"{OUTDIR}/{{sample}}/{{sample}}.dedup.flagstat",
    params:
        heap=lambda w, resources: f"-Xmx{int(resources.mem_mb * 0.85)}m",
    resources:
        mem_mb=16000,
    log:
        f"{LOGDIR}/remove_duplicates/{{sample}}.txt",
    benchmark:
        f"{BENCHDIR}/remove_duplicates/{{sample}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        picard {params.heap} MarkDuplicates \
            --INPUT {input.bam:q} \
            --OUTPUT {output.bam:q} \
            --METRICS_FILE {output.metrics:q} \
            --REMOVE_DUPLICATES true \
            --DUPLICATE_SCORING_STRATEGY SUM_OF_BASE_QUALITIES \
            --TMP_DIR {resources.tmpdir:q}
        samtools index {output.bam:q}
        samtools flagstat {output.bam:q} > {output.flagstat:q}

        sync {output:q}
        """


rule het_alignments:
    input:
        bam=lambda w: final_alignment(w.sample),
        hap1=lambda w: f"{pgenome(SAMPLES[w.sample]['individual_id'])}/hets.hap1.bed",
        hap2=lambda w: f"{pgenome(SAMPLES[w.sample]['individual_id'])}/hets.hap2.bed",
    output:
        bam=f"{OUTDIR}/{{sample}}/{{sample}}.hets.{{kind}}.bam",
        bai=f"{OUTDIR}/{{sample}}/{{sample}}.hets.{{kind}}.bam.bai",
        flagstat=f"{OUTDIR}/{{sample}}/{{sample}}.hets.{{kind}}.flagstat",
    params:
        # STAR gives MAPQ 255 to uniquely mapped reads only.
        mapq=lambda w: ["--min-MQ", "255"] if w.kind == "uniq" else ["--expr", "mapq!=255"],
    log:
        f"{LOGDIR}/het_alignments/{{sample}}.{{kind}}.txt",
    benchmark:
        f"{BENCHDIR}/het_alignments/{{sample}}.{{kind}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        cat {input.hap1:q} {input.hap2:q} \
        | samtools view \
            --bam \
            --with-header \
            --target-file - \
            {params.mapq:q} \
            --output {output.bam:q} \
            {input.bam:q}
        samtools index {output.bam:q}
        samtools flagstat {output.bam:q} > {output.flagstat:q}

        sync {output:q}
        """


rule pileup:
    input:
        bam=f"{OUTDIR}/{{sample}}/{{sample}}.hets.{{kind}}.bam",
        bai=f"{OUTDIR}/{{sample}}/{{sample}}.hets.{{kind}}.bam.bai",
        fasta=lambda w: f"{pgenome(SAMPLES[w.sample]['individual_id'])}/hap{w.hap}.fa",
        fai=lambda w: f"{pgenome(SAMPLES[w.sample]['individual_id'])}/hap{w.hap}.fa.fai",
        positions=lambda w: f"{pgenome(SAMPLES[w.sample]['individual_id'])}/hets.hap{w.hap}.bed",
    output:
        f"{OUTDIR}/{{sample}}/{{sample}}.hap{{hap}}.{{kind}}.mpileup",
    log:
        f"{LOGDIR}/pileup/{{sample}}.hap{{hap}}.{{kind}}.txt",
    benchmark:
        f"{BENCHDIR}/pileup/{{sample}}.hap{{hap}}.{{kind}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        # --ff UNMAP keeps secondary, QC-failed and duplicate reads, as upstream.
        samtools mpileup \
            --no-BAQ \
            --min-BQ 0 \
            --max-depth 999999 \
            --ff UNMAP \
            --fasta-ref {input.fasta:q} \
            --positions {input.positions:q} \
            --output {output:q} \
            {input.bam:q}

        sync {output:q}
        """
