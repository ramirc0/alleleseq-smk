# Run-independent references: the primary reference, gene regions and each
# individual's diploid personal genome, STAR index and WGS read depth.


rule primary_reference:
    input:
        fasta=REFERENCE,
    output:
        primary=f"{REFDIR}/genome/primary.fa",
        fai=f"{REFDIR}/genome/primary.fa.fai",
        extra=f"{REFDIR}/genome/extra_contigs.fa",
    params:
        primary=PRIMARY_CONTIGS,
        extra=EXTRA_CONTIGS,
    log:
        f"{LOGDIR}/primary_reference.txt",
    benchmark:
        f"{BENCHDIR}/primary_reference.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        samtools faidx {input.fasta:q} {params.primary:q} > {output.primary:q}
        samtools faidx {output.primary:q}
        : > {output.extra:q}
        if [ -n "{params.extra}" ]; then
            samtools faidx {input.fasta:q} -r <(printf '%s\n' {params.extra:q}) > {output.extra:q}
        fi

        sync {output:q}
        """


rule genes_bed:
    input:
        gtf=ANNOTATION,
    output:
        bed=f"{REFDIR}/genome/genes.bed",
        names=f"{REFDIR}/genome/genes.tsv",
    log:
        f"{LOGDIR}/genes_bed.txt",
    benchmark:
        f"{BENCHDIR}/genes_bed.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/genes_bed.py \
            --gtf {input.gtf:q} \
            --bed {output.bed:q} \
            --names {output.names:q}

        sync {output:q}
        """


rule subset_vcf:
    input:
        vcf=lambda w: INDIVIDUALS[w.individual]["vcf"],
    output:
        vcf=f"{REFDIR}/{{individual}}/variants.vcf.gz",
        index=f"{REFDIR}/{{individual}}/variants.vcf.gz.tbi",
    params:
        sample=lambda w: INDIVIDUALS[w.individual]["vcf_sample"],
        types=",".join(VARIANT_TYPES),
    log:
        f"{LOGDIR}/subset_vcf/{{individual}}.txt",
    benchmark:
        f"{BENCHDIR}/subset_vcf/{{individual}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bcftools view \
            --no-version \
            --samples {params.sample:q} \
            --apply-filters PASS \
            --min-alleles 2 \
            --max-alleles 2 \
            --types {params.types:q} \
            {input.vcf:q} \
        | bcftools view \
            --no-version \
            --exclude 'GT="./." || GT="mis" || GT="ref"' \
            --output-type z \
            --output {output.vcf:q}
        bcftools index --tbi {output.vcf:q}

        sync {output:q}
        """


rule phased_hets:
    input:
        vcf=lambda w: INDIVIDUALS[w.individual]["phased_vcf"],
    output:
        vcf=f"{REFDIR}/{{individual}}/phased_hets.vcf.gz",
    params:
        sample=lambda w: INDIVIDUALS[w.individual]["vcf_sample"],
    threads: 4
    log:
        f"{LOGDIR}/phased_hets/{{individual}}.txt",
    benchmark:
        f"{BENCHDIR}/phased_hets/{{individual}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bcftools view \
            --no-version \
            --threads {threads} \
            --samples {params.sample:q} \
            --output-type u \
            {input.vcf:q} \
        | bcftools view \
            --no-version \
            --include 'GT="het"' \
            --output-type z \
            --output {output.vcf:q}
        echo "Phased hets: $(bcftools view --no-header {output.vcf:q} | wc -l)"

        sync {output:q}
        """


rule orient_hets:
    input:
        vcf=f"{REFDIR}/{{individual}}/variants.vcf.gz",
        phased=lambda w: (
            [f"{REFDIR}/{w.individual}/phased_hets.vcf.gz"]
            if INDIVIDUALS[w.individual].get("phased_vcf")
            else []
        ),
    output:
        vcf=f"{REFDIR}/{{individual}}/oriented.vcf.gz",
        index=f"{REFDIR}/{{individual}}/oriented.vcf.gz.csi",
        orientations=f"{REFDIR}/{{individual}}/orientations.tsv.gz",
        summary=f"{REFDIR}/{{individual}}/phasing.tsv",
    params:
        phased=lambda w, input: ["--phased", input.phased[0]] if input.phased else [],
        min_confidence=config["personal_genome"]["min_phase_confidence"],
        seed=config["personal_genome"]["orientation_seed"],
        plain=f"{REFDIR}/{{individual}}/oriented.vcf",
    log:
        f"{LOGDIR}/orient_hets/{{individual}}.txt",
    benchmark:
        f"{BENCHDIR}/orient_hets/{{individual}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/orient_hets.py \
            --vcf {input.vcf:q} \
            {params.phased:q} \
            --min-phase-confidence {params.min_confidence:q} \
            --seed {params.seed:q} \
            --output {params.plain:q} \
            --orientations {output.orientations:q} \
            --summary {output.summary:q}
        bcftools view --no-version --output-type z --output {output.vcf:q} {params.plain:q}
        rm {params.plain:q}
        bcftools index --csi {output.vcf:q}
        cat {output.summary:q}

        sync {output:q}
        """


rule haplotype:
    input:
        fasta=f"{REFDIR}/genome/primary.fa",
        vcf=f"{REFDIR}/{{individual}}/oriented.vcf.gz",
        index=f"{REFDIR}/{{individual}}/oriented.vcf.gz.csi",
    output:
        fasta=f"{REFDIR}/{{individual}}/hap{{hap}}.fa.gz",
        fai=f"{REFDIR}/{{individual}}/hap{{hap}}.fa.gz.fai",
        gzi=f"{REFDIR}/{{individual}}/hap{{hap}}.fa.gz.gzi",
        chain=f"{REFDIR}/{{individual}}/hap{{hap}}.chain",
    params:
        sample=lambda w: INDIVIDUALS[w.individual]["vcf_sample"],
        # Mitochondria are inherited from the mother only, so hap1 (the
        # paternal slot in vcf2diploid) leaves chrM out.
        drop=lambda w: MITO_CONTIGS if w.hap == "1" else [],
        raw_chain=f"{REFDIR}/{{individual}}/hap{{hap}}.raw.chain",
    log:
        f"{LOGDIR}/haplotype/{{individual}}.hap{{hap}}.txt",
    benchmark:
        f"{BENCHDIR}/haplotype/{{individual}}.hap{{hap}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        bcftools consensus \
            --fasta-ref {input.fasta:q} \
            --samples {params.sample:q} \
            --haplotype {wildcards.hap} \
            --chain {params.raw_chain:q} \
            {input.vcf:q} \
        | python {workflow.basedir:q}/scripts/rename_contigs.py \
            --format fasta \
            --suffix _hap{wildcards.hap} \
            --drop {params.drop:q} \
        | bgzip > {output.fasta:q}
        python {workflow.basedir:q}/scripts/rename_contigs.py \
            --format chain \
            --suffix _hap{wildcards.hap} \
            < {params.raw_chain:q} \
            > {output.chain:q}
        rm {params.raw_chain:q}
        samtools faidx {output.fasta:q}

        sync {output:q}
        """


rule het_snvs:
    input:
        vcf=f"{REFDIR}/{{individual}}/oriented.vcf.gz",
        fai=f"{REFDIR}/genome/primary.fa.fai",
        hap1_chain=f"{REFDIR}/{{individual}}/hap1.chain",
        hap2_chain=f"{REFDIR}/{{individual}}/hap2.chain",
    output:
        ref=f"{REFDIR}/{{individual}}/hets.ref.bed",
        hap1=f"{REFDIR}/{{individual}}/hets.hap1.bed",
        hap2=f"{REFDIR}/{{individual}}/hets.hap2.bed",
        hap1_unlifted=f"{REFDIR}/{{individual}}/hets.hap1.not_lifted.bed",
        hap2_unlifted=f"{REFDIR}/{{individual}}/hets.hap2.not_lifted.bed",
    params:
        sample=lambda w: INDIVIDUALS[w.individual]["vcf_sample"],
    log:
        f"{LOGDIR}/het_snvs/{{individual}}.txt",
    benchmark:
        f"{BENCHDIR}/het_snvs/{{individual}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/hetsnv_bed.py \
            --vcf {input.vcf:q} \
            --sample {params.sample:q} \
            --fai {input.fai:q} \
            --output {output.ref:q}
        liftOver {output.ref:q} {input.hap1_chain:q} {output.hap1:q} {output.hap1_unlifted:q}
        liftOver {output.ref:q} {input.hap2_chain:q} {output.hap2:q} {output.hap2_unlifted:q}

        sync {output:q}
        """


rule diploid_annotation:
    input:
        gtf=ANNOTATION,
        hap1_chain=f"{REFDIR}/{{individual}}/hap1.chain",
        hap2_chain=f"{REFDIR}/{{individual}}/hap2.chain",
    output:
        gtf=f"{REFDIR}/{{individual}}/diploid.gtf",
        hap1_unlifted=f"{REFDIR}/{{individual}}/annotation.hap1.not_lifted",
        hap2_unlifted=f"{REFDIR}/{{individual}}/annotation.hap2.not_lifted",
        hap1=temp(f"{REFDIR}/{{individual}}/annotation.hap1.gtf"),
        hap2=temp(f"{REFDIR}/{{individual}}/annotation.hap2.gtf"),
    log:
        f"{LOGDIR}/diploid_annotation/{{individual}}.txt",
    benchmark:
        f"{BENCHDIR}/diploid_annotation/{{individual}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        liftOver -gff {input.gtf:q} {input.hap1_chain:q} {output.hap1:q} {output.hap1_unlifted:q}
        liftOver -gff {input.gtf:q} {input.hap2_chain:q} {output.hap2:q} {output.hap2_unlifted:q}
        # Upstream concatenates the maternal (hap2) annotation first.
        cat {output.hap2:q} {output.hap1:q} > {output.gtf:q}

        sync {output:q}
        """


rule star_index:
    input:
        hap1=f"{REFDIR}/{{individual}}/hap1.fa.gz",
        hap2=f"{REFDIR}/{{individual}}/hap2.fa.gz",
        extra=f"{REFDIR}/genome/extra_contigs.fa",
    output:
        # 54 GB per individual. Snakemake deletes it once the individual's samples are aligned.
        temp(directory(f"{REFDIR}/{{individual}}/star_index")),
    params:
        ram=config["star"]["genome_generate_ram"],
        sa_index_nbases=config["star"]["sa_index_nbases"],
    threads: 32
    log:
        f"{LOGDIR}/star_index/{{individual}}.txt",
    benchmark:
        f"{BENCHDIR}/star_index/{{individual}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        # STAR reads only plain FASTA.
        mkdir -p {output:q}/fasta
        bgzip --decompress --stdout {input.hap1:q} > {output:q}/fasta/hap1.fa
        bgzip --decompress --stdout {input.hap2:q} > {output:q}/fasta/hap2.fa
        fastas=({output:q}/fasta/hap1.fa {output:q}/fasta/hap2.fa)
        # Unplaced and decoy contigs soak up reads that would otherwise mismap.
        if [ -s {input.extra:q} ]; then fastas+=({input.extra:q}); fi
        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --limitGenomeGenerateRAM {params.ram} \
            --genomeSAindexNbases {params.sa_index_nbases} \
            --genomeDir {output:q} \
            --genomeFastaFiles "${{fastas[@]}}" \
            --outFileNamePrefix {output:q}/
        rm -r {output:q}/fasta

        sync {output:q}
        """


rule wgs_read_depth:
    input:
        hets=f"{REFDIR}/{{individual}}/hets.ref.bed",
        alignment=lambda w: INDIVIDUALS[w.individual]["wgs_alignment"],
        reference=REFERENCE,
    output:
        depth=f"{REFDIR}/{{individual}}/wgs_read_depth.tsv",
        windows=temp(f"{REFDIR}/{{individual}}/wgs_windows.bed"),
        genome=temp(f"{REFDIR}/{{individual}}/wgs_genome.txt"),
        header=temp(f"{REFDIR}/{{individual}}/wgs_header.sam"),
    threads: 4
    log:
        f"{LOGDIR}/wgs_read_depth/{{individual}}.txt",
    benchmark:
        f"{BENCHDIR}/wgs_read_depth/{{individual}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        samtools view --header-only {input.alignment:q} > {output.header:q}
        python {workflow.basedir:q}/scripts/wgs_windows.py \
            --hets {input.hets:q} \
            --header {output.header:q} \
            --windows {output.windows:q} \
            --genome {output.genome:q}
        bedtools intersect \
            -a <(samtools view --threads {threads} -u --reference {input.reference:q} {input.alignment:q}) \
            -b {output.windows:q} \
            -g {output.genome:q} \
            -sorted \
            -wb \
            -bed \
        | python {workflow.basedir:q}/scripts/read_depth.py \
            --windows {output.windows:q} \
            --output {output.depth:q}

        sync {output:q}
        """
