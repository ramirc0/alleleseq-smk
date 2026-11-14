# Per-unit het SNV counts and filters. A unit is a sample or a pooled group.


def unit_pileups(unit, kind):
    """Return the mpileups of every sample in a unit.

    Parameters
    ----------
    unit : str
        Sample or group ID.
    kind : str
        ``uniq`` or ``mmap``.

    Returns
    -------
    list of str
        Hap1 and hap2 pileups, sample by sample.
    """
    return [
        f"{sample_prefix(s)}.hap{hap}.{kind}.mpileup"
        for s in samples_of(unit)
        for hap in (1, 2)
    ]


rule count_alleles:
    input:
        pileups=lambda w: unit_pileups(w.unit, "uniq"),
        ref=lambda w: f"{pgenome(individual_of(w.unit))}/hets.ref.bed",
        hap1=lambda w: f"{pgenome(individual_of(w.unit))}/hets.hap1.bed",
        hap2=lambda w: f"{pgenome(individual_of(w.unit))}/hets.hap2.bed",
    output:
        counts=f"{OUTDIR}/{{unit}}/snvs/{{unit}}.counts.raw.tsv",
        discarded=f"{OUTDIR}/{{unit}}/snvs/{{unit}}.discarded.log",
    log:
        f"{LOGDIR}/count_alleles/{{unit}}.txt",
    benchmark:
        f"{BENCHDIR}/count_alleles/{{unit}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/count_alleles.py \
            --ref-bed {input.ref:q} \
            --hap1-bed {input.hap1:q} \
            --hap2-bed {input.hap2:q} \
            --discarded {output.discarded:q} \
            --output {output.counts:q} \
            --pileups {input.pileups:q}

        sync {output:q}
        """


rule filter_counts:
    input:
        counts=f"{OUTDIR}/{{unit}}/snvs/{{unit}}.counts.raw.tsv",
        depth=lambda w: f"{pgenome(individual_of(w.unit))}/wgs_read_depth.tsv",
        pileups=lambda w: unit_pileups(w.unit, "mmap"),
    output:
        counts=f"{OUTDIR}/{{unit}}/snvs/{{unit}}.counts.filtered.tsv",
        cnv=f"{OUTDIR}/{{unit}}/snvs/{{unit}}.discarded.cnv.log",
        warnings=f"{OUTDIR}/{{unit}}/snvs/{{unit}}.discarded.warnings.log",
        mmap=f"{OUTDIR}/{{unit}}/snvs/{{unit}}.mmap.log",
    log:
        f"{LOGDIR}/filter_counts/{{unit}}.txt",
    benchmark:
        f"{BENCHDIR}/filter_counts/{{unit}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/filter_counts.py \
            --counts {input.counts:q} \
            --read-depth {input.depth:q} \
            --mmap-pileups {input.pileups:q} \
            --output {output.counts:q} \
            --cnv-log {output.cnv:q} \
            --warning-log {output.warnings:q} \
            --mmap-log {output.mmap:q}

        sync {output:q}
        """
