# Statistics shared by het SNVs and genes: testable rows, the binomial FDR,
# the AlleleDB beta-binomial FDR and allele ratio plots.


def raw_counts(unit, level):
    """Return the raw counts table of a unit at one level.

    Parameters
    ----------
    unit : str
        Sample or group ID.
    level : str
        ``snvs`` or ``genes``.

    Returns
    -------
    str
        Path of the raw counts TSV.
    """
    suffix = "raw" if level == "snvs" else "raw_uniq"
    return f"{unit_prefix(unit, level)}.counts.{suffix}.tsv"


rule select_testable:
    input:
        counts=f"{OUTDIR}/{{unit}}/{{level}}/{{unit}}.counts.filtered.tsv",
        regions=lambda w: [f"{REFDIR}/genome/genes.bed"] if w.level == "genes" else [],
    output:
        f"{OUTDIR}/{{unit}}/{{level}}/{{unit}}.counts.testable.tsv",
    params:
        min_total=config["counts"]["min_total"],
        min_per_allele=config["counts"]["min_per_allele"],
        keep=config["counts"]["keep_chromosomes"],
        regions=lambda w, input: ["--regions", *input.regions] if input.regions else [],
    log:
        f"{LOGDIR}/select_testable/{{unit}}.{{level}}.txt",
    benchmark:
        f"{BENCHDIR}/select_testable/{{unit}}.{{level}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/select_testable.py \
            --counts {input.counts:q} \
            --output {output:q} \
            --min-total {params.min_total} \
            --min-per-allele {params.min_per_allele} \
            --keep-chromosomes {params.keep:q} \
            {params.regions:q}

        sync {output:q}
        """


rule binomial_fdr:
    input:
        f"{OUTDIR}/{{unit}}/{{level}}/{{unit}}.counts.testable.tsv",
    output:
        report=f"{OUTDIR}/{{unit}}/{{level}}/{{unit}}.binomial.fdr.txt",
        significant=f"{OUTDIR}/{{unit}}/{{level}}/{{unit}}.binomial.significant.tsv",
    params:
        fdr=config["stats"]["fdr"],
        simulations=config["stats"]["fdr_simulations"],
    log:
        f"{LOGDIR}/binomial_fdr/{{unit}}.{{level}}.txt",
    benchmark:
        f"{BENCHDIR}/binomial_fdr/{{unit}}.{{level}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/binomial_fdr.py \
            --counts {input:q} \
            --report {output.report:q} \
            --output {output.significant:q} \
            --fdr {params.fdr} \
            --simulations {params.simulations}

        sync {output:q}
        """


rule overdispersion:
    input:
        f"{OUTDIR}/{{unit}}/{{level}}/{{unit}}.counts.testable.tsv",
    output:
        summary=f"{OUTDIR}/{{unit}}/{{level}}/{{unit}}.overdispersion.tsv",
        search=f"{OUTDIR}/{{unit}}/{{level}}/{{unit}}.overdispersion.search.tsv",
        plot=f"{OUTDIR}/{{unit}}/{{level}}/{{unit}}.overdispersion.svg",
        plot_png=f"{OUTDIR}/{{unit}}/{{level}}/{{unit}}.overdispersion.png",
    params:
        min_total=config["counts"]["min_total"],
    log:
        f"{LOGDIR}/overdispersion/{{unit}}.{{level}}.txt",
    benchmark:
        f"{BENCHDIR}/overdispersion/{{unit}}.{{level}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/overdispersion.py \
            --counts {input:q} \
            --min-total {params.min_total} \
            --summary {output.summary:q} \
            --search {output.search:q} \
            --plot {output.plot:q}

        sync {output:q}
        """


rule betabinomial_fdr:
    input:
        counts=f"{OUTDIR}/{{unit}}/{{level}}/{{unit}}.counts.testable.tsv",
        overdispersion=f"{OUTDIR}/{{unit}}/{{level}}/{{unit}}.overdispersion.tsv",
    output:
        tested=f"{OUTDIR}/{{unit}}/{{level}}/{{unit}}.betabinomial.tsv",
        significant=f"{OUTDIR}/{{unit}}/{{level}}/{{unit}}.betabinomial.significant.tsv",
        report=f"{OUTDIR}/{{unit}}/{{level}}/{{unit}}.betabinomial.fdr.txt",
    params:
        fdr=config["stats"]["fdr"],
        min_total=config["counts"]["min_total"],
    log:
        f"{LOGDIR}/betabinomial_fdr/{{unit}}.{{level}}.txt",
    benchmark:
        f"{BENCHDIR}/betabinomial_fdr/{{unit}}.{{level}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/betabinomial_fdr.py \
            --counts {input.counts:q} \
            --overdispersion {input.overdispersion:q} \
            --fdr {params.fdr} \
            --min-total {params.min_total} \
            --all {output.tested:q} \
            --output {output.significant:q} \
            --report {output.report:q}

        sync {output:q}
        """


rule plot_allelic_ratio:
    input:
        lambda w: raw_counts(w.unit, w.level)
        if w.stage == "raw"
        else f"{unit_prefix(w.unit, w.level)}.counts.{w.stage}.tsv",
    output:
        plot=f"{OUTDIR}/{{unit}}/{{level}}/{{unit}}.allelic_ratio.{{stage}}.svg",
        plot_png=f"{OUTDIR}/{{unit}}/{{level}}/{{unit}}.allelic_ratio.{{stage}}.png",
    params:
        title=lambda w: f"{w.unit} {w.level}, {w.stage}",
    log:
        f"{LOGDIR}/plot_allelic_ratio/{{unit}}.{{level}}.{{stage}}.txt",
    benchmark:
        f"{BENCHDIR}/plot_allelic_ratio/{{unit}}.{{level}}.{{stage}}.tsv"
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/plot_allelic_ratio.py \
            --counts {input:q} \
            --title {params.title:q} \
            --plot {output.plot:q}

        sync {output:q}
        """
