# Config, sample sheets, paths and helpers shared by every rule file.

import re
from pathlib import Path

import polars as pl
from snakemake.exceptions import WorkflowError


configfile: "config/config.yaml"


RESULTS = config["outdir"]
RUN_ID = config.get("run_id", "default")
OUTDIR = f"{RESULTS}/{RUN_ID}"  # per-run
REFDIR = f"{RESULTS}/refs"  # run-independent
LOGDIR = f"logs/{RUN_ID}"
BENCHDIR = f"benchmarks/{RUN_ID}"
REFERENCE = config["reference"]["fasta"]
ANNOTATION = config["reference"]["annotation"]
VARIANT_TYPES = config["personal_genome"]["variant_types"]
VARIANT_SET = "-".join(VARIANT_TYPES)

# The workspace path is relative to the rule file. Every rule file sits in rules/.
SOFTWARE_ENV = pixi(workspace="../envs", env=config["pixi_env"], locked=True)


# The pixi plugin installs the env inside each job. Parallel SLURM jobs then race
# on a fresh or relocked env, so install it once before any job starts.
onstart:
    shell(
        "pixi install --manifest-path {workflow.basedir:q}/envs/pixi.toml --all --locked"
    )


def _read_sheet(path, key, required):
    """Read a TSV sheet with every cell as a string, keyed by its ID column.

    Parameters
    ----------
    path : str
        Sheet path.
    key : str
        ID column. Its values MUST be unique.
    required : list of str
        Columns the sheet must have besides `key`.

    Returns
    -------
    dict of str to dict
        One dict per row, keyed by its `key` value. Empty cells are None.
    """
    sheet = Path(path)
    if not sheet.exists():
        raise WorkflowError(f"Sheet not found: {sheet}. See config/*.example.tsv.")
    table = pl.read_csv(sheet, separator="\t", infer_schema_length=0, comment_prefix="#")
    missing = [c for c in [key, *required] if c not in table.columns]
    if missing:
        raise WorkflowError(f"{sheet} lacks columns: {', '.join(missing)}")
    if table.height == 0:
        raise WorkflowError(f"{sheet} has no rows.")
    repeated = table[key].filter(table[key].is_duplicated()).unique(maintain_order=True)
    if repeated.len():
        raise WorkflowError(f"{sheet} has duplicate {key}: {', '.join(repeated)}")
    return {
        row[key]: {k: (v or None) for k, v in row.items()}
        for row in table.iter_rows(named=True)
    }


INDIVIDUALS = _read_sheet(
    config["individuals"], "individual_id", ["vcf", "vcf_sample", "wgs_alignment"]
)
SAMPLES = _read_sheet(config["samples"], "sample_id", ["individual_id", "r1"])
GROUPS = {}
for _sample, _row in SAMPLES.items():
    if _row["individual_id"] not in INDIVIDUALS:
        raise WorkflowError(f"Sample {_sample}: unknown individual {_row['individual_id']}.")
    if _row.get("group"):
        GROUPS.setdefault(_row["group"], []).append(_sample)
for _group, _members in GROUPS.items():
    if _group in SAMPLES:
        raise WorkflowError(f"Group {_group} has the same ID as a sample.")
    if len({SAMPLES[s]["individual_id"] for s in _members}) > 1:
        raise WorkflowError(f"Group {_group} mixes individuals.")
UNITS = [*SAMPLES, *GROUPS]


def _alternatives(names):
    """Return a regex that matches exactly one of the names.

    Parameters
    ----------
    names : iterable of str
        Allowed values.

    Returns
    -------
    str
        Alternation of the escaped names.
    """
    return "|".join(re.escape(n) for n in sorted(names, key=len, reverse=True))


wildcard_constraints:
    individual=_alternatives(INDIVIDUALS),
    sample=_alternatives(SAMPLES),
    unit=_alternatives(UNITS),
    hap="1|2",
    kind="uniq|mmap",
    level="snvs|genes",
    stage="raw|filtered|testable",
    test="binomial|betabinomial",


def samples_of(unit):
    """Return the samples pooled in a unit.

    Parameters
    ----------
    unit : str
        Sample or group ID.

    Returns
    -------
    list of str
        The sample itself, or the group's members.
    """
    return [unit] if unit in SAMPLES else GROUPS[unit]


def individual_of(unit):
    """Return the individual a sample or group comes from.

    Parameters
    ----------
    unit : str
        Sample or group ID.

    Returns
    -------
    str
        Individual ID.
    """
    return SAMPLES[samples_of(unit)[0]]["individual_id"]


def pgenome(individual):
    """Return the personal genome directory of an individual.

    Parameters
    ----------
    individual : str
        Individual ID.

    Returns
    -------
    str
        ``results/refs/<individual>/<variant set>``.
    """
    return f"{REFDIR}/{individual}/{VARIANT_SET}"


def sample_prefix(sample):
    """Return the output prefix of a sample.

    Parameters
    ----------
    sample : str
        Sample ID.

    Returns
    -------
    str
        ``<outdir>/<run_id>/<sample>/<sample>``.
    """
    return f"{OUTDIR}/{sample}/{sample}"


def unit_prefix(unit, level):
    """Return the output prefix of a unit at one level.

    Parameters
    ----------
    unit : str
        Sample or group ID.
    level : str
        ``snvs`` or ``genes``.

    Returns
    -------
    str
        ``<outdir>/<run_id>/<unit>/<level>/<unit>``.
    """
    return f"{OUTDIR}/{unit}/{level}/{unit}"


def final_alignment(sample):
    """Return the alignment that later steps read.

    Parameters
    ----------
    sample : str
        Sample ID.

    Returns
    -------
    str
        The deduplicated BAM, or STAR's BAM when deduplication is off.
    """
    if config["remove_duplicates"]:
        return f"{sample_prefix(sample)}.dedup.bam"
    return f"{sample_prefix(sample)}.Aligned.sortedByCoord.out.bam"


def _contigs():
    """Return the reference's primary and extra contig names.

    Returns
    -------
    tuple of (list of str, list of str)
        Primary chromosomes (1-22, X, Y, M in either naming) and the rest.
    """
    fai = Path(f"{REFERENCE}.fai")
    if not fai.exists():
        raise WorkflowError(f"Index the reference first: samtools faidx {REFERENCE}")
    names = [line.split("\t")[0] for line in fai.read_text().splitlines()]
    canonical = {f"chr{c}" for c in [*range(1, 23), "X", "Y", "M"]}
    canonical |= {str(c) for c in [*range(1, 23), "X", "Y", "MT"]}
    primary = [n for n in names if n in canonical]
    return primary, [n for n in names if n not in canonical]


PRIMARY_CONTIGS, EXTRA_CONTIGS = _contigs()
MITO_CONTIGS = [c for c in PRIMARY_CONTIGS if c in ("chrM", "MT")]


def star_read_files(sample):
    """Return STAR's FASTQ arguments and decompression command.

    Parameters
    ----------
    sample : str
        Sample ID.

    Returns
    -------
    list of str
        ``--readFilesIn`` values followed by ``--readFilesCommand`` and its value.
    """
    row = SAMPLES[sample]
    reads = [row["r1"]] + ([row["r2"]] if row.get("r2") else [])
    command = "zcat" if row["r1"].endswith(".gz") else "cat"
    return [*reads, "--readFilesCommand", command]
