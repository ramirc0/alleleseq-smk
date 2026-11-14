"""Primary human chromosome names in the two common conventions."""

UCSC = [f"chr{c}" for c in [*range(1, 23), "X", "Y", "M"]]
ENSEMBL = [str(c) for c in [*range(1, 23), "X", "Y", "MT"]]
AUTOSOMES = [f"chr{c}" for c in range(1, 23)] + [str(c) for c in range(1, 23)]


def primary_names(names):
    """Return the primary chromosomes among contig names, in input order.

    Parameters
    ----------
    names : iterable of str
        Contig names, for example from a FASTA index or a BAM header.

    Returns
    -------
    list of str
        Names found in either convention.
    """
    return [name for name in names if name in UCSC or name in ENSEMBL]


def harmonize(name, target):
    """Rename a chromosome to the convention used by `target`.

    Parameters
    ----------
    name : str
        Chromosome name to rename.
    target : list of str
        Primary chromosome names of the destination, as from `primary_names`.

    Returns
    -------
    str or None
        The matching name in `target`, or None when there is none.
    """
    if name in target:
        return name
    if not target:
        return None
    if name in ENSEMBL and target[0] in UCSC:
        return UCSC[ENSEMBL.index(name)]
    if name in UCSC and target[0] in ENSEMBL:
        return ENSEMBL[UCSC.index(name)]
    if "M" in name:
        mito = [t for t in target if "M" in t]
        return mito[-1] if mito else None
    return None


def is_testable(name, keep):
    """Tell whether a chromosome enters the allele-specific tests.

    Parameters
    ----------
    name : str
        Chromosome name.
    keep : list of str
        Extra chromosomes to keep, such as ``["X"]``. Both conventions match.

    Returns
    -------
    bool
        True for autosomes and for chromosomes listed in `keep`.
    """
    extra = [f"chr{c}" for c in keep] + list(keep)
    return name in AUTOSOMES or name in extra
