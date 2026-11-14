"""Base counts from ``samtools mpileup`` text output.

Each line holds a contig, a 1-based position, the reference base, the depth,
the read bases and their qualities. Matches show as ``.`` or ``,``,
mismatches as ``ACGTN`` in either case, deletions as ``*`` (or ``#``), and
reference skips as ``<`` or ``>``. ``^`` plus one MAPQ character marks a read
start and ``$`` a read end. ``+3ACG`` and ``-3ACG`` describe an indel after
the position.
"""

from dataclasses import dataclass

BASES = ("A", "C", "G", "T", "N")
DELETION_MARKS = "*#"
SKIP_MARKS = "<>"


@dataclass
class SiteCounts:
    """Base counts at one haplotype position.

    Attributes
    ----------
    counts : dict of str to int
        Reads per base, keyed by ``A``, ``C``, ``G``, ``T`` and ``N``.
    warning : str
        ``.`` when only the haplotype base was seen. Otherwise
        ``<n>_other_alleles`` or ``hap_allele_not_largest_cnt``, joined by
        commas across pooled pileup files.
    allele : str
        The haplotype's base at this position, upper case.
    """

    counts: dict
    warning: str
    allele: str


def read_bases(column):
    """Strip read starts, read ends and indels from a pileup base column.

    Parameters
    ----------
    column : str
        The fifth mpileup column.

    Returns
    -------
    str
        One character per read covering the position.
    """
    kept = []
    i = 0
    while i < len(column):
        char = column[i]
        if char == "^":
            i += 2
        elif char == "$":
            i += 1
        elif char in "+-":
            j = i + 1
            while column[j].isdigit():
                j += 1
            i = j + int(column[i + 1 : j])
        else:
            kept.append(char)
            i += 1
    return "".join(kept)


def parse_line(line):
    """Count the bases of one mpileup line.

    Parameters
    ----------
    line : str
        One mpileup line.

    Returns
    -------
    key : str
        ``<contig>_<position>``.
    site : SiteCounts
        Base counts, warning and haplotype base.

    Raises
    ------
    ValueError
        If the read characters do not add up to the reported depth.
    """
    contig, position, allele, depth, column = line.split("\t")[:5]
    allele = allele.upper()
    bases = read_bases(column)
    upper = bases.upper()
    counts = {base: upper.count(base) for base in BASES}
    counts[allele] = bases.count(".") + bases.count(",")
    accounted = (
        sum(counts.values())
        + sum(bases.count(mark) for mark in DELETION_MARKS)
        + sum(bases.count(mark) for mark in SKIP_MARKS)
    )
    if accounted != int(depth):
        raise ValueError(f"unexpected base counts or symbols in pileup line:\n{line}")
    warning = "."
    if sum(counts.values()) > counts[allele]:
        warning = f"{sum(counts.values()) - counts[allele]}_other_alleles"
    if max(counts.values()) > counts[allele]:
        warning = "hap_allele_not_largest_cnt"
    return f"{contig}_{position}", SiteCounts(counts, warning, allele)


def pileup_counts(paths):
    """Pool base counts from mpileup files by haplotype position.

    Parameters
    ----------
    paths : list of str
        Mpileup files, typically one per sample and haplotype.

    Returns
    -------
    dict of str to SiteCounts
        Pooled counts keyed by ``<contig>_<position>``. Warnings from every
        file are joined with commas. The base comes from the last file.
    """
    pooled = {}
    for path in paths:
        with open(path) as handle:
            for line in handle:
                key, site = parse_line(line.rstrip("\n"))
                if key not in pooled:
                    pooled[key] = site
                    continue
                seen = pooled[key]
                for base in BASES:
                    seen.counts[base] += site.counts[base]
                seen.warning += "," + site.warning
                seen.allele = site.allele
    return pooled
