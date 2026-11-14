"""Tag a haplotype's contig names, such as ``chr1`` to ``chr1_hap1``.

Reads a FASTA or a chain file from ``bcftools consensus`` on stdin and
writes it to stdout. FASTA headers lose their description. In a chain, only
the query (haplotype) name changes.
"""

import argparse
import sys


def build_parser():
    """Build the rename_contigs.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the rename_contigs.py flags.
    """
    p = argparse.ArgumentParser(description="Add a haplotype suffix to contig names.")
    p.add_argument(
        "--format", choices=["fasta", "chain"], required=True, help="Input type."
    )
    p.add_argument("--suffix", required=True, help="Suffix to append, e.g. _hap1.")
    p.add_argument("--drop", nargs="*", default=[], help="FASTA contigs to leave out.")
    return p


def main(argv=None):
    """Copy stdin to stdout with renamed contigs.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    rename = rename_fasta if args.format == "fasta" else rename_chain
    sys.stdout.writelines(rename(sys.stdin, args.suffix, set(args.drop)))


def rename_fasta(lines, suffix, drop):
    """Rename FASTA headers and skip dropped contigs.

    Parameters
    ----------
    lines : iterable of str
        FASTA lines.
    suffix : str
        Suffix for each name.
    drop : set of str
        Contig names to leave out.

    Yields
    ------
    str
        Output lines.
    """
    keep = True
    for line in lines:
        if line.startswith(">"):
            name = line[1:].split()[0]
            keep = name not in drop
            line = f">{name}{suffix}\n"
        if keep:
            yield line


def rename_chain(lines, suffix, drop):
    """Rename the query contig of each chain header.

    Parameters
    ----------
    lines : iterable of str
        Chain file lines.
    suffix : str
        Suffix for each query name.
    drop : set of str
        Unused. Chains keep every contig, so lifted coordinates exist even for
        contigs left out of the FASTA.

    Yields
    ------
    str
        Output lines.
    """
    for line in lines:
        if line.startswith("chain"):
            fields = line.split()
            fields[7] += suffix
            line = " ".join(fields) + "\n"
        yield line


if __name__ == "__main__":
    main()
