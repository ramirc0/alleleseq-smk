"""Add gene names and types after the region column of a gene table."""

import argparse


def build_parser():
    """Build the name_genes.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the name_genes.py flags.
    """
    p = argparse.ArgumentParser(description="Join gene names into a region table.")
    p.add_argument("--table", required=True, help="Region table keyed by gene_id.")
    p.add_argument("--names", required=True, help="Names TSV from genes_bed.py.")
    p.add_argument("--output", required=True, help="Named table to write.")
    return p


def main(argv=None):
    """Write the table with ``gene_name`` and ``gene_type`` columns.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with open(args.names) as handle:
        handle.readline()
        names = {}
        for line in handle:
            gene, name, gene_type = line.rstrip("\n").split("\t")
            names[gene] = (name, gene_type)
    with open(args.table) as table, open(args.output, "w") as out:
        header = table.readline().rstrip("\n").split("\t")
        out.write("\t".join([header[0], "gene_name", "gene_type", *header[1:]]) + "\n")
        for line in table:
            row = line.rstrip("\n").split("\t")
            name, gene_type = names.get(row[0], (".", "."))
            out.write("\t".join([row[0], name, gene_type, *row[1:]]) + "\n")


if __name__ == "__main__":
    main()
