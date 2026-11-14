"""Build the gene regions BED from a GENCODE GTF.

Each gene becomes one BED4 line named by its ``gene_id``, which is unique in
GENCODE (PAR copies carry a ``_PAR_Y`` suffix). Gene names repeat, so a
separate table maps IDs to names and types for the final reports.
"""

import argparse
import gzip
import re

ATTRIBUTE = re.compile(r'(\S+) "([^"]*)"')


def build_parser():
    """Build the genes_bed.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the genes_bed.py flags.
    """
    p = argparse.ArgumentParser(description="Gene regions from a GTF.")
    p.add_argument("--gtf", required=True, help="GENCODE GTF, plain or gzipped.")
    p.add_argument("--bed", required=True, help="BED4 of genes to write.")
    p.add_argument("--names", required=True, help="TSV of gene IDs, names and types.")
    return p


def main(argv=None):
    """Write the genes BED and the names table.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    seen = set()
    opener = gzip.open if args.gtf.endswith(".gz") else open
    with (
        opener(args.gtf, "rt") as gtf,
        open(args.bed, "w") as bed,
        open(args.names, "w") as names,
    ):
        names.write("gene_id\tgene_name\tgene_type\n")
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[2] != "gene":
                continue
            attributes = dict(ATTRIBUTE.findall(fields[8]))
            gene = attributes["gene_id"]
            if gene in seen:
                raise SystemExit(f"gene_id {gene} appears twice in {args.gtf}")
            seen.add(gene)
            bed.write(f"{fields[0]}\t{int(fields[3]) - 1}\t{fields[4]}\t{gene}\n")
            name = attributes.get("gene_name", gene)
            names.write(f"{gene}\t{name}\t{attributes.get('gene_type', '.')}\n")


if __name__ == "__main__":
    main()
