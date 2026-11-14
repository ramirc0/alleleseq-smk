"""List a sample's PASS het SNVs as a BED in the reference's naming.

Port of AlleleSeq2's ``get_hetSNV_bed.py``. The name column is
``<chr>_<pos>_<ref>_<allele1>_<allele2>``, with the alleles in genotype
order, so allele 1 sits on hap1 once the genotype is phased.
"""

import argparse
import gzip

from chromosomes import harmonize, primary_names


def build_parser():
    """Build the hetsnv_bed.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the hetsnv_bed.py flags.
    """
    p = argparse.ArgumentParser(description="Write het SNVs as BED.")
    p.add_argument("--vcf", required=True, help="VCF, plain or gzipped.")
    p.add_argument("--sample", required=True, help="Sample column to read.")
    p.add_argument("--fai", required=True, help="FASTA index of the reference.")
    p.add_argument("--output", required=True, help="BED to write.")
    return p


def main(argv=None):
    """Write the het SNV BED.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with open(args.fai) as fai:
        reference = primary_names(line.split("\t")[0] for line in fai)
    if not reference:
        raise SystemExit(f"no primary chromosomes in {args.fai}")
    opener = gzip.open if args.vcf.endswith(".gz") else open
    with opener(args.vcf, "rt") as vcf, open(args.output, "w") as out:
        sample_column = None
        for line in vcf:
            if line.startswith("#CHROM"):
                sample_column = line.rstrip("\n").split("\t").index(args.sample)
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            record = het_snv(fields, sample_column, reference)
            if record is not None:
                out.write("\t".join(record) + "\n")


def het_snv(fields, sample_column, reference):
    """Return the BED fields of a PASS het SNV, or None.

    Parameters
    ----------
    fields : list of str
        VCF record columns.
    sample_column : int
        Index of the sample column.
    reference : list of str
        Primary chromosome names of the reference.

    Returns
    -------
    list of str or None
        Chromosome, start, end and name, or None when the record is not a
        PASS het SNV on a primary chromosome.
    """
    chrom, pos, _, ref, alt, _, status = fields[:7]
    if status != "PASS":
        return None
    genotype = fields[sample_column].split(":")[0].replace("/", "|").split("|")
    if len(genotype) != 2 or "." in genotype:
        return None
    alleles = [ref, *alt.split(",")]
    first, second = alleles[int(genotype[0])], alleles[int(genotype[1])]
    if first == second or not len(first) == len(second) == len(ref) == 1:
        return None
    chrom = harmonize(chrom, reference)
    if chrom is None:
        return None
    return [chrom, str(int(pos) - 1), pos, f"{chrom}_{pos}_{ref}_{first}_{second}"]


if __name__ == "__main__":
    main()
