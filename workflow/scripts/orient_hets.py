"""Phase every genotype so that ``bcftools consensus -H 1/-H 2`` can split them.

Genotypes come from the callset. A het takes its phase from a statistically
phased VCF, such as SHAPEIT5's, when that VCF calls the same het with a
phasing confidence (``PP``) at or above the cutoff, or with none. Common
variants carry no ``PP``. Other phased genotypes without a phase set
(``PS``) keep their order. A phase set,
such as GATK's physical phasing, only fixes phase within its block, so each
block gets one random orientation. Each unphased het gets a random
orientation from a seeded generator, as the AlleleSeq paper does for hets the
parents cannot phase. The orientations and the phased fraction are written
out, because gene-level aggregation is only meaningful for phased hets.
"""

import argparse
import gzip
import io
import random


def build_parser():
    """Build the orient_hets.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the orient_hets.py flags.
    """
    p = argparse.ArgumentParser(description="Orient unphased hets at random.")
    p.add_argument("--vcf", required=True, help="Single-sample VCF, plain or gzipped.")
    p.add_argument("--phased", help="Single-sample phased VCF of the sample's hets.")
    p.add_argument(
        "--min-phase-confidence",
        type=float,
        default=0.99,
        help="Smallest PP at which a het from --phased is trusted.",
    )
    p.add_argument("--seed", type=int, default=0, help="Seed for the orientations.")
    p.add_argument(
        "--output", required=True, help="Oriented VCF to write (plain text)."
    )
    p.add_argument(
        "--orientations", required=True, help="Gzipped TSV of oriented hets."
    )
    p.add_argument("--summary", required=True, help="TSV with the phased fraction.")
    return p


def main(argv=None):
    """Write the oriented VCF, the orientation table and the summary.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    phase = (
        confident_phase(args.phased, args.min_phase_confidence) if args.phased else {}
    )
    rng = random.Random(args.seed)
    flips = {}
    hets = phased = 0
    with (
        open_vcf(args.vcf) as vcf,
        open(args.output, "w") as out,
        open_deterministic(args.orientations) as table,
    ):
        table.write("chrom\tpos\tref\talt\tgenotype\toriented\n")
        for line in vcf:
            if line.startswith("#"):
                out.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            genotype, oriented, chromosome_wide = orient(fields, rng, flips, phase)
            alleles = oriented.split("|")
            if len(alleles) == 2 and alleles[0] != alleles[1]:
                hets += 1
                if chromosome_wide:
                    phased += 1
                else:
                    row = [*fields[:2], fields[3], fields[4], genotype, oriented]
                    table.write("\t".join(row) + "\n")
            out.write("\t".join(fields) + "\n")
    with open(args.summary, "w") as out:
        fraction = phased / hets if hets else float("nan")
        out.write("het_sites\tphased\tfraction_phased\n")
        out.write(f"{hets}\t{phased}\t{fraction:.4f}\n")


def open_vcf(path):
    """Open a plain or gzipped VCF as text.

    Parameters
    ----------
    path : str
        VCF path. A ``.gz`` suffix means gzip or BGZF.

    Returns
    -------
    io.TextIOBase
        Text stream over the VCF.
    """
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def confident_phase(path, min_confidence):
    """Read the confidently phased hets of a single-sample phased VCF.

    Parameters
    ----------
    path : str
        Phased VCF, plain or gzipped.
    min_confidence : float
        Smallest ``PP`` kept. Hets without ``PP`` are always kept.

    Returns
    -------
    dict
        Phased genotype, ``"0|1"`` or ``"1|0"``, keyed by
        ``(chrom, pos, ref, alt)``.
    """
    phase = {}
    with open_vcf(path) as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            calls = dict(zip(fields[8].split(":"), fields[9].split(":"), strict=False))
            confidence = calls.get("PP", ".")
            if calls["GT"] in ("0|1", "1|0") and (
                confidence == "." or float(confidence) >= min_confidence
            ):
                phase[(fields[0], fields[1], fields[3], fields[4])] = calls["GT"]
    return phase


def open_deterministic(path):
    """Open a gzipped text file without a timestamp, so reruns are byte-identical.

    Parameters
    ----------
    path : str
        Output path.

    Returns
    -------
    io.TextIOWrapper
        Text stream over a gzip file with mtime 0.
    """
    return io.TextIOWrapper(gzip.GzipFile(path, "wb", mtime=0))


def orient(fields, rng, flips, phase):
    """Phase the genotype of one VCF record in place.

    Parameters
    ----------
    fields : list of str
        VCF record columns. The sample is the tenth column.
    rng : random.Random
        Source of orientations for unphased hets and phase sets.
    flips : dict
        Orientation drawn for each ``(chrom, PS)`` phase set so far. Updated
        in place.
    phase : dict
        Confident phase from a phased VCF, as `confident_phase` returns.

    Returns
    -------
    genotype : str
        The genotype as read.
    oriented : str
        The phased genotype written back into `fields`.
    chromosome_wide : bool
        Whether the genotype was phased along the whole chromosome.
    """
    keys = fields[8].split(":")
    values = fields[9].split(":")
    calls = dict(zip(keys, values, strict=False))
    genotype = calls["GT"]
    alleles = genotype.replace("/", "|").split("|")
    het = len(alleles) == 2 and alleles[0] != alleles[1]
    phase_set = calls.get("PS", ".")
    chromosome_wide = "|" in genotype and phase_set == "."
    phased = phase.get((fields[0], fields[1], fields[3], fields[4])) if het else None
    if phased:
        alleles = phased.split("|")
        chromosome_wide = True
    elif het and not chromosome_wide:
        if "|" in genotype:
            block = (fields[0], phase_set)
            if block not in flips:
                flips[block] = rng.random() < 0.5
            flip = flips[block]
        else:
            flip = rng.random() < 0.5
        if flip:
            alleles.reverse()
    oriented = "|".join(alleles)
    values[keys.index("GT")] = oriented
    fields[9] = ":".join(values)
    return genotype, oriented, chromosome_wide


if __name__ == "__main__":
    main()
