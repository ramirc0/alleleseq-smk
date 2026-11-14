"""Check that bcftools consensus rebuilds vcf2diploid's haplotypes.

vcf2diploid orients unphased hets at random, so its output can only be
matched with the same orientations. This script reads each het's
orientation from the golden haplotype FASTAs, writes them as a phased VCF,
runs ``bcftools consensus`` for both haplotypes and compares every contig.

Usage (local only, inside the rule env)::

    python tests/golden/haplotypes.py <golden dir> <scratch dir>
"""

import subprocess
import sys
from pathlib import Path

PRIMARY = [f"chr{c}" for c in [*range(1, 23), "X", "Y", "M"]]


def read_fasta(path):
    """Yield contig names and sequences from a FASTA file.

    Parameters
    ----------
    path : Path
        FASTA file.

    Yields
    ------
    tuple of (str, str)
        Name (first word of the header) and sequence.
    """
    name, chunks = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.rstrip("\n"))
    if name is not None:
        yield name, "".join(chunks)


def golden_orientation(golden, sample):
    """Map each het to the allele vcf2diploid put on hap1.

    Parameters
    ----------
    golden : Path
        Worked-example directory.
    sample : str
        VCF sample ID.

    Returns
    -------
    dict of tuple to str
        Hap1 base, upper case, keyed by (chromosome, 1-based position).
    """
    pgenome = golden / "pgenome_IND1"
    wanted = {}
    with open(pgenome / f"{sample}_hetSNVs_hap1.bed") as bed:
        for line in bed:
            contig, _, end, name = line.rstrip("\n").split("\t")
            wanted.setdefault(contig, []).append((int(end), name))
    hap1_base = {}
    for contig, sequence in read_fasta(pgenome / f"{sample}_hap1.fa"):
        for end, name in wanted.get(contig, []):
            chrom, pos = name.split("_")[:2]
            hap1_base[(chrom, int(pos))] = sequence[end - 1].upper()
    return hap1_base


def write_phased_vcf(vcf, out, hap1_base):
    """Rewrite genotypes so that allele 1 is the golden hap1 allele.

    Parameters
    ----------
    vcf : Path
        Original VCF.
    out : Path
        Phased VCF to write.
    hap1_base : dict of tuple to str
        Golden hap1 allele per het.
    """
    with open(vcf) as source, open(out, "w") as target:
        for line in source:
            if line.startswith("#"):
                target.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            alleles = [fields[3], *fields[4].split(",")]
            keys = fields[8].split(":")
            values = fields[9].split(":")
            gt = values[keys.index("GT")].replace("/", "|").split("|")
            base = hap1_base.get((fields[0], int(fields[1])))
            if (
                base is not None
                and gt[0] != gt[1]
                and alleles[int(gt[0])].upper() != base
            ):
                gt.reverse()
            values[keys.index("GT")] = "|".join(gt)
            fields[9] = ":".join(values)
            target.write("\t".join(fields) + "\n")


def main(golden, scratch, sample="IND1"):
    """Run the check and report per-contig agreement.

    Parameters
    ----------
    golden : str
        Worked-example directory.
    scratch : str
        Directory for intermediate files.
    sample : str
        VCF sample ID.
    """
    golden, scratch = Path(golden), Path(scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    reference = golden / "resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    phased = scratch / "phased.vcf"
    write_phased_vcf(
        golden / f"resources/{sample}.snvs.vcf",
        phased,
        golden_orientation(golden, sample),
    )
    subprocess.run(
        ["bcftools", "view", "-Oz", "-o", f"{phased}.gz", str(phased)], check=True
    )
    subprocess.run(["bcftools", "index", "-f", f"{phased}.gz"], check=True)
    primary = scratch / "primary.fa"
    with open(primary, "w") as out:
        subprocess.run(
            ["samtools", "faidx", str(reference), *PRIMARY], stdout=out, check=True
        )
    failures = 0
    for hap in (1, 2):
        rebuilt = scratch / f"hap{hap}.fa"
        with open(rebuilt, "w") as out:
            subprocess.run(
                [
                    "bcftools",
                    "consensus",
                    "-f",
                    str(primary),
                    "-s",
                    sample,
                    "-H",
                    str(hap),
                    f"{phased}.gz",
                ],
                stdout=out,
                check=True,
            )
        ours = dict(read_fasta(rebuilt))
        theirs = dict(
            read_fasta(golden / f"pgenome_IND1/{sample}_hap{hap}.fa")
        )
        for contig, sequence in theirs.items():
            mine = ours[contig.rsplit("_", 1)[0]]
            status = "same" if mine == sequence else "DIFFERENT"
            if mine != sequence and mine.upper() == sequence.upper():
                status = "same except case"
            if status != "same":
                failures += 1
            print(f"hap{hap}\t{contig}\t{len(sequence)}\t{status}", flush=True)
    print("ALL IDENTICAL" if failures == 0 else f"{failures} contigs differ")


if __name__ == "__main__":
    main(*sys.argv[1:])
