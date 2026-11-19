"""Generate the synthetic inputs of the integration test.

The toy genome has all 25 primary contigs (genes sit on chr1 to chr3) and one
unplaced contig. A two-sample joint VCF carries the individual's het and
homozygous SNVs plus a few indels. A WGS alignment covers the primary
contigs, with a doubled stretch on chr2 that the CNV filter must remove. A
SHAPEIT5-style phased VCF phases the unphased gene hets, two thirds of them
confidently, and adds two records the callset lacks. Two RNA samples read both
haplotypes of each gene. Eight genes, listed in
``planted.tsv``, express hap1 at 85%.

Run from the repository root inside the rule env::

    python .tests/make_toy_data.py .tests/integration/resources
"""

import gzip
import subprocess
import sys
from pathlib import Path

import numpy as np

RNG = np.random.default_rng(7)
BASES = np.array(list("ACGT"))
CONTIG_LENGTHS = {"chr1": 90_000, "chr2": 90_000, "chr3": 90_000}
CONTIG_LENGTHS |= {f"chr{c}": 6_000 for c in [*range(4, 23), "X", "Y"]}
CONTIG_LENGTHS |= {"chrM": 4_000, "chrUn_toy": 3_000}
GENES_PER_CONTIG = 14
PLANTED = 8
READ_LENGTH = 75
FRAGMENT = 220
ERROR_RATE = 0.002
SAMPLE = "IND1"


def random_sequence(length):
    """Return a random DNA string."""
    return "".join(RNG.choice(BASES, length))


def make_genes(genome):
    """Lay out two- or three-exon genes on chr1 to chr3."""
    genes = []
    for contig in ("chr1", "chr2", "chr3"):
        step = CONTIG_LENGTHS[contig] // GENES_PER_CONTIG
        for i in range(GENES_PER_CONTIG):
            start = i * step + 1500
            n_exons = 2 + i % 2
            exons, pos = [], start
            for _ in range(n_exons):
                length = int(RNG.integers(500, 800))
                exons.append((pos, pos + length - 1))
                pos += length + int(RNG.integers(300, 900))
            strand = "+" if i % 3 else "-"
            genes.append(
                {
                    "id": f"ENSGTOY{len(genes):05d}.1",
                    "name": f"TOY{len(genes)}",
                    "contig": contig,
                    "strand": strand,
                    "exons": exons,
                }
            )
    mito = {
        "id": "ENSGTOY99999.1",
        "name": "MT-TOY",
        "contig": "chrM",
        "strand": "+",
        "exons": [(200, 1400)],
    }
    return [*genes, mito]


def write_gtf(genes, path):
    """Write the genes as a GENCODE-style GTF."""
    with open(path, "w") as out:
        for gene in genes:
            start, end = gene["exons"][0][0], gene["exons"][-1][1]
            attrs = f'gene_id "{gene["id"]}"; gene_type "protein_coding"; gene_name "{gene["name"]}";'
            tx = f'{attrs} transcript_id "{gene["id"].replace("G", "T", 1)}";'
            out.write(
                f"{gene['contig']}\ttoy\tgene\t{start}\t{end}\t.\t{gene['strand']}\t.\t{attrs}\n"
            )
            out.write(
                f"{gene['contig']}\ttoy\ttranscript\t{start}\t{end}\t.\t{gene['strand']}\t.\t{tx}\n"
            )
            out.writelines(
                f"{gene['contig']}\ttoy\texon\t{a}\t{b}\t.\t{gene['strand']}\t.\t{tx} exon_number {n};\n"
                for n, (a, b) in enumerate(gene["exons"], 1)
            )


def make_variants(genome, genes):
    """Place het and homozygous variants, mostly in exons."""
    variants = {}
    planted_ids = {g["id"] for g in genes[:PLANTED]}
    for gene in genes:
        phased_gene = gene["id"] in planted_ids
        for a, b in gene["exons"]:
            pos = a + int(RNG.integers(20, 60))
            while pos < b - 20:
                ref = genome[gene["contig"]][pos - 1]
                alt = str(RNG.choice([x for x in "ACGT" if x != ref]))
                kind = RNG.random()
                if kind < 0.1:
                    gt = "1/1"
                else:
                    phased = phased_gene or RNG.random() < 0.3
                    order = "0|1" if RNG.random() < 0.5 else "1|0"
                    gt = order if phased else "0/1"
                variants[(gene["contig"], pos)] = (ref, alt, gt)
                # Every fifth het gets a neighbour within a read length.
                pos += (
                    int(RNG.integers(8, 30))
                    if RNG.random() < 0.2
                    else int(RNG.integers(60, 140))
                )
    for contig in ("chr4", "chrX"):
        for pos in range(500, CONTIG_LENGTHS[contig] - 500, 400):
            ref = genome[contig][pos - 1]
            alt = str(RNG.choice([x for x in "ACGT" if x != ref]))
            variants[(contig, pos)] = (ref, alt, "0/1")
    # Indels in gene-free flanks shift haplotype coordinates downstream of them.
    for contig in ("chr1", "chr2"):
        for pos in (700, 89_000):
            ref = genome[contig][pos - 1 : pos + 3]
            variants[(contig, pos)] = (ref, ref[0], "0|1")
            ins = genome[contig][pos + 99]
            variants[(contig, pos + 100)] = (ins, ins + "GATT", "1|0")
    return variants


def haplotypes(genome, variants):
    """Apply the true phase: allele 1 of each genotype goes to haplotype A."""
    haps = {"A": {}, "B": {}}
    for contig, seq in genome.items():
        for hap, index in (("A", 0), ("B", 1)):
            pieces, last = [], 0
            for (c, pos), (ref, alt, gt) in sorted(variants.items()):
                if c != contig:
                    continue
                allele = [ref, alt][int(gt.replace("/", "|").split("|")[index])]
                pieces.append(seq[last : pos - 1])
                pieces.append(allele)
                last = pos - 1 + len(ref)
            pieces.append(seq[last:])
            haps[hap][contig] = "".join(pieces)
    return haps


def lift(variants, contig, pos, hap_index):
    """Map a reference position onto a haplotype through the indels."""
    shift = 0
    for (c, p), (ref, alt, gt) in variants.items():
        if c == contig and p < pos and len(ref) != len(alt):
            allele = [ref, alt][int(gt.replace("/", "|").split("|")[hap_index])]
            shift += len(allele) - len(ref)
    return pos + shift


def transcript(hap_genome, variants, gene, hap_index):
    """Return a gene's spliced sequence on one haplotype."""
    seq = hap_genome[gene["contig"]]
    return "".join(
        seq[
            lift(variants, gene["contig"], a, hap_index) - 1 : lift(
                variants, gene["contig"], b, hap_index
            )
        ]
        for a, b in gene["exons"]
    )


def mutate(read):
    """Add sequencing errors."""
    chars = list(read)
    for i in np.flatnonzero(RNG.random(len(chars)) < ERROR_RATE):
        chars[i] = str(RNG.choice([x for x in "ACGT" if x != chars[i]]))
    return "".join(chars)


def revcomp(seq):
    """Reverse-complement a DNA string."""
    return seq[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))


def write_rna(prefix, haps, variants, genes, pairs_per_gene, planted_ids):
    """Write paired FASTQs sampled from both haplotypes of every gene."""
    with (
        gzip.open(f"{prefix}_R1.fastq.gz", "wt") as r1,
        gzip.open(f"{prefix}_R2.fastq.gz", "wt") as r2,
    ):
        n = 0
        for gene in genes:
            share = 0.85 if gene["id"] in planted_ids else 0.5
            tx = {h: transcript(haps[h], variants, gene, i) for i, h in enumerate("AB")}
            for _ in range(pairs_per_gene):
                hap = "A" if RNG.random() < share else "B"
                seq = tx[hap]
                start = int(RNG.integers(0, len(seq) - FRAGMENT))
                fragment = seq[start : start + FRAGMENT]
                if RNG.random() < 0.5:
                    fragment = revcomp(fragment)
                quality = "I" * READ_LENGTH
                r1.write(f"@r{n}/1\n{mutate(fragment[:READ_LENGTH])}\n+\n{quality}\n")
                r2.write(
                    f"@r{n}/2\n{mutate(revcomp(fragment)[:READ_LENGTH])}\n+\n{quality}\n"
                )
                n += 1


def write_wgs(path, genome, reference):
    """Write a coordinate-sorted WGS CRAM with a doubled stretch on chr2."""
    sam = Path(f"{path}.sam")
    with open(sam, "w") as out:
        out.writelines(
            f"@SQ\tSN:{contig}\tLN:{len(seq)}\n" for contig, seq in genome.items()
        )
        n = 0
        for contig in genome:
            if contig == "chrUn_toy":
                continue
            seq = genome[contig]
            depth = 20
            for start in range(0, len(seq) - READ_LENGTH, max(1, READ_LENGTH // depth)):
                copies = 2 if contig == "chr2" and 30_000 <= start < 40_000 else 1
                for _ in range(copies):
                    out.write(
                        f"w{n}\t0\t{contig}\t{start + 1}\t60\t{READ_LENGTH}M\t*\t0\t0\t{seq[start : start + READ_LENGTH]}\t*\n"
                    )
                    n += 1
    subprocess.run(
        [
            "samtools",
            "sort",
            "-O",
            "cram",
            "--reference",
            str(reference),
            "-o",
            str(path),
            str(sam),
        ],
        check=True,
    )
    subprocess.run(["samtools", "index", str(path)], check=True)
    sam.unlink()


def write_vcf(path, variants, genome):
    """Write a bgzipped two-sample joint VCF."""
    plain = Path(str(path).removesuffix(".gz"))
    with open(plain, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.writelines(
            f"##contig=<ID={contig},length={len(seq)}>\n"
            for contig, seq in genome.items()
        )
        out.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        out.write(
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tOTHER\t{SAMPLE}\n"
        )
        out.writelines(
            f"{contig}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t0/0\t{gt}\n"
            for (contig, pos), (ref, alt, gt) in sorted(
                variants.items(),
                key=lambda kv: (list(genome).index(kv[0][0]), kv[0][1]),
            )
        )
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(path), str(plain)], check=True)
    subprocess.run(["bcftools", "index", "-t", str(path)], check=True)
    plain.unlink()


def write_phased_vcf(path, variants, genome):
    """Write a bgzipped phased VCF of the unphased gene hets, with PP.

    Unphased hets are 0|1 in the truth. A third of them carry no PP, as
    SHAPEIT5's common variants do, a third PP 0.995 and a third PP 0.6 with the
    wrong order. One record is an imputed het at a site the callset lacks, one
    a split multiallelic with another ALT. Neither may reach the haplotypes.
    """
    unphased = [
        (contig, pos, ref, alt)
        for (contig, pos), (ref, alt, gt) in variants.items()
        if gt == "0/1" and contig in ("chr1", "chr2", "chr3")
    ]
    records = []
    for i, (contig, pos, ref, alt) in enumerate(sorted(unphased)):
        gt, pp = [("0|1", "."), ("0|1", "0.995"), ("1|0", "0.6")][i % 3]
        records.append((contig, pos, ref, alt, gt, pp))
    contig, pos, ref, alt = unphased[0]
    other = next(b for b in "ACGT" if b not in (ref, alt))
    records.append((contig, pos, ref, other, "1|0", "."))
    free = next(p for p in range(1000, 2000) if ("chr3", p) not in variants)
    ref = genome["chr3"][free - 1]
    alt = next(b for b in "ACGT" if b != ref)
    records.append(("chr3", free, ref, alt, "0|1", "."))
    plain = Path(str(path).removesuffix(".gz"))
    with open(plain, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.writelines(
            f"##contig=<ID={contig},length={len(seq)}>\n"
            for contig, seq in genome.items()
        )
        out.write(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotypes">\n'
        )
        out.write(
            '##FORMAT=<ID=PP,Number=1,Type=Float,Description="Phasing confidence">\n'
        )
        out.write(
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tOTHER\t{SAMPLE}\n"
        )
        out.writelines(
            f"{c}\t{p}\t.\t{r}\t{a}\t.\t.\t.\tGT:PP\t0|0:.\t{gt}:{pp}\n"
            for c, p, r, a, gt, pp in sorted(
                records, key=lambda x: (list(genome).index(x[0]), x[1], x[3])
            )
        )
    subprocess.run(["bcftools", "view", "-Oz", "-o", str(path), str(plain)], check=True)
    plain.unlink()


def main(outdir):
    """Write every toy input under `outdir`."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    genome = {c: random_sequence(n) for c, n in CONTIG_LENGTHS.items()}
    reference = out / "toy.fa"
    with open(reference, "w") as fa:
        for contig, seq in genome.items():
            fa.write(f">{contig}\n")
            fa.writelines(seq[i : i + 60] + "\n" for i in range(0, len(seq), 60))
    subprocess.run(["samtools", "faidx", str(reference)], check=True)
    genes = make_genes(genome)
    write_gtf(genes, out / "toy.gtf")
    variants = make_variants(genome, genes)
    write_vcf(out / "joint.vcf.gz", variants, genome)
    write_phased_vcf(out / "phased.vcf.gz", variants, genome)
    haps = haplotypes(genome, variants)
    write_wgs(out / "wgs.cram", genome, reference)
    planted_ids = {g["id"] for g in genes[:PLANTED]}
    for sample in ("RNA1", "RNA2"):
        write_rna(str(out / sample), haps, variants, genes, 700, planted_ids)
    with open(out / "planted.tsv", "w") as planted:
        planted.write("gene_id\n" + "".join(f"{g}\n" for g in sorted(planted_ids)))


if __name__ == "__main__":
    main(sys.argv[1])
