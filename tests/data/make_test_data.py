#!/usr/bin/env python3
"""Simulate a tiny paired-end CAGE (CAGEscan-like) data set with a known truth.

Writes, into the output directory:
  genome.fa            two random chromosomes with GT-AG introns
  annotation.gtf       reference gene models
  truth_tss.bed        true TSS of every simulated isoform
  samplesheet.csv      two samples (paths relative to the pipeline base dir)
  SAMPLE_R{1,2}.fastq.gz

Besides six hand-designed genes (geneA-geneF, the truth checked by the tests)
it adds --random-genes single-TSS genes with log-normal expression on chr3, so
that replicate-level statistics (RECLU / IDR) have enough clusters to work on.
The two samples are replicates of one group.

Read structure: READ1 starts at the capped 5' end (a few bases of TSS jitter),
with an extra non-templated G in half of the reads; one gene is trans-spliced
to a spliced leader (SL) that precedes its READ1s.  READ2 is the reverse
complement of the 3' end of a random-length fragment.  A few percent of pairs
are random genomic background.
"""

import argparse
import gzip
import os
import random

SL = "GGTTTAATTACCCAAGTTTGAG"  # C. elegans SL1
READ_LEN = 75
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


# Gene models: chrom, gene, strand, isoforms {name: (exons as 1-based inclusive (start, end) in genomic order, weight)}, SL
GENES = [
    ("chr1", "geneA", "+", {"geneA.1": ([(5001, 5300), (6001, 6400), (7001, 8200)], 1.0)}, False),
    ("chr1", "geneB", "-", {
        "geneB.1": ([(20001, 21200), (22001, 22150), (23001, 23250), (24001, 24300)], 0.6),
        "geneB.2": ([(20001, 21200), (23001, 23250), (24001, 24300)], 0.4),  # skips exon 22001-22150
    }, False),
    ("chr1", "geneC", "+", {
        "geneC.1": ([(40001, 40250), (42001, 42300), (43001, 44500)], 0.5),
        "geneC.2": ([(41001, 41200), (42001, 42300), (43001, 44500)], 0.5),  # alternative first exon / TSS
    }, False),
    ("chr1", "geneD", "-", {"geneD.1": ([(60001, 62000)], 0.7)}, False),
    ("chr2", "geneE", "+", {"geneE.1": ([(10001, 10400), (11001, 12500)], 0.8)}, True),  # SL trans-spliced
    ("chr2", "geneF", "-", {"geneF.1": ([(30001, 31500), (32001, 32200), (33001, 33300)], 0.6)}, False),
]
CHROM_LEN = {"chr1": 80000, "chr2": 50000}
TSS_JITTER = [(-2, 0.05), (-1, 0.15), (0, 0.6), (1, 0.15), (2, 0.05)]


def add_random_genes(rng, n, slot=6000):
    """Append n random single-TSS genes (1-4 exons) on chr3, one per `slot` bases."""
    if n <= 0:
        return
    CHROM_LEN["chr3"] = n * slot + slot
    for i in range(n):
        base = i * slot + rng.randint(500, 1500)
        exons, pos = [], base
        n_exons = rng.randint(1, 4)
        for k in range(n_exons):
            length = rng.randint(800, 1500) if k == n_exons - 1 else rng.randint(100, 400)
            exons.append((pos, pos + length - 1))
            pos += length + rng.randint(300, 900)
        strand = rng.choice("+-")
        weight = rng.lognormvariate(-2.5, 1.2)
        GENES.append(("chr3", f"rnd{i + 1:03d}", strand, {f"rnd{i + 1:03d}.1": (exons, weight)}, False))


def make_genome(rng):
    genome = {c: [rng.choice("AACCGTTAGT") for _ in range(n)] for c, n in CHROM_LEN.items()}
    # Put canonical splice sites into every intron (GT...AG on the transcript strand).
    for chrom, _, strand, isoforms, _ in GENES:
        for exons, _ in isoforms.values():
            for (_, e1), (s2, _) in zip(exons, exons[1:]):
                left, right = ("GT", "AG") if strand == "+" else ("CT", "AC")
                seq = genome[chrom]
                seq[e1:e1 + 2] = list(left)  # first two intron bases (0-based e1, e1+1)
                seq[s2 - 3:s2 - 1] = list(right)  # last two intron bases
    return {c: "".join(s) for c, s in genome.items()}


def transcript_seq(genome, chrom, strand, exons):
    seq = "".join(genome[chrom][s - 1:e] for s, e in exons)
    return seq if strand == "+" else revcomp(seq)


def mutate(seq, rng, rate=0.002):
    return "".join(rng.choice("ACGT".replace(b, "")) if rng.random() < rate else b for b in seq)


def jitter(rng):
    r, acc = rng.random(), 0.0
    for offset, p in TSS_JITTER:
        acc += p
        if r < acc:
            return offset
    return 0


def simulate_sample(name, genome, n_pairs, rng, outdir, expression_scale):
    isoforms = []
    for chrom, gene, strand, isos, sl in GENES:
        for iso, (exons, weight) in isos.items():
            isoforms.append((chrom, strand, exons, weight * expression_scale.get(gene, 1.0), sl))
    total = sum(i[3] for i in isoforms)
    r1 = gzip.open(os.path.join(outdir, f"{name}_R1.fastq.gz"), "wt", compresslevel=9)
    r2 = gzip.open(os.path.join(outdir, f"{name}_R2.fastq.gz"), "wt", compresslevel=9)
    qual = "I" * READ_LEN
    for n in range(n_pairs):
        if rng.random() < 0.04:  # genomic background
            chrom = rng.choice(list(CHROM_LEN))
            start = rng.randrange(0, CHROM_LEN[chrom] - 600)
            frag = genome[chrom][start:start + rng.randint(200, 500)]
            if rng.random() < 0.5:
                frag = revcomp(frag)
            prefix = ""
        else:
            x, acc = rng.random() * total, 0.0
            for chrom, strand, exons, weight, sl in isoforms:
                acc += weight
                if x < acc:
                    break
            tx = transcript_seq(genome, chrom, strand, exons)
            # TSS jitter: negative offsets extend into the genomic upstream sequence.
            off = jitter(rng)
            if off < 0:
                if strand == "+":
                    up = genome[chrom][exons[0][0] - 1 + off:exons[0][0] - 1]
                else:
                    up = revcomp(genome[chrom][exons[-1][1]:exons[-1][1] - off])
                tx = up + tx
            else:
                tx = tx[off:]
            length = max(120, min(int(rng.gauss(380, 120)), len(tx), 900))
            frag = tx[:length]
            prefix = SL if sl else ("G" if rng.random() < 0.5 else "")
        read1 = mutate((prefix + frag)[:READ_LEN], rng)
        read2 = mutate(revcomp(frag)[:READ_LEN], rng)
        rid = f"{name}.{n + 1}"
        r1.write(f"@{rid}/1\n{read1}\n+\n{qual[:len(read1)]}\n")
        r2.write(f"@{rid}/2\n{read2}\n+\n{qual[:len(read2)]}\n")
    r1.close()
    r2.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default=os.path.dirname(os.path.abspath(__file__)))
    parser.add_argument("--pairs", type=int, default=12000)
    parser.add_argument("--random-genes", type=int, default=150)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    add_random_genes(rng, args.random_genes)
    genome = make_genome(rng)
    with open(os.path.join(args.outdir, "genome.fa"), "w") as out:
        for chrom, seq in genome.items():
            out.write(f">{chrom}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")

    with open(os.path.join(args.outdir, "annotation.gtf"), "w") as gtf, open(os.path.join(args.outdir, "truth_tss.bed"), "w") as bed:
        for chrom, gene, strand, isoforms, _ in GENES:
            for iso, (exons, _) in isoforms.items():
                attrs = f'gene_id "{gene}"; transcript_id "{iso}";'
                gtf.write(f"{chrom}\tsim\ttranscript\t{exons[0][0]}\t{exons[-1][1]}\t.\t{strand}\t.\t{attrs}\n")
                for s, e in exons:
                    gtf.write(f"{chrom}\tsim\texon\t{s}\t{e}\t.\t{strand}\t.\t{attrs}\n")
                tss0 = exons[0][0] - 1 if strand == "+" else exons[-1][1] - 1
                bed.write(f"{chrom}\t{tss0}\t{tss0 + 1}\t{iso}\t0\t{strand}\n")

    simulate_sample("SAMPLE1", genome, args.pairs, rng, args.outdir, {})
    simulate_sample("SAMPLE2", genome, args.pairs, rng, args.outdir, {"geneA": 0.5, "geneF": 2.0})

    with open(os.path.join(args.outdir, "samplesheet.csv"), "w") as out:
        out.write("sample,fastq_1,fastq_2,group\n")
        for name in ("SAMPLE1", "SAMPLE2"):
            out.write(f"{name},tests/data/{name}_R1.fastq.gz,tests/data/{name}_R2.fastq.gz,A\n")


if __name__ == "__main__":
    main()
