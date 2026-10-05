#!/usr/bin/env python3
"""Find enhancer-RNA-like divergent TSS pairs.

Enhancer RNAs (eRNAs) are transcribed bidirectionally from enhancers: CAGE
shows a pair of TSS clusters on opposite strands, facing away from each other
(divergent: minus-strand TSS upstream of the plus-strand TSS) and typically
100-400 bp apart, away from gene promoters; the RNAs are short and unspliced
and the two strands are similarly expressed (Andersson et al. 2014, Nature
507:455).  This script finds such pairs among the TSS clusters of one cluster
set:

  1. usable clusters: >= --min-reads READ1s (all samples) and a spliced-leader
     fraction <= --max-sl-fraction (SL-dominated clusters are trans-splice
     acceptors, not TSSs)
  2. pairs: minus-strand peak p-, plus-strand peak p+, with
     --min-distance <= p+ - p- <= --max-distance (divergent)
  3. distal: both peaks > --distal-distance from any annotated transcript 5'
     end (only with --annotation)
  4. balanced: directionality |F - R| / (F + R) < --max-directionality
  5. unspliced: no multi-exon transcript starts at either cluster

Convergent pairs (plus-strand peak upstream of a minus-strand peak, same
distances) are counted through the same filters as a chance / overlapping-
transcription control; the divergent/convergent ratio at each step is the
evidence that the candidates are more than chance.

Outputs (PREFIX = --prefix):
  PREFIX.erna_candidates.bed   candidate pairs (filters 1-5): the region from p- to p+
  PREFIX.erna_candidates.tsv   the same with read counts, directionality, distances
  PREFIX.erna_stats_mqc.tsv    divergent vs convergent pairs after each filter (MultiQC)
"""

import argparse
import bisect
import collections
import os
import re


def read_clusters(path):
    """tss_clusters.tsv: tss_id, chrom:start1-end, strand, peak (1-based), total, width."""
    rows = []
    with open(path) as handle:
        handle.readline()
        for line in handle:
            tss_id, location, strand, peak, _, _ = line.rstrip("\n").split("\t")
            rows.append((tss_id, location.rsplit(":", 1)[0], strand, int(peak) - 1))
    return rows


def read_matrix_sums(path):
    """Sum of the per-sample columns of a merge_counts.py matrix (columns 5...)."""
    sums = {}
    if path and os.path.exists(path):
        with open(path) as handle:
            handle.readline()
            for line in handle:
                f = line.rstrip("\n").split("\t")
                sums[f[0]] = sum(float(x) for x in f[4:])
    return sums


def spliced_tss(gtf):
    """TSS cluster IDs (gene_id) with at least one multi-exon transcript."""
    exons = collections.Counter()
    gene = {}
    with open(gtf) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            exons[tid] += 1
            gene[tid] = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
    return {gene[t] for t, n in exons.items() if n > 1}


def annotated_five_prime(gtf):
    """Sorted annotated transcript 5' ends per chromosome (both strands)."""
    tx = collections.defaultdict(list)
    info = {}
    with open(gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            if not m:
                continue
            tx[m.group(1)].append((int(f[3]) - 1, int(f[4])))
            info[m.group(1)] = (f[0], f[6])
    five = collections.defaultdict(list)
    for tid, ex in tx.items():
        chrom, strand = info[tid]
        five[chrom].append(min(s for s, _ in ex) if strand == "+" else max(e for _, e in ex) - 1)
    for v in five.values():
        v.sort()
    return five


def distance_to(five, chrom, pos):
    f = five.get(chrom, [])
    i = bisect.bisect_left(f, pos)
    return min([abs(pos - f[j]) for j in (i - 1, i) if 0 <= j < len(f)] or [float("inf")])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--clusters", required=True, help="tss_clusters.tsv")
    parser.add_argument("--counts", required=True, help="tss_cluster_counts.tsv (READ1s per cluster and sample)")
    parser.add_argument("--sl-counts", default=None, help="tss_cluster_sl_counts.tsv (optional)")
    parser.add_argument("--transcripts", required=True, help="merged TSS-anchored transcripts GTF of the same cluster set")
    parser.add_argument("--annotation", default=None, help="reference annotation GTF (optional; enables the distal filter)")
    parser.add_argument("--min-reads", type=float, default=10, help="(default: %(default)s)")
    parser.add_argument("--max-sl-fraction", type=float, default=0.1, help="(default: %(default)s)")
    parser.add_argument("--min-distance", type=int, default=50, help="minimum peak-to-peak distance (default: %(default)s)")
    parser.add_argument("--max-distance", type=int, default=400, help="maximum peak-to-peak distance (default: %(default)s)")
    parser.add_argument("--distal-distance", type=int, default=500, help="minimum distance to annotated 5' ends (default: %(default)s)")
    parser.add_argument("--max-directionality", type=float, default=0.8, help="(default: %(default)s)")
    args = parser.parse_args()

    counts = read_matrix_sums(args.counts)
    sl = read_matrix_sums(args.sl_counts)
    spliced = spliced_tss(args.transcripts)
    five = annotated_five_prime(args.annotation) if args.annotation else None

    usable = collections.defaultdict(list)
    for tss_id, chrom, strand, peak in read_clusters(args.clusters):
        n = counts.get(tss_id, 0)
        if n >= args.min_reads and sl.get(tss_id, 0) / n <= args.max_sl_fraction:
            usable[chrom].append((peak, strand, tss_id, n))
    for v in usable.values():
        v.sort()

    pairs = {"divergent": [], "convergent": []}
    for chrom, rows in usable.items():
        peaks = [r[0] for r in rows]
        for i, (p, s, tid, n) in enumerate(rows):
            j = bisect.bisect_left(peaks, p + args.min_distance)
            while j < len(rows) and rows[j][0] - p <= args.max_distance:
                q, s2, tid2, n2 = rows[j]
                if s == "-" and s2 == "+":
                    pairs["divergent"].append((chrom, p, q, tid, tid2, n, n2))
                elif s == "+" and s2 == "-":
                    pairs["convergent"].append((chrom, p, q, tid, tid2, n, n2))
                j += 1

    def distal(x):
        return five is None or min(distance_to(five, x[0], x[1]), distance_to(five, x[0], x[2])) > args.distal_distance

    def directionality(x):
        return abs(x[6] - x[5]) / (x[5] + x[6])

    def balanced(x):
        return directionality(x) < args.max_directionality

    def unspliced(x):
        return x[3] not in spliced and x[4] not in spliced

    steps = [("opposite-strand pairs", lambda x: True)]
    if five is not None:
        steps.append((f"distal (> {args.distal_distance} bp from annotated 5' ends)", distal))
    steps += [(f"balanced (directionality < {args.max_directionality})", balanced), ("unspliced", unspliced)]

    kept = {k: list(v) for k, v in pairs.items()}
    table = []
    for label, test in steps:
        for k in kept:
            kept[k] = [x for x in kept[k] if test(x)]
        d, c = len(kept["divergent"]), len(kept["convergent"])
        table.append((label, d, c, f"{d / c:.1f}" if c else "n/a"))

    cand = sorted(kept["divergent"])
    with open(f"{args.prefix}.erna_candidates.bed", "w") as bed, open(f"{args.prefix}.erna_candidates.tsv", "w") as tsv:
        tsv.write("chrom\tstart\tend\ttss_minus\ttss_plus\treads_minus\treads_plus\tdirectionality\tpeak_distance\tdistance_to_annotated_5p\n")
        for x in cand:
            name = f"{x[3]}|{x[4]}"
            dist = min(distance_to(five, x[0], x[1]), distance_to(five, x[0], x[2])) if five is not None else "NA"
            bed.write(f"{x[0]}\t{x[1]}\t{x[2] + 1}\t{name}\t{min(1000, int(x[5] + x[6]))}\t.\n")
            tsv.write(f"{x[0]}\t{x[1]}\t{x[2] + 1}\t{x[3]}\t{x[4]}\t{x[5]:g}\t{x[6]:g}\t{directionality(x):.3f}\t{x[2] - x[1]}\t{dist}\n")

    with open(f"{args.prefix}.erna_stats_mqc.tsv", "w") as out:
        out.write(f"# id: 'erna_stats_{args.prefix.replace('.', '_')}'\n")
        out.write(f"# section_name: 'eRNA-like divergent TSS pairs ({args.prefix})'\n")
        out.write(f"# description: 'Pairs of opposite-strand non-SL TSS clusters ({args.min_distance}-{args.max_distance} bp, >= {args.min_reads:g} reads each) after each filter. Divergent pairs are eRNA candidates; convergent pairs are the chance control.'\n")
        out.write("# plot_type: 'table'\n")
        out.write("Filter\tDivergent pairs\tConvergent pairs (control)\tRatio\n")
        for label, d, c, r in table:
            out.write(f"{label}\t{d}\t{c}\t{r}\n")


if __name__ == "__main__":
    main()
