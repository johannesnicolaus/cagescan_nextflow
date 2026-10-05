#!/usr/bin/env python3
"""RECLU-style reproducible TSS clusters (Ohmiya et al. 2014, BMC Genomics 15:269).

  reclu.py hierarchy --prefix S --ctss S.ctss.bed S.paraclu.txt
      From raw (uncut) paraclu output of one replicate, keep clusters with
      >= --min-tpm-per-base TPM per base, and score each cluster by its
      hierarchical stability: its own paraclu stability (max density / min
      density) plus that of every cluster containing it (RECLU).  With
      --score log-stability (default) the log stabilities are summed instead:
      a stability depends on the distance to the nearest stray CTSS and so
      varies by large factors between replicates, which makes the plain sum
      poorly reproducible.  --score tpm uses expression instead.  Singleton
      and top-level clusters have no defined stability and contribute nothing;
      singletons are not reported (as in paraclu-cut).
      -> S.reclu.bed (chrom, start, end, name, hierarchical stability, strand,
         TPM, own stability)

  reclu.py idr --prefix G.S1_S2 S1.reclu.bed S2.reclu.bed
      Pair clusters of two replicates with >= --min-overlap reciprocal
      overlap (1:1, best overlap first; coordinates = the intersection), then
      estimate the irreproducible discovery rate of the pairs from their
      hierarchical stabilities with the IDR model of Li et al. 2011 (python
      idr package; its own peak merging is not used because RECLU clusters
      are nested).  With fewer than --min-pairs pairs no model is fitted and
      IDR is reported as NA (overlap-only reproducibility).
      -> G.S1_S2.reclu_idr.tsv

  reclu.py combine --prefix P *.reclu_idr.tsv
      Keep pairs with IDR < --max-idr (or NA) and length <= --max-length,
      unite all replicate pairs and groups, keep the innermost clusters (those
      not containing another reproducible cluster), and merge any that still
      overlap.  These are the "sharp" TSS clusters.
      -> P.reclu_clusters.bed (paraclu-like: chrom, first, last, name, score,
         strand; first/last are 0-based inclusive CTSS positions)
         P.reclu_stats_mqc.tsv
"""

import argparse
import bisect
import collections
import os

INF = 1e99  # paraclu prints +/-1e100 for undefined densities


# ---------------------------------------------------------------- hierarchy

def total_ctss(path):
    total = 0
    with open(path) as handle:
        for line in handle:
            fields = line.split("\t")
            if len(fields) >= 6 and not line.startswith(("#", "track")):
                total += int(float(fields[4]))
    return total


def read_paraclu(path):
    """Yield (chrom, strand, first, last, n_sites, total, stability) from raw paraclu output."""
    with open(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            first, last = int(f[2]), int(f[3])
            min_d, max_d = float(f[6]), float(f[7])
            stability = max_d / min_d if (0 < min_d and max_d < INF) else 0.0
            yield f[0], f[1], first, last, int(f[4]), float(f[5]), stability


def hierarchy(args):
    import math
    total = total_ctss(args.ctss)
    per_strand = collections.defaultdict(list)
    for chrom, strand, first, last, n_sites, value, stability in read_paraclu(args.paraclu):
        per_strand[(chrom, strand)].append((first, -last, value, stability, n_sites))

    rows = []
    for (chrom, strand), clusters in per_strand.items():
        clusters.sort()  # by start, then longest first: parents precede children
        stack = []  # (last, cumulative stability)
        for first, neg_last, value, stability, n_sites in clusters:
            last = -neg_last
            while stack and stack[-1][0] < first:
                stack.pop()
            parent_sum = stack[-1][1] if stack and stack[-1][0] >= last else 0.0
            own = math.log(stability) if (args.score == "log-stability" and stability > 0) else stability
            cumulative = parent_sum + own
            stack.append((last, cumulative))
            if first == last:
                continue  # singleton
            tpm = value * 1e6 / total if total else 0.0
            if tpm / (last - first + 1) < args.min_tpm_per_base:
                continue
            score = tpm if args.score == "tpm" else cumulative
            rows.append((chrom, first, last + 1, strand, score, tpm, stability))

    rows.sort()
    with open(f"{args.prefix}.reclu.bed", "w") as out:
        for i, (chrom, start, end, strand, cumulative, tpm, stability) in enumerate(rows, 1):
            out.write(f"{chrom}\t{start}\t{end}\t{args.prefix}_{i}\t{cumulative:.6g}\t{strand}\t{tpm:.6g}\t{stability:.6g}\n")


# ---------------------------------------------------------------- idr

def read_reclu_bed(path):
    clusters = collections.defaultdict(list)
    with open(path) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            clusters[(f[0], f[5])].append((int(f[1]), int(f[2]), float(f[4])))
    for value in clusters.values():
        value.sort()
    return clusters


def pair_clusters(rep1, rep2, min_overlap):
    """1:1 pairs with reciprocal overlap >= min_overlap: (chrom, strand, start, end, score1, score2)."""
    candidates = []
    for key, c1 in rep1.items():
        c2 = rep2.get(key, [])
        starts2 = [c[0] for c in c2]
        for i, (s1, e1, v1) in enumerate(c1):
            length1 = e1 - s1
            # Reciprocal overlap >= f implies |start difference| <= (1 - f) * longer length
            # and lengths within a factor f, so the longer one is <= length1 / f.
            slack = (1 - min_overlap) * length1 / min_overlap + 1
            lo = bisect.bisect_left(starts2, s1 - slack)
            hi = bisect.bisect_right(starts2, s1 + slack)
            for j in range(lo, hi):
                s2, e2, v2 = c2[j]
                inter = min(e1, e2) - max(s1, s2)
                if inter <= 0:
                    continue
                overlap = min(inter / length1, inter / (e2 - s2))
                if overlap >= min_overlap:
                    candidates.append((overlap, key, i, j))
    used1, used2, pairs = set(), set(), []
    for overlap, key, i, j in sorted(candidates, key=lambda c: (-c[0], c[1], c[2], c[3])):
        if (key, i) in used1 or (key, j) in used2:
            continue
        used1.add((key, i))
        used2.add((key, j))
        s1, e1, v1 = rep1[key][i]
        s2, e2, v2 = rep2[key][j]
        pairs.append((key[0], key[1], max(s1, s2), min(e1, e2), v1, v2))
    pairs.sort()
    return pairs


def run_idr(pairs, seed):
    import numpy
    from idr.idr import calc_global_IDR, fit_model_and_calc_local_idr

    rng = numpy.random.RandomState(seed)
    s1 = numpy.array([p[4] for p in pairs], dtype=float)
    s2 = numpy.array([p[5] for p in pairs], dtype=float)
    # Same ranking as idr.build_rank_vectors, with a seeded tie breaker.
    r1 = numpy.lexsort((rng.random_sample(len(s1)), s1)).argsort().astype(int)
    r2 = numpy.lexsort((rng.random_sample(len(s2)), s2)).argsort().astype(int)
    local = fit_model_and_calc_local_idr(r1, r2)
    return local, calc_global_IDR(local)


def idr(args):
    pairs = pair_clusters(read_reclu_bed(args.rep1), read_reclu_bed(args.rep2), args.min_overlap)
    if len(pairs) >= args.min_pairs:
        local, glob = run_idr(pairs, args.seed)
        local, glob = [f"{x:.6g}" for x in local], [f"{x:.6g}" for x in glob]
    else:
        local = glob = ["NA"] * len(pairs)
    with open(f"{args.prefix}.reclu_idr.tsv", "w") as out:
        out.write("chrom\tstart\tend\tstrand\tstability_rep1\tstability_rep2\tlocal_idr\tidr\n")
        for (chrom, strand, start, end, v1, v2), l, g in zip(pairs, local, glob):
            out.write(f"{chrom}\t{start}\t{end}\t{strand}\t{v1:.6g}\t{v2:.6g}\t{l}\t{g}\n")


# ---------------------------------------------------------------- combine

def combine(args):
    reproducible = collections.defaultdict(set)
    stats = []
    for path in sorted(args.files):
        n_pairs = n_pass = 0
        no_model = False
        with open(path) as handle:
            handle.readline()
            for line in handle:
                chrom, start, end, strand, _, _, _, value = line.rstrip("\n").split("\t")
                n_pairs += 1
                start, end = int(start), int(end)
                if value == "NA":
                    no_model = True
                elif float(value) >= args.max_idr:
                    continue
                if end - start > args.max_length:
                    continue
                n_pass += 1
                reproducible[(chrom, strand)].add((start, end))
        name = os.path.basename(path).replace(".reclu_idr.tsv", "")
        stats.append((name, n_pairs, n_pass, "no (too few pairs)" if no_model else "yes"))

    rows = []
    for (chrom, strand), spans in reproducible.items():
        spans = sorted(spans, key=lambda se: (se[0], -se[1]))
        # Innermost: spans containing no other span.
        innermost = []
        for i, (s, e) in enumerate(spans):
            contains = False
            for s2, e2 in spans[i + 1:]:
                if s2 >= e:
                    break
                if e2 <= e and (s2, e2) != (s, e):
                    contains = True
                    break
            if not contains:
                innermost.append([s, e])
        innermost.sort()
        merged = []
        for s, e in innermost:
            if merged and s < merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        rows.extend((chrom, s, e, strand) for s, e in merged)

    rows.sort()
    with open(f"{args.prefix}.reclu_clusters.bed", "w") as out:
        for i, (chrom, start, end, strand) in enumerate(rows, 1):
            out.write(f"{chrom}\t{start}\t{end - 1}\treclu_{i}\t0\t{strand}\n")
    with open(f"{args.prefix}.reclu_stats_mqc.tsv", "w") as out:
        out.write("# id: 'reclu_stats'\n")
        out.write("# section_name: 'Sharp TSS clusters (RECLU)'\n")
        out.write(f"# description: 'Replicate pairs: clusters paired by >= 90% reciprocal overlap, and pairs passing IDR < {args.max_idr} and length <= {args.max_length} bp. {len(rows)} sharp clusters in total.'\n")
        out.write("# plot_type: 'table'\n")
        out.write("Replicate pair\tPaired clusters\tReproducible\tIDR model fitted\n")
        for name, n_pairs, n_pass, fitted in stats:
            out.write(f"{name}\t{n_pairs}\t{n_pass}\t{fitted}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("hierarchy", help="filter and score raw paraclu clusters of one replicate")
    p.add_argument("--prefix", required=True)
    p.add_argument("--ctss", required=True, help="CTSS bed of the replicate (for TPM normalisation)")
    p.add_argument("--min-tpm-per-base", type=float, default=0.1, help="(default: %(default)s)")
    p.add_argument("--score", choices=["stability", "log-stability", "tpm"], default="log-stability",
                   help="score used for IDR: hierarchical stability (RECLU), hierarchical log-stability, or TPM (default: %(default)s)")
    p.add_argument("paraclu", help="raw paraclu output (not paraclu-cut)")
    p.set_defaults(func=hierarchy)

    p = sub.add_parser("idr", help="pair two replicates and estimate IDR")
    p.add_argument("--prefix", required=True)
    p.add_argument("--min-overlap", type=float, default=0.9, help="minimum reciprocal overlap (default: %(default)s)")
    p.add_argument("--min-pairs", type=int, default=50, help="minimum pairs to fit the IDR model (default: %(default)s)")
    p.add_argument("--seed", type=int, default=1, help="seed for rank tie breaking (default: %(default)s)")
    p.add_argument("rep1")
    p.add_argument("rep2")
    p.set_defaults(func=idr)

    p = sub.add_parser("combine", help="unite reproducible clusters and keep the innermost")
    p.add_argument("--prefix", required=True)
    p.add_argument("--max-idr", type=float, default=0.1, help="(default: %(default)s)")
    p.add_argument("--max-length", type=int, default=200, help="(default: %(default)s)")
    p.add_argument("files", nargs="+")
    p.set_defaults(func=combine)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
