#!/usr/bin/env python3
"""Helpers for building consensus TSS clusters from CTSS files.

  tss_clusters.py pool --prefix P a.ctss.bed b.ctss.bed ...
      Sum CTSS counts across samples -> P.ctss.bed (same 6-column format).

  tss_clusters.py annotate --prefix P --ctss pooled.ctss.bed paraclu.bed
      Turn paraclu output into named, correctly bounded TSS clusters, with
      sparse tails trimmed (--trim-fraction):
      P.tss_clusters.bed (BED9: chrom, start, end, id, score, strand,
      peak, peak+1, rgb) and P.tss_clusters.tsv (id, location, strand,
      peak position, total CTSS count, width).

paraclu reports the first and last CTSS position of each cluster; since the
CTSS positions are 0-based, the BED end coordinate is last + 1.
"""

import argparse
import bisect
import collections


def read_ctss(path):
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6 or line.startswith(("#", "track")):
                continue
            yield fields[0], int(fields[1]), fields[5], int(float(fields[4]))


def pool(args):
    counts = collections.Counter()
    for path in args.ctss:
        for chrom, pos, strand, count in read_ctss(path):
            counts[(chrom, pos, strand)] += count
    with open(f"{args.prefix}.ctss.bed", "w") as out:
        for (chrom, pos, strand), count in sorted(counts.items()):
            out.write(f"{chrom}\t{pos}\t{pos + 1}\t.\t{count}\t{strand}\n")


def trim_tails(members, fraction):
    """Drop sparse CTSS from both cluster ends, removing at most `fraction` of the counts per end.

    paraclu-cut keeps the largest qualifying cluster, which can include isolated
    background CTSS tens of bases away from the core.  The dominant CTSS is never removed.
    """
    if not members or fraction <= 0:
        return members
    total = sum(c for _, c in members)
    peak = max(range(len(members)), key=lambda i: members[i][1])
    lo, dropped = 0, 0
    while lo < peak and dropped + members[lo][1] <= fraction * total:
        dropped += members[lo][1]
        lo += 1
    hi, dropped = len(members), 0
    while hi - 1 > peak and dropped + members[hi - 1][1] <= fraction * total:
        dropped += members[hi - 1][1]
        hi -= 1
    return members[lo:hi]


def annotate(args):
    # Index pooled CTSS per (chrom, strand) for peak lookup.
    positions = collections.defaultdict(list)
    for chrom, pos, strand, count in read_ctss(args.ctss):
        positions[(chrom, strand)].append((pos, count))
    for value in positions.values():
        value.sort()
    starts = {key: [p for p, _ in value] for key, value in positions.items()}

    clusters = []
    with open(args.paraclu) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            chrom, first, last, strand = fields[0], int(fields[1]), int(fields[2]), fields[5]
            key = (chrom, strand)
            lo = bisect.bisect_left(starts.get(key, []), first)
            hi = bisect.bisect_right(starts.get(key, []), last)
            members = trim_tails(positions.get(key, [])[lo:hi], args.trim_fraction)
            if not members:
                continue
            first, last = members[0][0], members[-1][0]
            total = sum(c for _, c in members)
            # Dominant CTSS; ties broken towards the 5'-most position.
            best = max(members, key=lambda pc: (pc[1], -pc[0] if strand == "+" else pc[0]))
            clusters.append((chrom, first, last + 1, strand, best[0], total))

    clusters.sort()
    width = max(6, len(str(len(clusters))))
    rgb = {"+": "0,127,0", "-": "127,0,127"}
    with open(f"{args.prefix}.tss_clusters.bed", "w") as bed, open(f"{args.prefix}.tss_clusters.tsv", "w") as tsv:
        tsv.write("tss_id\tlocation\tstrand\tpeak\ttotal_count\twidth\n")
        for i, (chrom, start, end, strand, peak, total) in enumerate(clusters, 1):
            tss_id = f"{args.id_prefix}{i:0{width}d}"
            score = min(1000, total)
            bed.write(f"{chrom}\t{start}\t{end}\t{tss_id}\t{score}\t{strand}\t{peak}\t{peak + 1}\t{rgb.get(strand, '0,0,0')}\n")
            tsv.write(f"{tss_id}\t{chrom}:{start + 1}-{end}\t{strand}\t{peak + 1}\t{total}\t{end - start}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p_pool = sub.add_parser("pool", help="sum CTSS counts across samples")
    p_pool.add_argument("--prefix", required=True)
    p_pool.add_argument("ctss", nargs="+")
    p_pool.set_defaults(func=pool)

    p_annot = sub.add_parser("annotate", help="name paraclu clusters and find their peaks")
    p_annot.add_argument("--prefix", required=True)
    p_annot.add_argument("--ctss", required=True, help="pooled CTSS bed")
    p_annot.add_argument("--id-prefix", default="TC", help="cluster ID prefix (default: %(default)s)")
    p_annot.add_argument("--trim-fraction", type=float, default=0.01,
                         help="trim cluster ends by up to this fraction of the cluster's counts per end (default: %(default)s)")
    p_annot.add_argument("paraclu", help="paraclu bed (chrom, first, last, name, count, strand)")
    p_annot.set_defaults(func=annotate)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
