#!/usr/bin/env python3
"""Merge per-sample CAGEscan clusters into meta-clusters (FANTOM5 style).

  merge_cagescan_clusters.py --prefix P S1.cagescan_clusters.bed12 S2.cagescan_clusters.bed12 ...
      -> P.cagescan_clusters.bed12

Clusters with the same name (TSS cluster ID) are combined: their blocks are
united, the thick region spans all input thick regions, and the score is the
number of samples (libraries) contributing to the meta-cluster, as in
Bertin et al. 2017 (Sci Data 4:170147).
"""

import argparse
import collections


def read_bed12(path):
    with open(path) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or line.startswith(("#", "track")):
                continue
            start = int(f[1])
            sizes = [int(x) for x in f[10].rstrip(",").split(",")]
            offsets = [int(x) for x in f[11].rstrip(",").split(",")]
            blocks = [(start + o, start + o + s) for o, s in zip(offsets, sizes)]
            yield f[3], f[0], f[5], int(f[6]), int(f[7]), blocks


def merge_blocks(blocks):
    merged = []
    for start, end in sorted(blocks):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("files", nargs="+")
    args = parser.parse_args()

    meta = collections.OrderedDict()
    for path in sorted(args.files):
        for name, chrom, strand, thick_start, thick_end, blocks in read_bed12(path):
            m = meta.setdefault(name, {"chrom": chrom, "strand": strand, "blocks": [], "thick": [], "samples": 0})
            m["blocks"].extend(blocks)
            if thick_start < thick_end:
                m["thick"].append((thick_start, thick_end))
            m["samples"] += 1

    rgb = {"+": "0,127,0", "-": "127,0,127"}
    rows = []
    for name, m in meta.items():
        blocks = merge_blocks(m["blocks"])
        start, end = blocks[0][0], blocks[-1][1]
        if m["thick"]:
            thick_start, thick_end = min(t[0] for t in m["thick"]), max(t[1] for t in m["thick"])
        else:
            thick_start = thick_end = start
        rows.append((m["chrom"], start, end, name, min(1000, m["samples"]), m["strand"], thick_start, thick_end, blocks))

    with open(f"{args.prefix}.cagescan_clusters.bed12", "w") as out:
        for chrom, start, end, name, score, strand, ts, te, blocks in sorted(rows):
            sizes = ",".join(str(e - s) for s, e in blocks)
            offsets = ",".join(str(s - start) for s, _ in blocks)
            out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{ts}\t{te}\t{rgb.get(strand, '0,0,0')}\t{len(blocks)}\t{sizes}\t{offsets}\n")


if __name__ == "__main__":
    main()
