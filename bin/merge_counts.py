#!/usr/bin/env python3
"""Merge per-sample counts into matrices (rows = features, columns = samples).

  merge_counts.py tss --clusters tss_clusters.tsv --prefix P S1.tss_counts.tsv S2.tss_counts.tsv ...
      -> P.tss_cluster_counts.tsv, P.tss_cluster_cpm.tsv, and
         P.tss_cluster_sl_counts.tsv (spliced-leader READ1s) if any were found

  merge_counts.py transcripts --prefix P S1.transcripts.gtf S2.transcripts.gtf ...
      (StringTie -e output) -> P.transcript_tpm.tsv, P.transcript_cov.tsv

The sample name is the file name up to the first '.'.
"""

import argparse
import collections
import os
import re

ATTR_RE = re.compile(r'\s*([^\s]+)\s+"?([^";]*)"?;?')


def sample_name(path):
    return os.path.basename(path).split(".")[0]


def write_matrix(path, header_cols, rows, samples, values, fmt):
    with open(path, "w") as out:
        out.write("\t".join(header_cols + samples) + "\n")
        for row in rows:
            out.write("\t".join(row + [fmt(values[s].get(row[0], 0)) for s in samples]) + "\n")


def tss(args):
    info = collections.OrderedDict()
    with open(args.clusters) as handle:
        handle.readline()
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            info[fields[0]] = fields[1:4]  # location, strand, peak
    samples, counts, sl_counts = [], {}, {}
    for path in sorted(args.files):
        name = sample_name(path)
        samples.append(name)
        counts[name], sl_counts[name] = {}, {}
        with open(path) as handle:
            handle.readline()
            for line in handle:
                tss_id, count, sl_count = line.rstrip("\n").split("\t")
                counts[name][tss_id] = int(count)
                sl_counts[name][tss_id] = int(sl_count)
    rows = [[tss_id] + values for tss_id, values in info.items()]
    header = ["tss_id", "location", "strand", "peak"]
    write_matrix(f"{args.prefix}.tss_cluster_counts.tsv", header, rows, samples, counts, str)
    if any(any(v.values()) for v in sl_counts.values()):
        write_matrix(f"{args.prefix}.tss_cluster_sl_counts.tsv", header, rows, samples, sl_counts, str)
    cpm = {}
    for name in samples:
        total = sum(counts[name].values()) or 1
        cpm[name] = {k: v * 1e6 / total for k, v in counts[name].items()}
    write_matrix(f"{args.prefix}.tss_cluster_cpm.tsv", header, rows, samples, cpm, lambda x: f"{x:.3f}")


def transcripts(args):
    samples, tpm, cov, info = [], {}, {}, collections.OrderedDict()
    for path in sorted(args.files):
        name = sample_name(path)
        samples.append(name)
        tpm[name], cov[name] = {}, {}
        with open(path) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if line.startswith("#") or len(fields) < 9 or fields[2] != "transcript":
                    continue
                attrs = {m.group(1): m.group(2) for m in ATTR_RE.finditer(fields[8])}
                tid = attrs.get("reference_id") or attrs["transcript_id"]
                gid = attrs.get("ref_gene_id") or attrs.get("gene_id", "")
                info.setdefault(tid, [gid, f"{fields[0]}:{fields[3]}-{fields[4]}", fields[6]])
                tpm[name][tid] = float(attrs.get("TPM", 0))
                cov[name][tid] = float(attrs.get("cov", 0))
    rows = [[tid] + values for tid, values in sorted(info.items())]
    header = ["transcript_id", "tss_id", "location", "strand"]
    write_matrix(f"{args.prefix}.transcript_tpm.tsv", header, rows, samples, tpm, lambda x: f"{x:.4f}")
    write_matrix(f"{args.prefix}.transcript_cov.tsv", header, rows, samples, cov, lambda x: f"{x:.4f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p_tss = sub.add_parser("tss")
    p_tss.add_argument("--clusters", required=True)
    p_tss.add_argument("--prefix", required=True)
    p_tss.add_argument("files", nargs="+")
    p_tss.set_defaults(func=tss)
    p_tx = sub.add_parser("transcripts")
    p_tx.add_argument("--prefix", required=True)
    p_tx.add_argument("files", nargs="+")
    p_tx.set_defaults(func=transcripts)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
