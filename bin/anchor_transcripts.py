#!/usr/bin/env python3
"""Anchor assembled transcripts to CAGE TSS clusters.

Each transcript's 5' end is matched to a TSS cluster on the same strand
(cluster extended by --window bases).  Anchored transcripts get their 5' end
moved to the cluster's dominant CTSS, are renamed <tss_id>.<n> (n ordered by
coverage) and grouped under gene_id <tss_id>, so a TSS with several isoforms
yields several transcripts.  Transcripts that become identical after snapping
are collapsed.  Unanchored transcripts are dropped unless --keep-unanchored.

Outputs (PREFIX = --prefix):
  PREFIX.transcripts.gtf          anchored transcripts (GTF)
  PREFIX.transcripts.bed12        same, BED12; the 5' block is thick for
                                  TSS-anchored transcripts (CAGEscan convention)
  PREFIX.anchor_stats_mqc.tsv     one-row summary table for MultiQC
"""

import argparse
import bisect
import collections
import re

ATTR_RE = re.compile(r'\s*([^\s]+)\s+"?([^";]*)"?;?')
COPY_ATTRS = ("cov", "FPKM", "TPM", "ref_gene_id", "ref_gene_name", "reference_id", "class_code")


def parse_attrs(text):
    return {m.group(1): m.group(2) for m in ATTR_RE.finditer(text)}


def read_gtf(path):
    transcripts = collections.OrderedDict()
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in ("transcript", "exon"):
                continue
            attrs = parse_attrs(fields[8])
            tid = attrs.get("transcript_id")
            if tid is None:
                continue
            t = transcripts.setdefault(tid, {"chrom": fields[0], "strand": fields[6], "exons": [], "attrs": {}})
            if fields[2] == "transcript":
                t["attrs"].update(attrs)
            else:
                t["exons"].append((int(fields[3]) - 1, int(fields[4])))  # 0-based half-open
                for key, value in attrs.items():
                    t["attrs"].setdefault(key, value)
    for t in transcripts.values():
        t["exons"] = merge_exons(t["exons"])
    return [dict(t, id=tid) for tid, t in transcripts.items() if t["exons"]]


def merge_exons(exons):
    merged = []
    for start, end in sorted(exons):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def read_clusters(path):
    """Read tss_clusters.tsv: tss_id, chrom:start1-end, strand, peak (1-based), total_count, width."""
    clusters = collections.defaultdict(list)
    with open(path) as handle:
        header = handle.readline()
        assert header.startswith("tss_id"), f"unexpected header in {path}"
        for line in handle:
            tss_id, location, strand, peak, total, _ = line.rstrip("\n").split("\t")
            chrom, span = location.rsplit(":", 1)
            start1, end = span.split("-")
            clusters[(chrom, strand)].append((int(start1) - 1, int(end), tss_id, int(peak) - 1, int(total)))
    for value in clusters.values():
        value.sort()
    return clusters


def find_cluster(clusters, starts, chrom, pos, strand, window):
    key = (chrom, strand)
    if key not in clusters:
        return None
    i = bisect.bisect_right(starts[key], pos + window)
    best, best_dist = None, None
    for j in range(i - 1, -1, -1):
        cluster = clusters[key][j]
        start, end = cluster[0], cluster[1]
        if end + window <= pos:
            break
        dist = 0 if start <= pos < end else min(abs(pos - start), abs(pos - (end - 1)))
        if best is None or dist < best_dist:
            best, best_dist = cluster, dist
    return best


def five_prime(t):
    return t["exons"][0][0] if t["strand"] == "+" else t["exons"][-1][1] - 1


def snap(t, peak):
    """Move the 5' end to peak if that keeps the first exon non-empty; return the shift (bases added upstream)."""
    exons = list(t["exons"])
    if t["strand"] == "+":
        start, end = exons[0]
        if peak >= end:
            return 0
        exons[0] = (peak, end)
        shift = start - peak
    else:
        start, end = exons[-1]
        if peak + 1 <= start:
            return 0
        exons[-1] = (start, peak + 1)
        shift = peak + 1 - end
    t["exons"] = exons
    return shift


def coverage(t):
    for key in ("cov", "TPM", "FPKM"):
        try:
            return float(t["attrs"][key])
        except (KeyError, ValueError):
            continue
    return 0.0


def gtf_attrs(pairs):
    return " ".join(f'{k} "{v}";' for k, v in pairs if v is not None and v != "")


def write_outputs(prefix, records):
    rgb = {"+": "0,127,0", "-": "127,0,127"}
    records.sort(key=lambda t: (t["chrom"], t["exons"][0][0], t["exons"][-1][1], t["new_id"]))
    with open(f"{prefix}.transcripts.gtf", "w") as gtf, open(f"{prefix}.transcripts.bed12", "w") as bed:
        for t in records:
            chrom, strand, exons = t["chrom"], t["strand"], t["exons"]
            start, end = exons[0][0], exons[-1][1]
            base = [("gene_id", t["gene_id"]), ("transcript_id", t["new_id"])]
            extra = [
                ("tss_id", t.get("tss_id")),
                ("tss_peak", t.get("tss_peak")),
                ("tss_count", t.get("tss_count")),
                ("tss_shift", t.get("tss_shift")),
                ("anchored", "yes" if t.get("tss_id") else "no"),
            ] + [(k, t["attrs"].get(k)) for k in COPY_ATTRS]
            gtf.write(f"{chrom}\tcagescan\ttranscript\t{start + 1}\t{end}\t.\t{strand}\t.\t{gtf_attrs(base + extra)}\n")
            ordered = exons if strand == "+" else exons[::-1]
            for n, (s, e) in enumerate(ordered, 1):
                gtf.write(f"{chrom}\tcagescan\texon\t{s + 1}\t{e}\t.\t{strand}\t.\t{gtf_attrs(base + [('exon_number', n)])}\n")

            head = exons[0] if strand == "+" else exons[-1]
            thick_start, thick_end = head if t.get("tss_id") else (start, start)
            score = min(1000, int(t.get("tss_count") or 0))
            sizes = ",".join(str(e - s) for s, e in exons)
            offsets = ",".join(str(s - start) for s, _ in exons)
            bed.write(
                f"{chrom}\t{start}\t{end}\t{t['new_id']}\t{score}\t{strand}\t{thick_start}\t{thick_end}\t"
                f"{rgb.get(strand, '0,0,0')}\t{len(exons)}\t{sizes}\t{offsets}\n"
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("gtf", help="assembled transcripts (StringTie GTF)")
    parser.add_argument("clusters", help="tss_clusters.tsv from tss_clusters.py annotate")
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--window", type=int, default=50, help="max distance between transcript 5' end and cluster (default: %(default)s)")
    parser.add_argument("--keep-unanchored", action="store_true", help="also output transcripts without a TSS cluster")
    args = parser.parse_args()

    clusters = read_clusters(args.clusters)
    starts = {key: [c[0] for c in value] for key, value in clusters.items()}
    transcripts = read_gtf(args.gtf)

    anchored, unanchored = [], []
    for t in transcripts:
        hit = find_cluster(clusters, starts, t["chrom"], five_prime(t), t["strand"], args.window)
        if hit is None:
            unanchored.append(t)
            continue
        _, _, tss_id, peak, total = hit
        t["tss_shift"] = snap(t, peak)
        t["tss_id"] = tss_id
        t["tss_peak"] = f"{t['chrom']}:{peak + 1}:{t['strand']}"
        t["tss_count"] = total
        anchored.append(t)

    # Collapse transcripts that became identical after snapping (keep the best covered).
    unique = {}
    for t in sorted(anchored, key=coverage, reverse=True):
        unique.setdefault((t["tss_id"], t["chrom"], t["strand"], tuple(t["exons"])), t)
    n_collapsed = len(anchored) - len(unique)

    by_tss = collections.defaultdict(list)
    for t in unique.values():
        by_tss[t["tss_id"]].append(t)
    records = []
    for tss_id, group in by_tss.items():
        group.sort(key=lambda t: (-coverage(t), t["exons"]))
        for n, t in enumerate(group, 1):
            t["gene_id"] = tss_id
            t["new_id"] = f"{tss_id}.{n}"
            records.append(t)
    if args.keep_unanchored:
        for t in unanchored:
            t["gene_id"] = "U_" + t["attrs"].get("gene_id", t["id"])
            t["new_id"] = "U_" + t["id"]
            records.append(t)

    write_outputs(args.prefix, records)

    lengths = sorted(sum(e - s for s, e in t["exons"]) for t in records)
    median = lengths[len(lengths) // 2] if lengths else 0
    with open(f"{args.prefix}.anchor_stats_mqc.tsv", "w") as out:
        out.write("# id: 'anchor_stats'\n")
        out.write("# section_name: 'Transcript anchoring'\n")
        out.write("# description: 'Assembled transcripts whose 5-prime end falls in a TSS cluster; their 5-prime end is moved to the dominant CTSS.'\n")
        out.write("# plot_type: 'table'\n")
        out.write("Sample\tAssembled\tAnchored\tCollapsed duplicates\tUnanchored\tUnanchored kept\tTSS clusters with transcripts\tTranscripts output\tMedian length\n")
        kept = len(unanchored) if args.keep_unanchored else 0
        out.write(f"{args.prefix}\t{len(transcripts)}\t{len(anchored)}\t{n_collapsed}\t{len(unanchored)}\t{kept}\t{len(by_tss)}\t{len(records)}\t{median}\n")


if __name__ == "__main__":
    main()
