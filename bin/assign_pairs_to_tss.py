#!/usr/bin/env python3
"""Group CAGE read pairs by TSS cluster.

Without UMIs, a CAGEscan "molecule" can only be approximated by the TSS
cluster that its READ1 5' end falls into.  This script keeps the read pairs
whose READ1 5' end lies inside a TSS cluster (extended by --window bases on
either side, same strand), tags every kept record with TC:Z:<tss_id>, and
counts READ1s per cluster.  Paired READ1s must be in a proper pair.

Pass 1 assigns READ1s to clusters; pass 2 writes all primary records of the
assigned templates, so mates are kept even when they map far downstream.

Pass 2 also builds FANTOM5-style CAGEscan clusters (Bertin et al. 2017,
CAGEscan-Clustering): for each TSS cluster, the union of the aligned blocks of
all its pairs, as one BED12 line (score = number of pairs, capped at 1000;
thickStart/thickEnd = the TSS cluster).  Unlike the assembled transcripts,
these are not isoform-resolved.

Outputs (PREFIX = --prefix):
  PREFIX.tss.bam                  TSS-anchored pairs (coordinate-sorted, indexed)
  PREFIX.tss_counts.tsv           tss_id, READ1 count, SL-clipped READ1 count
  PREFIX.cagescan_clusters.bed12  union of pair blocks per TSS cluster
  PREFIX.tss_assign_stats_mqc.tsv one-row summary table for MultiQC
"""

import argparse
import bisect
import collections

import pysam

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def five_prime(read):
    """Return (pos0, strand, softclipped 5' sequence in read orientation)."""
    cigar = read.cigartuples
    seq = read.query_sequence or ""
    if read.is_reverse:
        clip = cigar[-1][1] if cigar and cigar[-1][0] == 4 else 0
        return read.reference_end - 1, "-", revcomp(seq[len(seq) - clip:]) if clip else ""
    clip = cigar[0][1] if cigar and cigar[0][0] == 4 else 0
    return read.reference_start, "+", seq[:clip]


def is_sl(clip, sl_sequence, min_overlap, max_aligned=3):
    """True if the soft clip ends with (a suffix of) the spliced-leader sequence.

    The last few SL bases can align to the genome by chance, so the clip may end
    up to `max_aligned` bases before the SL's own 3' end.
    """
    if not sl_sequence or len(clip) < min_overlap:
        return False
    clip = clip.upper()
    for k in range(max_aligned + 1):
        sl = sl_sequence[:len(sl_sequence) - k]
        n = min(len(clip), len(sl))
        if n >= min_overlap and clip[-n:] == sl[-n:]:
            return True
    return False


def is_primary(read):
    return not (read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail or read.is_duplicate)


class ClusterIndex:
    """Look up the TSS cluster nearest to a position, per chromosome and strand."""

    def __init__(self, bed, window):
        self.window = window
        self.clusters = collections.defaultdict(list)
        self.ids = []
        self.info = {}  # tss_id -> (chrom, start, end, strand)
        with open(bed) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6 or line.startswith(("#", "track")):
                    continue
                chrom, start, end, tss_id, strand = fields[0], int(fields[1]), int(fields[2]), fields[3], fields[5]
                self.clusters[(chrom, strand)].append((start, end, tss_id))
                self.ids.append(tss_id)
                self.info[tss_id] = (chrom, start, end, strand)
        for value in self.clusters.values():
            value.sort()
        self.starts = {key: [c[0] for c in value] for key, value in self.clusters.items()}

    def find(self, chrom, pos, strand):
        key = (chrom, strand)
        clusters = self.clusters.get(key)
        if not clusters:
            return None
        # Candidates: clusters starting at or before pos + window.  paraclu-cut
        # clusters on one strand are disjoint, so ends are sorted like starts and
        # we can stop at the first cluster that ends too far upstream.
        i = bisect.bisect_right(self.starts[key], pos + self.window)
        best, best_dist = None, None
        for j in range(i - 1, -1, -1):
            start, end, tss_id = clusters[j]
            if end + self.window <= pos:
                break
            dist = 0 if start <= pos < end else min(abs(pos - start), abs(pos - (end - 1)))
            if best is None or dist < best_dist:
                best, best_dist = tss_id, dist
        return best


def add_blocks(union, blocks):
    """Record aligned blocks as a sparse map: block start -> furthest block end."""
    for start, end in blocks:
        if end > union.get(start, start):
            union[start] = end


def merge_blocks(union):
    merged = []
    for start in sorted(union):
        end = union[start]
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def write_cagescan_clusters(path, index, unions, counts):
    rgb = {"+": "0,127,0", "-": "127,0,127"}
    rows = []
    for tss_id, union in unions.items():
        chrom, tc_start, tc_end, strand = index.info[tss_id]
        blocks = merge_blocks(union)
        start, end = blocks[0][0], blocks[-1][1]
        thick_start, thick_end = max(start, tc_start), min(end, tc_end)
        if thick_start >= thick_end:
            thick_start = thick_end = start
        rows.append((chrom, start, end, tss_id, min(1000, counts[tss_id]), strand, thick_start, thick_end, blocks))
    with open(path, "w") as out:
        for chrom, start, end, tss_id, score, strand, ts, te, blocks in sorted(rows):
            sizes = ",".join(str(e - s) for s, e in blocks)
            offsets = ",".join(str(s - start) for s, _ in blocks)
            out.write(f"{chrom}\t{start}\t{end}\t{tss_id}\t{score}\t{strand}\t{ts}\t{te}\t{rgb.get(strand, '0,0,0')}\t{len(blocks)}\t{sizes}\t{offsets}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bam", help="coordinate-sorted BAM file")
    parser.add_argument("clusters", help="TSS clusters BED (chrom, start, end, id, score, strand)")
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--window", type=int, default=50, help="extend clusters by this many bases (default: %(default)s)")
    parser.add_argument("--min-mapq", type=int, default=10, help="minimum READ1 mapping quality (default: %(default)s)")
    parser.add_argument("--sl-sequence", default="", help="spliced-leader sequence, 5'->3' (default: none)")
    parser.add_argument("--sl-min-overlap", type=int, default=8)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    sl_sequence = args.sl_sequence.upper()

    index = ClusterIndex(args.clusters, args.window)
    assigned = {}
    counts = collections.Counter()
    sl_counts = collections.Counter()
    n_read1 = 0

    # Pass 1: assign READ1s to TSS clusters.
    with pysam.AlignmentFile(args.bam, "rb", threads=args.threads) as bam:
        for read in bam.fetch(until_eof=True):
            if not is_primary(read) or (read.is_paired and not read.is_read1):
                continue
            if read.mapping_quality < args.min_mapq:
                continue
            if read.is_paired and not read.is_proper_pair:
                continue
            n_read1 += 1
            pos, strand, clip = five_prime(read)
            tss_id = index.find(read.reference_name, pos, strand)
            if tss_id is None:
                continue
            assigned[read.query_name] = tss_id
            counts[tss_id] += 1
            if is_sl(clip, sl_sequence, args.sl_min_overlap):
                sl_counts[tss_id] += 1

    # Pass 2: write every primary record of the assigned templates.
    out_bam = f"{args.prefix}.tss.bam"
    n_records = 0
    unions = collections.defaultdict(dict)
    with pysam.AlignmentFile(args.bam, "rb", threads=args.threads) as bam:
        with pysam.AlignmentFile(out_bam, "wb", template=bam, threads=args.threads) as out:
            for read in bam.fetch(until_eof=True):
                if not is_primary(read):
                    continue
                tss_id = assigned.get(read.query_name)
                if tss_id is None:
                    continue
                read.set_tag("TC", tss_id, value_type="Z")
                out.write(read)
                n_records += 1
                if read.reference_name == index.info[tss_id][0]:
                    add_blocks(unions[tss_id], read.get_blocks())
    pysam.index(out_bam)
    write_cagescan_clusters(f"{args.prefix}.cagescan_clusters.bed12", index, unions, counts)

    with open(f"{args.prefix}.tss_counts.tsv", "w") as out:
        out.write("tss_id\tcount\tsl_count\n")
        for tss_id in index.ids:
            out.write(f"{tss_id}\t{counts[tss_id]}\t{sl_counts[tss_id]}\n")

    n_assigned = len(assigned)
    with open(f"{args.prefix}.tss_assign_stats_mqc.tsv", "w") as out:
        out.write("# id: 'tss_assign_stats'\n")
        out.write("# section_name: 'TSS cluster assignment'\n")
        out.write("# description: 'READ1s whose 5-prime end falls in a consensus TSS cluster; their pairs are used to build transcripts.'\n")
        out.write("# plot_type: 'table'\n")
        out.write("Sample\tREAD1 passing filters\tREAD1 in TSS clusters\tIn TSS clusters (%)\tClusters with reads\tBAM records kept\n")
        frac = 100 * n_assigned / max(n_read1, 1)
        detected = sum(1 for c in counts.values() if c > 0)
        out.write(f"{args.prefix}\t{n_read1}\t{n_assigned}\t{frac:.2f}\t{detected}\t{n_records}\n")


if __name__ == "__main__":
    main()
