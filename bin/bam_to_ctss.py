#!/usr/bin/env python3
"""Extract CAGE transcription start sites (CTSS) from a paired-end CAGE BAM.

The 5' end of each READ1 alignment marks a capped 5' end.  For every primary,
mapped READ1 that passes the MAPQ filter we record its 5' genomic position and
strand.  The soft-clipped 5' sequence (if any) is tallied, because it reveals
what precedes the cap in the read: nothing, the non-templated extra G added by
reverse transcriptase, a leftover linker, or a spliced-leader (SL) sequence.
When the last bases of the SL also match the genome (the acceptor ...AG / ...CAG
resembles the SL's 3' end), STAR aligns them and the read's 5' end lands inside
the acceptor; the CTSS of such SL reads is moved downstream past those bases, to
the trans-splice site.

Outputs (PREFIX = --prefix):
  PREFIX.ctss.bed            chrom, pos0, pos0+1, ., count, strand
  PREFIX.ctss.plus.bedgraph  CTSS counts on the plus strand
  PREFIX.ctss.minus.bedgraph CTSS counts on the minus strand
  PREFIX.softclip_5p.tsv     most frequent READ1 5' soft-clipped sequences
  PREFIX.ctss_stats_mqc.tsv  one-row summary table for MultiQC
"""

import argparse
import collections
import sys

import pysam

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def five_prime(read):
    """Return (pos0, strand, softclipped 5' sequence, first aligned bases), all in read orientation."""
    cigar = [op for op in (read.cigartuples or []) if op[0] != 5]  # hard clips are not in the sequence
    seq = read.query_sequence or ""
    if read.is_reverse:
        clip = cigar[-1][1] if cigar and cigar[-1][0] == 4 else 0
        oriented = revcomp(seq)
        return read.reference_end - 1, "-", oriented[:clip], oriented[clip:clip + 8]
    clip = cigar[0][1] if cigar and cigar[0][0] == 4 else 0
    return read.reference_start, "+", seq[:clip], seq[clip:clip + 8]


def sl_aligned(clip, aligned, sl_sequence, min_overlap, max_aligned=3):
    """Number of spliced-leader 3' bases aligned to the genome, or None if the read has no SL.

    The soft clip must end with (a suffix of) the SL.  The SL's last bases often
    match the genome, because the trans-splice acceptor (...AG, ...CAG) resembles
    the SL's own 3' end (...AACAG); then the clip stops up to `max_aligned` bases
    early and those SL bases start the alignment.
    """
    if not sl_sequence or len(clip) < min_overlap:
        return None
    clip, aligned = clip.upper(), aligned.upper()
    for k in range(max_aligned + 1):
        sl = sl_sequence[:len(sl_sequence) - k]
        n = min(len(clip), len(sl))
        if n >= min_overlap and clip[-n:] == sl[-n:] and aligned[:k] == sl_sequence[len(sl_sequence) - k:]:
            return k
    return None


def shift_downstream(pos, strand, k):
    return pos + k if strand == "+" else pos - k


def is_usable_read1(read, min_mapq):
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return False
    if read.is_qcfail or read.is_duplicate:
        return False
    if read.is_paired and not read.is_read1:
        return False
    return read.mapping_quality >= min_mapq


def write_bedgraph(path, ctss, strand):
    with open(path, "w") as out:
        for (chrom, pos, s), count in sorted(ctss.items()):
            if s == strand:
                out.write(f"{chrom}\t{pos}\t{pos + 1}\t{count}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bam", help="coordinate-sorted BAM file")
    parser.add_argument("--prefix", required=True, help="output file prefix (sample name)")
    parser.add_argument("--min-mapq", type=int, default=10, help="minimum READ1 mapping quality (default: %(default)s)")
    parser.add_argument("--sl-sequence", default="", help="spliced-leader sequence, 5'->3' (default: none)")
    parser.add_argument("--sl-min-overlap", type=int, default=8, help="minimum soft-clip length to call an SL match (default: %(default)s)")
    parser.add_argument("--max-clip-report", type=int, default=30, help="truncate reported soft clips to this many 3'-most bases (default: %(default)s)")
    parser.add_argument("--top", type=int, default=50, help="number of soft-clip sequences to report (default: %(default)s)")
    args = parser.parse_args()
    sl_sequence = args.sl_sequence.upper()

    ctss = collections.Counter()
    clips = collections.Counter()
    n_used = n_noclip = n_g = n_sl = n_sl_shifted = 0

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if not is_usable_read1(read, args.min_mapq):
                continue
            pos, strand, clip, aligned = five_prime(read)
            clip = clip.upper()
            k = sl_aligned(clip, aligned, sl_sequence, args.sl_min_overlap)
            if k is not None:
                n_sl += 1
                if k:
                    # SL bases aligned to the acceptor: the trans-splice site is k bases downstream
                    pos = shift_downstream(pos, strand, k)
                    n_sl_shifted += 1
            ctss[(read.reference_name, pos, strand)] += 1
            n_used += 1
            if not clip:
                n_noclip += 1
            elif clip == "G":
                n_g += 1
            clips[clip[-args.max_clip_report:] if clip else "-"] += 1

    with open(f"{args.prefix}.ctss.bed", "w") as out:
        for (chrom, pos, strand), count in sorted(ctss.items()):
            out.write(f"{chrom}\t{pos}\t{pos + 1}\t.\t{count}\t{strand}\n")
    write_bedgraph(f"{args.prefix}.ctss.plus.bedgraph", ctss, "+")
    write_bedgraph(f"{args.prefix}.ctss.minus.bedgraph", ctss, "-")

    with open(f"{args.prefix}.softclip_5p.tsv", "w") as out:
        out.write("softclip_5p\tlength\tcount\tfraction\n")
        for clip, count in clips.most_common(args.top):
            length = 0 if clip == "-" else len(clip)
            out.write(f"{clip}\t{length}\t{count}\t{count / max(n_used, 1):.4f}\n")

    def pct(x):
        return f"{100 * x / max(n_used, 1):.2f}"

    with open(f"{args.prefix}.ctss_stats_mqc.tsv", "w") as out:
        out.write("# id: 'ctss_stats'\n")
        out.write("# section_name: 'CTSS extraction'\n")
        out.write("# description: 'READ1 5-prime ends used as CAGE TSSs, and what was soft-clipped before them (extra G, linker, spliced leader).'\n")
        out.write("# plot_type: 'table'\n")
        out.write("Sample\tREAD1 used\tDistinct CTSS\tNo 5' clip (%)\tExtra G clip (%)\tOther clip (%)\tSL clip (%)\tSL end aligned, CTSS moved (%)\n")
        other = n_used - n_noclip - n_g
        out.write(f"{args.prefix}\t{n_used}\t{len(ctss)}\t{pct(n_noclip)}\t{pct(n_g)}\t{pct(other)}\t{pct(n_sl)}\t{pct(n_sl_shifted)}\n")

    if n_used == 0:
        print(f"WARNING: no usable READ1 alignments found in {args.bam}", file=sys.stderr)


if __name__ == "__main__":
    main()
