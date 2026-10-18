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
the trans-splice site.  Read pairs whose name starts with 'SL:' had the SL found
and trimmed before alignment (cutadapt); their READ1 starts at the trans-splice site.

Outputs (PREFIX = --prefix):
  PREFIX.ctss.bed            chrom, pos0, pos0+1, ., count, strand
  PREFIX.ctss.plus.bedgraph  CTSS counts on the plus strand
  PREFIX.ctss.minus.bedgraph CTSS counts on the minus strand
  PREFIX.sl_ctss.bed         CTSS of spliced-leader READ1s only (with --sl-sequence)
  PREFIX.nosl_ctss.bed       CTSS of the other READ1s (with --sl-sequence)
  PREFIX.softclip_5p.tsv     most frequent READ1 5' soft-clipped sequences
  PREFIX.ctss_stats_mqc.tsv  one-row summary table for MultiQC
"""

import argparse
import collections
import re
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

SL_TRIMMED = ("SL:", "SLG:")  # read-name prefixes of pairs whose SL was trimmed before alignment (SPLIT_SPLICED_LEADER)
MAX_SL_CLIP = 3  # longest 5' soft clip of an SL-trimmed READ1 that is added back to its 5' end


def sl_trimmed_5p(read, pos, strand, clip):
    """5' end of an SL-trimmed READ1, with a short 5' soft clip added back.

    After the SL (and the extra G before it) is trimmed, the read starts at the trans-splice site and holds no
    non-templated bases; STAR soft-clips first bases that differ from the genome (sequencing errors, strain
    variants), which would otherwise move the CTSS downstream.
    """
    if read.query_name.startswith(SL_TRIMMED) and 0 < len(clip) <= MAX_SL_CLIP:
        return max(0, pos - len(clip)) if strand == "+" else pos + len(clip)
    return pos


def sl_trimmed_unplaceable(read, clip):
    """SL-trimmed READ1 with a longer 5' soft clip: its 5' end cannot be placed (extra non-genomic bases after the
    leader, or an indel near the 5' end), so it gives no TSS."""
    return read.query_name.startswith(SL_TRIMMED) and len(clip) > MAX_SL_CLIP

MD_FIRST = re.compile(r"^0[ACGTN]")
MD_LAST = re.compile(r"[ACGTN]0$")


def first_g_mismatch(read):
    """READ1 aligned from its first 5' base, which is a G that differs from the genome (MD tag): the non-templated
    G added by the reverse transcriptase at the cap; the TSS is the next base (as CAGEr's removeFirstG)."""
    if read.query_name.startswith(SL_TRIMMED) or not read.has_tag("MD"):
        return False
    cigar = [op for op in (read.cigartuples or []) if op[0] != 5]
    seq = read.query_sequence or ""
    if not cigar or not seq:
        return False
    md = read.get_tag("MD")
    if read.is_reverse:
        return cigar[-1][0] in (0, 7, 8) and seq[-1] == "C" and bool(MD_LAST.search(md))
    return cigar[0][0] in (0, 7, 8) and seq[0] == "G" and bool(MD_FIRST.match(md))


def sl_status(read, clip, aligned, sl_sequence, min_overlap):
    """SL bases aligned to the genome (0-3) if READ1 carries the spliced leader, else None.

    A READ1 whose SL was trimmed before alignment starts at the trans-splice site.
    Otherwise the SL is recognised in the soft clip (see sl_aligned).
    """
    if read.query_name.startswith(SL_TRIMMED):
        return 0
    return sl_aligned(clip, aligned, sl_sequence, min_overlap)


def is_usable_read1(read, min_mapq):
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return False
    if read.is_qcfail or read.is_duplicate:
        return False
    if read.is_paired and not read.is_read1:
        return False
    return read.mapping_quality >= min_mapq


def write_ctss(path, ctss):
    with open(path, "w") as out:
        for (chrom, pos, strand), count in sorted(ctss.items()):
            out.write(f"{chrom}\t{pos}\t{pos + 1}\t.\t{count}\t{strand}\n")


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
    parser.add_argument("--remove-first-g", action="store_true",
                        help="READ1s aligned from their first base: a first G that differs from the genome is the non-templated G; the CTSS is the next base")
    parser.add_argument("--max-clip-report", type=int, default=30, help="truncate reported soft clips to this many 3'-most bases (default: %(default)s)")
    parser.add_argument("--top", type=int, default=50, help="number of soft-clip sequences to report (default: %(default)s)")
    args = parser.parse_args()
    sl_sequence = args.sl_sequence.upper()

    ctss = collections.Counter()
    sl_ctss = collections.Counter()
    clips = collections.Counter()
    n_used = n_noclip = n_g = n_g_removed = n_sl = n_sl_trimmed = n_sl_clip_added = n_sl_unplaced = n_sl_shifted = 0

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if not is_usable_read1(read, args.min_mapq):
                continue
            pos, strand, clip, aligned = five_prime(read)
            clip = clip.upper()
            if sl_trimmed_unplaceable(read, clip):
                n_sl_unplaced += 1
                continue
            k = sl_status(read, clip, aligned, sl_sequence, args.sl_min_overlap)
            if k is not None:
                n_sl += 1
                if read.query_name.startswith(SL_TRIMMED):
                    n_sl_trimmed += 1
                    if clip:
                        n_sl_clip_added += 1
                    pos = sl_trimmed_5p(read, pos, strand, clip)
                if k:
                    # SL bases aligned to the acceptor: the trans-splice site is k bases downstream
                    pos = shift_downstream(pos, strand, k)
                    n_sl_shifted += 1
            elif args.remove_first_g and first_g_mismatch(read):
                pos = shift_downstream(pos, strand, 1)
                n_g_removed += 1
            ctss[(read.reference_name, pos, strand)] += 1
            if k is not None:
                sl_ctss[(read.reference_name, pos, strand)] += 1
            n_used += 1
            if not clip:
                n_noclip += 1
            elif clip == "G":
                n_g += 1
            clips[clip[-args.max_clip_report:] if clip else "-"] += 1

    write_ctss(f"{args.prefix}.ctss.bed", ctss)
    if sl_sequence:
        # separate SL and non-SL CTSS, so that trans-splice sites and TSSs are clustered separately
        write_ctss(f"{args.prefix}.sl_ctss.bed", sl_ctss)
        write_ctss(f"{args.prefix}.nosl_ctss.bed", ctss - sl_ctss)
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
        out.write("Sample\tREAD1 used\tDistinct CTSS\tNo 5' clip (%)\tExtra G clip (%)\tExtra G removed, aligned (%)\tOther clip (%)\tSL (%)\tSL trimmed before alignment (%)\tSL trimmed, 5' clip added back (%)\tSL trimmed, 5' end unplaceable (%)\tSL end aligned, CTSS moved (%)\n")
        other = n_used - n_noclip - n_g
        out.write(f"{args.prefix}\t{n_used}\t{len(ctss)}\t{pct(n_noclip)}\t{pct(n_g)}\t{pct(n_g_removed)}\t{pct(other)}\t{pct(n_sl)}\t{pct(n_sl_trimmed)}\t{pct(n_sl_clip_added)}\t{100 * n_sl_unplaced / max(n_used + n_sl_unplaced, 1):.2f}\t{pct(n_sl_shifted)}\n")

    if n_used == 0:
        print(f"WARNING: no usable READ1 alignments found in {args.bam}", file=sys.stderr)


if __name__ == "__main__":
    main()
