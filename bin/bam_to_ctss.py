#!/usr/bin/env python3
"""Extract CAGE transcription start sites (CTSS) from a CAGE BAM (paired-end or single-end).

The 5' end of each READ1 (or single read) alignment marks a capped 5' end.  For every primary, mapped READ1 that
passes the MAPQ filter we record its 5' genomic position and strand:
- READ1 aligned from its first base (STAR --alignEndsType Extend5pOfRead1, the pipeline default): a first G that
  differs from the genome (MD tag) is the G added by the reverse transcriptase at the cap; with --remove-first-g the
  CTSS is the next base (as CAGEr's removeFirstG).
- Reads whose spliced leader (SL) was found and trimmed before alignment are named 'SL.<leader>:<name>' or
  'SLG.<leader>:<name>' ('SL:' / 'SLG:' in earlier versions); their READ1 starts at the trans-splice site. Up to 3
  soft-clipped 5' bases are added back; reads with longer 5' clips give no CTSS.
- With local alignment, the soft-clipped 5' sequence reveals what precedes the cap (nothing, the extra G, a linker
  or an SL) and is tallied. An SL is recognised when the clip ends with it; when its last bases also match the
  genome (the acceptor ...AG / ...CAG resembles the SL's 3' end), STAR aligns them and the CTSS is moved downstream
  past those bases, to the trans-splice site.

Outputs (PREFIX = --prefix):
  PREFIX.ctss.bed            chrom, pos0, pos0+1, ., count, strand
  PREFIX.ctss.plus.bedgraph  CTSS counts on the plus strand
  PREFIX.ctss.minus.bedgraph CTSS counts on the minus strand
  PREFIX.sl_ctss.bed         CTSS of spliced-leader READ1s only (with --sl-sequence)
  PREFIX.nosl_ctss.bed       CTSS of the other READ1s (with --sl-sequence)
  PREFIX.sl_leaders_mqc.tsv  share of SL READ1s per leader (with several leaders)
  PREFIX.softclip_5p.tsv     most frequent READ1 5' soft-clipped sequences
  PREFIX.ctss_stats_mqc.tsv  one-row summary table for MultiQC
"""

import argparse
import collections
import sys

import pysam

from cagescan_sl import (MAX_SL_CLIP, first_g_mismatch, five_prime, parse_leaders, shift_downstream,
                         sl_status, sl_trimmed_5p, sl_trimmed_unplaceable, trimmed_leader)

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
    parser.add_argument("--sl-sequence", default="", help="spliced leader(s), 5'->3': SEQ or NAME=SEQ,NAME=SEQ,... (default: none)")
    parser.add_argument("--sl-min-overlap", type=int, default=8, help="minimum soft-clip length to call an SL match (default: %(default)s)")
    parser.add_argument("--remove-first-g", action="store_true",
                        help="READ1s aligned from their first base: a first G that differs from the genome is the non-templated G; the CTSS is the next base")
    parser.add_argument("--max-clip-report", type=int, default=30, help="truncate reported soft clips to this many 3'-most bases (default: %(default)s)")
    parser.add_argument("--top", type=int, default=50, help="number of soft-clip sequences to report (default: %(default)s)")
    args = parser.parse_args()
    leaders = parse_leaders(args.sl_sequence)
    leader_reads = collections.Counter()

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
            k, leader = sl_status(read, clip, aligned, leaders, args.sl_min_overlap)
            if k is not None:
                n_sl += 1
                leader_reads[leader] += 1
                if trimmed_leader(read) is not None:
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
    if leaders:
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

    if len(leaders) > 1:
        names = [name for name, _ in leaders] + sorted(set(leader_reads) - {name for name, _ in leaders})
        n = sum(leader_reads.values()) or 1
        with open(f"{args.prefix}.sl_leaders_mqc.tsv", "w") as out:
            out.write("# id: 'sl_leaders'\n# section_name: 'Spliced leaders'\n")
            out.write("# description: 'Share of spliced-leader READ1s carrying each leader.'\n# plot_type: 'bargraph'\n")
            out.write("# pconfig:\n#     id: 'sl_leaders_plot'\n#     title: 'Spliced leaders'\n#     ylab: '% of SL READ1s'\n#     cpswitch_c_active: False\n")
            out.write("Sample\t" + "\t".join(names) + "\n")
            out.write(args.prefix + "\t" + "\t".join(f"{100 * leader_reads[x] / n:.2f}" for x in names) + "\n")

    if n_used == 0:
        print(f"WARNING: no usable READ1 alignments found in {args.bam}", file=sys.stderr)


if __name__ == "__main__":
    main()
