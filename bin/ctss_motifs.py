#!/usr/bin/env python3
"""Sequence at CAGE TSSs and trans-splice sites, as MultiQC custom-content tables (read-weighted).

  ctss_motifs.py --fasta genome.fa --prefix P --ctss S.ctss.bed [--sl-ctss S.sl_ctss.bed]

Outputs (one row per sample, merged by MultiQC):
  P.tss_initiator_mqc.tsv  TSS initiator: dinucleotide at -1/+1 of the CTSS (read-weighted); pyrimidine-purine
                           (YR, e.g. CA, TG) is the canonical initiator.  With --sl-ctss, --ctss should hold the
                           non-SL CTSS, since trans-splice sites are not transcription initiation sites.
  P.sl_acceptor_mqc.tsv    trans-splice sites (with --sl-ctss): the two bases before the CTSS (AG = splice acceptor)
  P.ctss_motif_stats_mqc.tsv  summary: % YR at TSSs, % AG before trans-splice sites
"""

import argparse
import collections

import pysam

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
BASES = "ACGT"
DINUCS = [a + b for a in BASES for b in BASES]


def oriented(genome, chrom, start, end, strand):
    """Genome sequence [start, end) on the given strand (5'->3'); '' if out of range."""
    if start < 0 or end > genome.get_reference_length(chrom):
        return ""
    seq = genome.fetch(chrom, start, end).upper()
    return seq if strand == "+" else seq.translate(COMPLEMENT)[::-1]


def count(genome, path, before, after):
    """Read-weighted counts of the sequence from `before` bases upstream to `after` bases downstream of each CTSS."""
    counts = collections.Counter()
    with open(path) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or f[0] not in genome.references:
                continue
            pos, strand, n = int(f[1]), f[5], int(float(f[4]))
            start, end = (pos - before, pos + after) if strand == "+" else (pos - after + 1, pos + before + 1)
            counts[oriented(genome, f[0], start, end, strand)] += n
    return counts


def bargraph(path, pid, section, description, ylab, sample, counts, categories):
    total = sum(counts[c] for c in categories) or 1
    with open(path, "w") as out:
        out.write(f"# id: '{pid}'\n# section_name: '{section}'\n# description: '{description}'\n")
        out.write(f"# plot_type: 'bargraph'\n# pconfig:\n#     id: '{pid}_plot'\n#     title: '{section}'\n#     ylab: '{ylab}'\n#     cpswitch_c_active: False\n")
        out.write("Sample\t" + "\t".join(categories) + "\n")
        out.write(sample + "\t" + "\t".join(f"{100 * counts[c] / total:.2f}" for c in categories) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True, help="indexed genome FASTA")
    parser.add_argument("--prefix", required=True, help="output prefix (sample name)")
    parser.add_argument("--ctss", required=True, help="CTSS bed of the TSSs (non-SL READ1s when --sl-ctss is given)")
    parser.add_argument("--sl-ctss", help="CTSS bed of spliced-leader READ1s (trans-splice sites)")
    args = parser.parse_args()
    genome = pysam.FastaFile(args.fasta)

    init = count(genome, args.ctss, 1, 1)  # bases -1 and +1
    bargraph(f"{args.prefix}.tss_initiator_mqc.tsv", "tss_initiator", "TSS initiator (-1/+1)",
             "Dinucleotide at positions -1/+1 of CAGE TSSs (non-SL READ1s), weighted by reads. Pyrimidine-purine (YR: CA, CG, TA, TG) is the canonical initiator.",
             "% of READ1s", args.prefix, init, DINUCS)
    yr = sum(init[d] for d in DINUCS if d[0] in "CT" and d[1] in "AG")
    n_init = sum(init[d] for d in DINUCS) or 1
    header, row = ["Sample", "TSS with YR initiator (%)"], [args.prefix, f"{100 * yr / n_init:.2f}"]
    if args.sl_ctss:
        acc = count(genome, args.sl_ctss, 2, 0)  # the two bases before the trans-splice site
        bargraph(f"{args.prefix}.sl_acceptor_mqc.tsv", "sl_acceptor", "Trans-splice sites (-2/-1)",
                 "The two bases before CAGE TSSs of spliced-leader READ1s, weighted by reads. Trans-splice sites follow an AG splice acceptor.",
                 "% of SL READ1s", args.prefix, acc, DINUCS)
        n_acc = sum(acc[d] for d in DINUCS) or 1
        header.append("SL sites after AG (%)")
        row.append(f"{100 * acc['AG'] / n_acc:.2f}")
    with open(f"{args.prefix}.ctss_motif_stats_mqc.tsv", "w") as out:
        out.write("# id: 'ctss_motif_stats'\n# section_name: 'TSS and trans-splice site sequence'\n")
        out.write("# description: 'Share of TSS READ1s with a pyrimidine-purine initiator, and of SL READ1s right after an AG acceptor.'\n# plot_type: 'table'\n")
        out.write("\t".join(header) + "\n" + "\t".join(row) + "\n")


if __name__ == "__main__":
    main()
