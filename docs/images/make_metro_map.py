#!/usr/bin/env python3
"""Generate docs/images/cagescan_metro_map.svg (subway-style overview of the pipeline)."""

import os
from xml.sax.saxutils import escape

W, H = 1400, 760
BLUE, PURPLE, RED, ORANGE, GREEN, GOLD, GREY = "#1c7ed6", "#9c36b5", "#e03131", "#f08c00", "#2f9e44", "#d6336c", "#495057"
out = []


def line(d, color, width=10, dash=None):
    extra = f' stroke-dasharray="{dash}"' if dash else ""
    out.append(f'<path d="{d}" fill="none" stroke="{color}" stroke-width="{width}" stroke-linecap="round" stroke-linejoin="round"{extra}/>')


def station(x, y, name, sub=None, above=True, optional=False, tall=0, sub_below=False):
    dash = ' stroke-dasharray="4 3"' if optional else ""
    if tall:
        out.append(f'<rect x="{x - 11}" y="{y - 11}" width="22" height="{22 + tall}" rx="11" fill="#fff" stroke="#212529" stroke-width="3"{dash}/>')
    else:
        out.append(f'<circle cx="{x}" cy="{y}" r="10" fill="#fff" stroke="#212529" stroke-width="3"{dash}/>')
    ty = y - 22 if above else y + 32 + tall
    out.append(f'<text x="{x}" y="{ty}" class="st">{escape(name)}</text>')
    if sub and sub_below:
        for i, s in enumerate(sub if isinstance(sub, list) else [sub]):
            out.append(f'<text x="{x}" y="{y + 30 + tall + 13 * i}" class="sub">{escape(s)}</text>')
        return
    if sub:
        sy = ty - 15 if above else ty + 15
        for i, s in enumerate(sub if isinstance(sub, list) else [sub]):
            yy = sy - 13 * i if above else sy + 13 * i
            out.append(f'<text x="{x}" y="{yy}" class="sub">{escape(s)}</text>')


def terminus(x, y, w, label, color, sub=None):
    out.append(f'<rect x="{x - w / 2}" y="{y - 16}" width="{w}" height="32" rx="8" fill="{color}"/>')
    out.append(f'<text x="{x}" y="{y + 5}" class="term">{escape(label)}</text>')
    if sub:
        out.append(f'<text x="{x}" y="{y + 32}" class="sub">{escape(sub)}</text>')


# ---------------------------------------------------------------- lines
# row A: preprocessing and alignment (blue), coverage (purple)
line("M 70 120 L 660 120", BLUE)
line("M 660 120 L 1230 120", PURPLE)
# row B: TSS (red) splitting into broad (red) and sharp (orange)
line("M 660 120 C 705 120 705 270 750 270 L 1330 270", RED)
line("M 900 270 C 940 270 940 350 980 350 L 1330 350", ORANGE)
# turn down into row C (right to left), both sets in parallel
line("M 1330 270 C 1385 270 1385 520 1330 520 L 330 520", RED)
line("M 1330 350 C 1360 350 1360 536 1330 536 L 330 536", ORANGE)
# eRNA spur (find_eRNA branch)
line("M 830 545 L 830 660", GOLD, width=8, dash="10 8")

# ---------------------------------------------------------------- stations
terminus(70, 120, 120, "FASTQ pairs", GREY, "samplesheet (+ group)")
station(185, 120, "cat · FastQC", "merge lanes, read QC", sub_below=True)
station(295, 120, "cutadapt", "5' linker", optional=True, sub_below=True)
station(395, 120, "fastp", "adapters · quality", sub_below=True)
station(485, 120, "SortMeRNA", "rRNA", optional=True, sub_below=True)
station(565, 120, "STAR", "spliced alignment", sub_below=True)
station(660, 120, "samtools", "stats", sub_below=True)

station(790, 120, "samtools view", "primary, MAPQ", sub_below=True)
station(930, 120, "bedtools genomecov", "per strand, per million pairs", sub_below=True)
station(1080, 120, "bedGraphToBigWig")
terminus(1230, 120, 150, "IGV session", PURPLE, "genome · annotation · all tracks")

station(790, 270, "CTSS", ["READ1 5' ends; soft-clip report", "extra G · linker · spliced leader"], above=False)
station(880, 270, "pool", "all samples", above=True)
station(1020, 270, "paraclu", "pooled CTSS", above=True)
station(1160, 270, "trim · merge", "tails, < 20 bp gaps", above=True)
station(1010, 350, "paraclu", "per replicate", above=False)
station(1110, 350, "TPM / base", "(log) stability", above=False)
station(1210, 350, "pair · IDR", "within a group", above=False)
station(1300, 350, "innermost", "≤ 200 bp", above=False)

for i, (x, name, sub) in enumerate([
    (1250, "assign pairs", ["READ1 in a TSS cluster", "counts · CAGEscan clusters"]),
    (1110, "StringTie", "per sample, stranded"),
    (970, "merge", "all samples"),
    (830, "anchor", ["5' end → cluster peak", "derive · collapse 3' ends"]),
    (690, "gffread", "GTF · GFF3 · FASTA · BED12"),
    (550, "StringTie -e", "quantification"),
    (420, "gffcompare", "vs --gtf"),
]):
    station(x, 520, name, sub, above=(i % 2 == 1), tall=16)
terminus(290, 528, 120, "transcripts", GREEN, "per TSS, per set")
station(830, 660, "find eRNA", ["divergent non-SL pairs · distal", "balanced · unspliced · vs convergent"], above=False)
out.append('<text x="845" y="625" class="subl">find_eRNA branch</text>')

# ---------------------------------------------------------------- titles and legend
out.append('<text x="40" y="44" class="title">luscombeu/cagescan</text>')
out.append('<text x="40" y="66" class="subl">paired-end CAGE (CAGEscan) without UMIs → TSS clusters → 5\'-anchored transcripts</text>')
legend = [(BLUE, "reads and alignment"), (PURPLE, "coverage and genome browser"), (RED, "broad TSS clusters (pooled)"),
          (ORANGE, "sharp TSS clusters (RECLU, replicates)"), (GOLD, "eRNA candidates (find_eRNA)")]
lx, ly = 40, 640
for i, (c, t) in enumerate(legend):
    y = ly + 20 * i
    dash = ' stroke-dasharray="8 6"' if c == GOLD else ""
    out.append(f'<line x1="{lx}" y1="{y}" x2="{lx + 34}" y2="{y}" stroke="{c}" stroke-width="8" stroke-linecap="round"{dash}/>')
    out.append(f'<text x="{lx + 46}" y="{y + 4}" class="leg">{t}</text>')
out.append(f'<circle cx="{lx + 300}" cy="{ly}" r="7" fill="#fff" stroke="#212529" stroke-width="2.5" stroke-dasharray="4 3"/>')
out.append(f'<text x="{lx + 314}" y="{ly + 4}" class="leg">optional step</text>')
out.append(f'<text x="{lx + 300}" y="{ly + 24}" class="leg" text-anchor="start">Downstream of the TSS clusters every step</text>')
out.append(f'<text x="{lx + 300}" y="{ly + 40}" class="leg" text-anchor="start">runs once per cluster set (red and orange).</text>')
out.append(f'<text x="{lx + 300}" y="{ly + 60}" class="leg" text-anchor="start">QC of every step → MultiQC.</text>')

style = """<style>
text { font-family: "Helvetica Neue", Helvetica, Arial, sans-serif; fill: #212529; }
.title { font-size: 26px; font-weight: 700; }
.st { font-size: 14px; font-weight: 700; text-anchor: middle; }
.sub { font-size: 11px; fill: #495057; text-anchor: middle; }
.term { font-size: 13px; font-weight: 700; fill: #fff; text-anchor: middle; }
.leg { font-size: 12px; }
.subl { font-size: 12px; fill: #495057; text-anchor: start; }
</style>"""
svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">\n{style}\n'
       f'<rect width="{W}" height="{H}" fill="#ffffff"/>\n' + "\n".join(out) + "\n</svg>\n")
path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cagescan_metro_map.svg")
with open(path, "w") as f:
    f.write(svg)
print(path)
