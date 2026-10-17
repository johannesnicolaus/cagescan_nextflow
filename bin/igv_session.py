#!/usr/bin/env python3
"""Write a self-contained IGV session for the pipeline results.

The genome FASTA, its index and the annotation are copied into igv/, so the
session only uses paths relative to the results directory: the whole results
directory (or just the session, igv/ and the track files) can be copied and
opened anywhere.

usage: igv_session.py --fasta genome.fa --fai genome.fa.fai [--gtf annotation.gtf]
                      --samples S1,S2 --sets broad,sharp [--erna]
  -> igv_session.xml, igv/genome.fa, igv/genome.fa.fai, igv/annotation.gtf
"""

import argparse
import os
import shutil
from xml.sax.saxutils import quoteattr

FEATURE = "org.broad.igv.track.FeatureTrack"
DATA = "org.broad.igv.track.DataSourceTrack"


def track(path, name, clazz, **attrs):
    attrs = {"clazz": clazz, "id": path, "name": name, "visible": "true", **attrs}
    return "        <Track " + " ".join(f"{k}={quoteattr(str(v))}" for k, v in attrs.items()) + "/>"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--gtf", default=None)
    parser.add_argument("--samples", required=True, help="comma-separated sample names")
    parser.add_argument("--sets", required=True, help="comma-separated TSS cluster sets")
    parser.add_argument("--erna", action="store_true", help="add the eRNA candidate tracks")
    args = parser.parse_args()

    os.makedirs("igv", exist_ok=True)
    shutil.copyfile(args.fasta, "igv/genome.fa")
    shutil.copyfile(args.fai, "igv/genome.fa.fai")
    tracks = []
    if args.gtf:
        shutil.copyfile(args.gtf, "igv/annotation.gtf")
        tracks.append(track("igv/annotation.gtf", "Gene models (--gtf)", FEATURE, displayMode="EXPANDED"))
    group = 0
    for s in args.sets.split(","):
        tracks.append(track(f"transcripts/{s}/merged/merged.{s}.transcripts.bed12", f"Transcripts ({s} TSS clusters)", FEATURE, displayMode="EXPANDED"))
        tracks.append(track(f"tss_clusters/{s}/all_samples.{s}.tss_clusters.bed", f"TSS clusters ({s})", FEATURE, displayMode="COLLAPSED"))
        if args.erna:
            tracks.append(track(f"erna/{s}/all_samples.{s}.erna_candidates.bed", f"eRNA-like divergent pairs ({s})", FEATURE, displayMode="EXPANDED", color="120,0,160"))
    group += 1
    for strand, color in (("plus", "200,0,0"), ("minus", "0,0,200")):
        tracks.append(track(f"ctss/bigwig/all_samples.ctss.{strand}.bigWig", f"TSS signal {'+' if strand == 'plus' else '-'} (all samples)", DATA,
                            autoScale="true", autoscaleGroup=group, color=color, renderer="BAR_CHART"))
    for sample in args.samples.split(","):
        group += 1
        for strand, color in (("plus", "200,0,0"), ("minus", "0,0,200")):
            tracks.append(track(f"coverage/{sample}.coverage.{strand}.bigWig", f"{sample} coverage {'+' if strand == 'plus' else '-'}", DATA,
                                autoScale="true", autoscaleGroup=group, color=color, renderer="BAR_CHART"))

    with open("igv_session.xml", "w") as out:
        out.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
        out.write('<Session genome="igv/genome.fa" locus="All" relativePath="true" version="8">\n')
        out.write("    <Resources>\n")
        for t in tracks:
            path = t.split('id="', 1)[1].split('"', 1)[0]
            out.write(f'        <Resource path="{path}"/>\n')
        out.write("    </Resources>\n")
        out.write('    <Panel name="DataPanel">\n')
        out.write("\n".join(tracks) + "\n")
        out.write("    </Panel>\n")
        out.write("</Session>\n")


if __name__ == "__main__":
    main()
