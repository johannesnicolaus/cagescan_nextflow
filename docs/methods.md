# luscombeu/cagescan: Methods

This document describes what the pipeline does, step by step, with the exact tool arguments and default values, so
that it can be cited and reproduced (e.g. in a methods section). Parameter names are given as `--name` (default).
User-facing documentation is in [usage](usage.md) and [output](output.md).

## Overview

luscombeu/cagescan processes paired-end CAGE libraries without unique molecular identifiers (UMIs): CAGEscan libraries
(READ2 random-primed downstream of the cap) and other paired-end CAGE libraries (READ1 at the cap).

- **TSS calling**: CAGE transcription start sites (CTSSs) are the 5' ends of READ1,
  clustered into TSS clusters per experiment, with trans-splice sites (reads carrying a spliced leader, SL) called
  separately from transcription start sites.
- **Transcripts** ("CAGEscan"; Plessy et al. 2010; Bertin et al. 2017): read pairs whose READ1 starts
  in a TSS cluster are assembled into isoforms, whose 5' ends are anchored to the cluster peak.

The pipeline is written in Nextflow
(>= 26.04) with the nf-core template (4.1.0) and runs every tool in a container.

## 1. Reads

1. **Lanes** of one sample are concatenated.
2. **Subsampling** (`--subsample_reads` 10,000,000 reads or read pairs per sample; 0 = all): `seqtk sample -2 -s
   <--subsample_seed 100>`, the same seed for READ1 and READ2 so that pairs stay together; samples with fewer reads
   are used whole.
3. **Quality control**: FastQC.
4. **Linker** (optional, `--r1_5p_linker`): `cutadapt -g <linker> [--discard-untrimmed] --minimum-length 20`
   removes the linker and everything before it from READ1.
5. **Adapters and quality** (unless `--skip_fastp`): `fastp --detect_adapter_for_pe --length_required 20
   [--trim_front1 <--r1_trim_front>] [--trim_front2 <--r2_trim_front>]`.
6. **rRNA** (optional, `--ribo_database_manifest`): SortMeRNA.
7. **Spliced leaders** (with `--sl_sequence`, unless `--sl_split false`). One or several leaders
   (`NAME=SEQ,NAME=SEQ`, or `SEQ`) are searched at the start of READ1 with cutadapt, as in splitSpliceLeaderPE
   (github.com/oist/plessy_splitspliceleaderpe; there with TagDust2): for every leader two anchored 5' patterns,
   `^SEQ` and `^GSEQ` (the non-templated G before the leader), `--no-indels`, `-e <--sl_max_error_rate 0.15>`
   mismatches, `--minimum-length 15`. Indels are not allowed so that the leader is cut exactly at its last base. Reads
   with a leader are trimmed and renamed `SL.<leader>:<name>` / `SLG.<leader>:<name>`, the others
   `no_adapter:<name>`; both sets are aligned together, the name carrying the SL status and the leader. On 500,000
   raw read pairs per library (O. dioica Osaka MGI and Barcelona Illumina CAGE), this step called 99.6% and 99.98% of
   the SL reads called by TagDust2, cut 100% of SL reads at the end of the leader (TagDust2: 99.2-99.5%), took 4 s
   instead of 290 s, and also kept SL reads that TagDust2 placed in its non-SL output.

## 2. Alignment

- **Reference**: `--fasta`; with `--gtf`, records with an empty `transcript_id` (gene rows of NCBI-style GTFs) are
  dropped first; all other records are used unchanged (STAR splice junctions, StringTie guide, gffcompare).
- **STAR** 2.7.11b (genome index built with the GTF when given): `--outSAMtype BAM SortedByCoordinate
  --outSAMstrandField intronMotif --outSAMattributes NH HI AS nM NM MD --outFilterMultimapNmax 20 --outFilterType
  BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignEndsType Extend5pOfRead1` (`--read1_end_to_end`;
  `Local` with `--read1_end_to_end false`, and always when a spliced leader is given but not split before
  alignment, so that it can be soft-clipped), plus `--star_extra_args`. READ1 is aligned from its first base;
  mates and READ1 3' ends are aligned locally.
- **Duplicates** (optional, `--dedup`): `samtools collate`, `fixmate`, `sort`, `markdup`; pairs with
  identical READ1 and READ2 positions are marked and ignored downstream.
- **Alignment QC**: `samtools stats`, `flagstat`, `idxstats`.

## 3. CTSS

For each sample, every READ1 that is primary, mapped, not QC-failed, not marked as duplicate and
has MAPQ >= `--min_mapq` (10) gives one CTSS at its 5' end: the first aligned base on the plus strand
(`reference_start`, 0-based), the last aligned base on the minus strand (`reference_end - 1`). Then:

- **Non-templated G** (READ1 aligned from its first base): if the first base is a G that differs from the genome
  (from the MD tag; on reverse-strand reads the last base, C in the BAM), the CTSS is the next base downstream
  (as `removeFirstG` in CAGEr; Haberle et al. 2015). A G matching a genomic G cannot be recognised and is kept.
  CAGEr's statistical correction of such cases (`correctSystematicG`) is not applied: on the Osaka data it lowered
  the share of TSS cluster peaks with a pyrimidine-purine initiator (90.6% to 88.8%).
- **SL-trimmed reads** start at the trans-splice site. If the aligner soft-clipped up to 3 bases at their 5' end
  (first exon bases that differ from the assembly, e.g. strain variants), those bases are added back; reads with
  longer 5' clips (non-genomic sequence after the leader, or an indel near the 5' end) give no CTSS.
- **Local alignment only**: the 5' soft clip is tallied (`softclip_5p.tsv`); a clip ending exactly with a leader
  marks an SL read, and when the leader's last k <= 3 bases also aligned to the genome (the acceptor `...AG`
  resembles the leader's 3' end) the CTSS is moved k bases downstream, to the trans-splice site.

CTSS are written per sample, and with `--sl_sequence` separately for SL and non-SL reads.

## 4. TSS clusters

### Broad clusters (default: CAGEfightR)

With `--cluster_method cagefightr` (default), CAGEfightR 1.26.0 (Thodberg et al. 2019) is run on the per-sample
CTSS, separately for non-SL and SL CTSS when `--sl_sequence` is given (`--sl_separate_clusters`):

1. counts are converted to tags per million (TPM) per sample and class: non-SL CTSS relative to the sample's non-SL
   READ1s, SL CTSS relative to its SL READ1s (so that a sample's share of trans-spliced reads does not change the
   threshold; normalising both to all READ1s removed 38% of the TSS clusters, of which 86.9% had a pyrimidine-purine
   initiator);
2. CTSSs with reads in fewer than `--ctss_min_samples` (2) samples are dropped (`subsetBySupport`);
3. the pooled signal (TPM summed across samples, `calcPooled`) is clustered by slice-reduce
   (`clusterUnidirectionally`, `pooledCutoff = 0`, `mergeDist = --tss_merge_distance` (20): positions with signal
   fewer than 20 bases apart are merged); the peak is the CTSS with the highest pooled TPM;
4. clusters are quantified (`quantifyClusters`) and kept if they have >= `--tc_min_tpm` (5) TPM in >=
   `--tc_min_samples` (1) samples.

Sample-count thresholds are capped at the number of samples. The 5 TPM in 1 sample default follows CAGEr's consensus
clusters (`aggregateTagClusters(tpmThreshold = 5)`) and keeps TSSs active in a single sample; lenient alternatives in
common use are 1 TPM (CAGEr tag clusters, the CAGEfightR vignette, FANTOM5 robust peaks with >= 11 reads at one
position and >= 1 TPM in one sample). On the Osaka data, 1 TPM in 2 samples gave 30,744 non-SL clusters (77.6% with a
pyrimidine-purine initiator at the peak), 5 TPM in 1 sample 17,908 (82.0%) (before the alignment change of section 2).

Non-SL clusters are named `TC<n>`, SL clusters `SL<n>`; the two kinds may overlap, so that a TSS next to a
trans-splice site keeps its own peak.

### Broad clusters with paraclu

With `--cluster_method paraclu`, the CTSS pooled across samples (per class) are clustered with paraclu 10 (Frith et
al. 2008) with a minimum of `--paraclu_min_cluster` (30) reads and simplified with `paraclu-cut` (patched to sort by
strand); sparse tails are trimmed (up to 1% of a cluster's reads per end) and clusters fewer than
`--tss_merge_distance` bases apart are merged; the peak is the CTSS with the most reads. paraclu does not report
single-position clusters, which can make it miss very sharp peaks (e.g. trans-splice sites).

### Sharp clusters (RECLU)

Unless `--skip_sharp_clusters`, reproducible clusters are built from replicates (samplesheet `group`) following RECLU
(Ohmiya et al. 2014): per replicate, the full paraclu hierarchy (`--reclu_paraclu_min` 2), clusters with >=
`--reclu_min_tpm_per_base` (0.1) TPM per base, scored by hierarchical stability (`--reclu_score log-stability`: the sum
of log stabilities, which was more reproducible than RECLU's plain sum on simulated replicates); replicate clusters
paired by >= `--reclu_min_overlap` (0.9) reciprocal overlap; irreproducible discovery rate (IDR; Li et al. 2011,
`idr` 2.0.4.2) < `--reclu_max_idr` (0.1) when >= `--reclu_min_pairs` (50) pairs, length <= `--reclu_max_length` (200);
the innermost reproducible clusters of all replicate pairs are kept (`TS<n>`). Everything downstream runs for both
cluster sets.

## 5. Counts

A READ1 is assigned to the cluster (of its own kind: SL or non-SL) whose span extended by `--tss_window` (50) bases
contains its 5' end, on the same strand (the nearest one if several). Paired-end READ1s must be in a proper pair.
Per cluster and sample: READ1 counts, counts per million assigned READ1s, SL READ1 counts and, with several leaders,
SL READ1 counts per leader (with the main leader of each cluster). Assigned reads (pairs) are written to
TSS-anchored BAMs with a `TC:Z:<cluster>` tag.

## 6. Transcripts

1. **CAGEscan clusters** (Bertin et al. 2017): per TSS cluster, the union of the aligned blocks of its read pairs, per
   sample and merged across samples.
2. **Assembly**: StringTie 3.0.3 on the TSS-anchored pairs of each sample (`--fr -m <--min_transcript_length 100>`,
   guided by `--gtf` with `--stringtie_guide`), merged across samples (`stringtie --merge -m 100`).
3. **Anchoring**: each transcript's 5' end is matched to a TSS cluster on the same strand within `--tss_window` (the
   cluster with the most reads if several) and moved to its peak; transcripts are renamed `<tss_id>.<n>` (by
   coverage); unanchored transcripts are dropped (unless `--keep_unanchored`). A transcript that is a 3' truncation
   of another of the same TSS is dropped. For a TSS cluster inside an exon of an anchored transcript, or up to
   `--tss_window` bases upstream of its start, a transcript is derived from it (cut or extended to the peak;
   `--derive_transcripts`).
4. **Formats**: GTF, GFF3 (`gffread --keep-genes`), BED12 (5' block thick when anchored), FASTA (`gffread -w`);
   gffread 0.12.7.
5. **Quantification**: StringTie `-e --fr` per sample on the merged anchored transcripts (TPM, coverage). CAGEscan
   coverage decays with distance from the TSS, so transcript abundances are approximate; TSS cluster counts are the
   recommended expression measure.
6. **Annotation comparison** (with `--gtf`): gffcompare 0.12.6.

CAGEscan transcripts span from the TSS to the end of the sequenced fragment (a few hundred bases to ~1 kb): they are
5'-anchored transcript fragments, not full-length transcripts.

## 7. Enhancer-RNA-like divergent TSS pairs

Unless `--find_erna false`, per cluster set (Andersson et al. 2014): non-SL clusters with >= `--erna_min_reads` (10)
READ1s and an SL fraction <= `--erna_max_sl_fraction` (0.1); pairs of a minus-strand peak upstream of a plus-strand peak
`--erna_min_distance` (50) to `--erna_max_distance` (400) bases apart (divergent); both peaks > `--erna_distal_distance`
(500) bases from annotated transcript 5' ends (with `--gtf`); directionality |F - R| / (F + R) <
`--erna_max_directionality` (0.8); and no multi-exon transcript starting at either cluster. Convergent
pairs (plus-strand peak upstream of a minus-strand peak) go through the same filters as a chance control; the
divergent/convergent ratio after each filter is reported.

## 8. Signal tracks and genome browser

- **TSS signal**: CTSS per strand, per sample and pooled, as bigWig (bedGraphToBigWig, UCSC 482).
- **Coverage**: `samtools view -b -q <min_mapq> -F 0x904` (`0xD04`, i.e. also without duplicates, with `--dedup`),
  `bedtools genomecov -bg -split -du -strand <+/->` (bedtools 2.31.1; both mates on the RNA strand), scaled per million
  READ1s used for CTSS, as bigWig.
- **IGV session** (`igv_session.xml`) with copies of the genome and annotation, the cluster sets, transcripts, eRNA
  candidates and every produced signal track.

## 9. Quality metrics (MultiQC 1.35)

Besides the tool reports (FastQC, cutadapt, fastp, SortMeRNA, STAR, samtools, gffcompare), per sample:

- **CTSS extraction**: READ1s used, distinct CTSS, 5' soft clips, extra G removed, SL reads (trimmed before
  alignment, clips added back, unplaceable, shifted at the acceptor).
- **TSS initiator**: the dinucleotide at positions -1/+1 of non-SL CTSS, read-weighted, and the share with a
  pyrimidine-purine (YR, "PyPu") initiator, the canonical initiator of RNA polymerase II.
- **Trans-splice sites**: the two bases before SL CTSS, read-weighted, and the share after an AG acceptor.
- **Spliced leaders** (several leaders): the share of SL READ1s per leader.
- TSS cluster assignment, transcript anchoring, sharp-cluster (RECLU) and eRNA statistics.

## 10. Validation on *O. dioica* CAGE

The defaults were evaluated on paired-end CAGE of *Oikopleura dioica* (Osaka strain; MGI, six libraries from three
developmental stages, ~470 M READ1s; genome O10), against the consensus clusters of the same libraries in
[oist/LuscombeU-CAGE_libraries](https://github.com/oist/LuscombeU-CAGE_libraries) (TagDust2 SL splitting, HISAT2,
CAGEr). As the identity of true TSSs is unknown, two sequence features independent of both pipelines were used: the
initiator dinucleotide at TSS peaks (pyrimidine-purine expected) and the two bases before trans-splice-site peaks (AG
expected). Annotated 5' ends were not used to compare peak positions, because the O10 v1.3 5' UTRs were derived from
the CAGEr clusters themselves.

| | this pipeline (commit b25ad8a) | CAGEr clusters |
| --- | --- | --- |
| TSS (non-SL) clusters / SL clusters | 17,203 / 6,572 | 13,969 / 5,939 |
| non-SL peaks with a PyPu initiator | 90.5% | 88.8% |
| SL peaks right after AG | 98.7% | 97.2% |
| CAGEr clusters found (same strand, <= 20 bp; SL / non-SL) | 99.7% / 89.1% (100% / 98.9% of their signal) | |
| where the peaks differ, PyPu (ours vs theirs) | 83.1% vs 63.5% | |
| READ1s in clusters that the CAGEr set also has | 99.2% | |

Full-depth run of commit b25ad8a (2026-10-19; 466 M READ1s after leader splitting and alignment, 380 M in TSS
clusters; eRNA candidates: 76 divergent vs 11 convergent pairs). The comparison scripts are in
`runs/osa_MGI_220209/compare_CAGE_libraries/` of the analysis directory, not in the pipeline. Individual changes were evaluated in the same way:

| change | effect |
| --- | --- |
| SL split before alignment (vs exact soft-clip match) | +35% SL READ1s recognised; their clusters 98.8% matched to CAGEr SL clusters |
| separate SL and non-SL clustering | removed 830 non-SL clusters whose peak was a trans-splice site (3.5% PyPu) |
| added-back 5' clips of SL reads | moved trans-splice peaks off the first exon base (e.g. a strain variant after the acceptor) |
| CAGEfightR (vs paraclu) | paraclu missed single-position trans-splice peaks |
| READ1 aligned from its first base + first-G removal | non-SL peaks with PyPu 82.0% -> 90.5%; 14% of READ1s unaligned (non-genomic 5' ends) |
| systematic G correction | not adopted (PyPu 90.6% -> 88.8%) |

## Software

Containers pinned in this version (BioContainers / Seqera Wave):

| tool | version |
| --- | --- |
| Nextflow | >= 26.04 (template nf-core 4.1.0, nf-schema 2.7.3) |
| cat / gzip (coreutils) | container `coreutils_grep_gzip_lbzip2_pruned` |
| seqtk | 1.4 |
| FastQC | 0.12.1 |
| cutadapt | 5.2 |
| fastp | 1.3.6 |
| SortMeRNA | 4.3.7 |
| STAR | 2.7.11b (with samtools 1.21) |
| samtools | 1.24 |
| Python / pysam | 3.12 / 0.24.1 (pipeline scripts in `bin/`) |
| CAGEfightR | 1.26.0 (R / Bioconductor) |
| paraclu | 10 |
| idr | 2.0.4.2 |
| StringTie | 3.0.3 |
| gffread / gffcompare | 0.12.7 / 0.12.6 |
| bedtools | 2.31.1 |
| UCSC bedGraphToBigWig | 482 |
| MultiQC | 1.35 |

The versions used by a run are written to `pipeline_info/cagescan_software_mqc_versions.yml`.

## Reproducibility

- Run a fixed version: `nextflow run johannesnicolaus/cagescan_nextflow -r <release or commit> -profile
  <singularity|docker|...> -params-file params.yaml`. `pipeline_info/params_<timestamp>.json` records every parameter
  value of a run.
- Defaults change between versions; results of runs made with different versions are not directly comparable.
- Test data: `tests/data/make_test_data.py` simulates paired-end CAGE of 156 genes on a 1 Mb genome, with known TSSs,
  isoforms and a trans-spliced gene; the nf-test suite runs the pipeline (default, paraclu,
  NCBI-style annotation) and every local module on it.

## References

See [CITATIONS.md](../CITATIONS.md); the MultiQC report's methods section lists the tools and references of each run.
