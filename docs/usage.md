# luscombeu/cagescan: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

The pipeline takes **paired-end CAGE** reads (CAGEscan-style: READ1 starts at the capped 5' end of the RNA,
READ2 is random-primed downstream) **without UMIs**, and reconstructs 5'-anchored transcripts:

1. read pairs are aligned to the genome with STAR (spliced, local alignment);
2. READ1 5' ends give CAGE TSSs (CTSS), which are clustered with CAGEfightR (TPM, support across samples) or paraclu;
3. read pairs whose READ1 starts in a TSS cluster are kept and tagged with the cluster ID (`TC:Z:`);
4. StringTie builds isoforms from those pairs, per sample and merged across samples;
5. every transcript is matched to a TSS cluster, its 5' end is moved to the cluster's dominant CTSS and it is
   renamed `<tss_id>.<n>`; one TSS can therefore have several isoforms.

### Why not assemble per TSS like the original CAGEscan pipeline?

The original pipeline assembles the reads of each *molecule* (defined by UMI) de novo and maps the contigs.
Without UMIs a group of reads can only be defined by its TSS cluster, which mixes many molecules and isoforms.
The original assembler gives up on inconsistent read layouts and the transcript builder rejects alternative
splicing, so per-TSS de novo assembly fails exactly where isoforms exist. Grouping the *aligned* pairs by TSS and
assembling them with a splice-aware, isoform-aware assembler avoids that.

### Broad and sharp TSS clusters

TSS clusters decide which read pairs build which transcripts, and there is no single right granularity. The pipeline
therefore builds two cluster sets and runs everything downstream for both (`transcripts/broad/`,
`transcripts/sharp/`, ...), so you can choose afterwards:

| | broad (`TC...`) | sharp (`TS...`) |
| --- | --- | --- |
| Method | [CAGEfightR](https://bioconductor.org/packages/CAGEfightR) (default): per-sample TPM, CTSSs seen in >= `--ctss_min_samples` samples, pooled TPM clustered by slice-reduce (positions within `--tss_merge_distance` bases merged), clusters with >= `--tc_min_tpm` TPM in >= `--tc_min_samples` samples; or paraclu on the pooled counts (`--cluster_method paraclu`) | [RECLU](https://doi.org/10.1186/1471-2164-15-269) (Ohmiya et al. 2014): full paraclu hierarchy per replicate, TPM-per-base filter, hierarchical stability, replicate pairing (>= 90% reciprocal overlap) and IDR, <= 200 bp, innermost reproducible clusters |
| Recruits | every TSS region expressed above the TPM threshold, as one cluster | only TSSs reproducible across replicates, split into their narrow cores (typically a few bp) |
| Needs replicates | no | yes, >= 2 per group |
| Use for | maximum sensitivity; one model per promoter region | high-confidence, fine-grained TSSs; alternative TSSs a few bases apart kept separate |

Replicates are defined by the optional `group` column of the samplesheet (default: all samples are replicates of one
group). RECLU runs every pair of replicates within a group and unites the reproducible clusters of all pairs and
groups, so a TSS that is reproducible in only one condition is kept. Groups with a single sample contribute nothing to
the sharp set; with no group of two or more samples, the sharp set is skipped with a warning. `--skip_sharp_clusters`
turns it off.

**READ1 5' ends.** STAR aligns READ1 to its first base (`--alignEndsType Extend5pOfRead1`; mates and READ1 3'
ends are aligned locally). Local alignment would soft-clip first bases that differ from the genome (sequencing errors,
strain variants) and move the TSS downstream. The G that the reverse transcriptase adds at the cap is then removed
when it differs from the genome: the TSS is the next base (as `removeFirstG` in CAGEr). `--read1_end_to_end false`
restores local alignment; it is also used automatically with `--sl_sequence` when `--sl_split` is off, because the
leader must then be soft-clipped. Extra STAR arguments go in `--star_extra_args`.

**Broad clusters with CAGEfightR.** The default caller follows the standard workflow of CAGEfightR
([Thodberg et al. 2019](https://doi.org/10.1186/s12859-019-3029-5)), one of the two Bioconductor packages commonly
used for CAGE, with CAGEr: counts are normalised to tags per million (TPM) per sample, single positions (CTSSs) seen
in fewer than `--ctss_min_samples` (2) samples are dropped, the mean TPM across samples is clustered by
slice-reduce (positions with signal closer than `--tss_merge_distance`, 20 bp, are merged), the peak of a cluster
is its CTSS with the highest pooled TPM, and clusters below `--tc_min_tpm` (1) TPM in at least `--tc_min_samples` (2)
samples are dropped (default: 5 TPM in at least 1 sample). Thresholds are capped at the number of samples.

How strict? Common practice ranges from 1 TPM (CAGEr's tag clusters, `clusterCTSS(threshold = 1)`; the CAGEfightR
vignette, >= 1 TPM in a minimum number of samples; FANTOM5 "robust" peaks, >= 1 TPM and >= 11 reads at one position
in at least one sample) to 5 TPM (CAGEr's consensus clusters, `aggregateTagClusters(tpmThreshold = 5)`;
nf-core/cageseq). On the full Osaka MGI data (6 libraries, ~470 M READ1s), raising the threshold from 1 TPM in 2
samples to 5 TPM in 1 sample reduced non-SL clusters from 30,744 to 17,908 and raised the share whose peak has a
pyrimidine-purine initiator from 78% to 82%; 5 TPM in 1 sample keeps TSSs active in a single stage or condition.
Lower `--tc_min_tpm` for more sensitivity (e.g. weak enhancer RNAs) and raise it, or `--tc_min_samples`, for fewer,
stronger clusters. With `--sl_sequence`, SL and non-SL READ1s are
clustered separately (see below). `--cluster_method paraclu` restores the earlier caller (paraclu on raw pooled
counts, `--paraclu_min_cluster` reads), which recruits many more weak clusters at high depth.

**Deviation from RECLU.** RECLU ranks clusters for IDR by their hierarchical stability: a cluster's paraclu
stability (max density / min density) plus that of every cluster containing it. A stability depends on the distance
to the nearest stray CTSS, so it varies by large factors between replicates, and the strongest TSSs can come out
irreproducible. By default (`--reclu_score log-stability`) the pipeline sums the log stabilities instead. On simulated
replicates this recovered 109 of 157 TSSs, against 24 with RECLU's plain sum (`stability`) and 78 with expression
(`tpm`); none of the three produced false clusters. The IDR model is fitted with the
[idr](https://github.com/nboley/idr) package's own functions; its command-line peak merging is not used, because it
would collapse RECLU's nested clusters. With fewer than `--reclu_min_pairs` paired clusters, no model is fitted and
pairs are kept on overlap alone (reported in the MultiQC table).

`--reclu_min_tpm_per_base` (0.1, as in RECLU) is relative to sequencing depth: in very small libraries one read
already exceeds it, so the filter has no effect there.

### Enhancer-RNA candidates

Enhancer RNAs are transcribed in both directions from enhancers, so CAGE shows them as a **divergent pair** of TSS
clusters on opposite strands, facing away from each other, typically 100-400 bp apart, away from gene promoters,
with similar expression on both strands and short unspliced RNAs ([Andersson et al. 2014](https://doi.org/10.1038/nature12787)).
For each cluster set, `erna/` lists such pairs:

1. both clusters have >= `--erna_min_reads` READ1s and a spliced-leader fraction <= `--erna_max_sl_fraction`
   (SL-dominated clusters are trans-splice acceptors, not TSSs);
2. the minus-strand peak lies `--erna_min_distance` to `--erna_max_distance` bases upstream of the plus-strand peak;
3. with `--gtf`, both peaks are more than `--erna_distal_distance` bases from any annotated transcript 5' end;
4. directionality |F - R| / (F + R) < `--erna_max_directionality`;
5. neither cluster starts a spliced transcript.

Convergent pairs (plus-strand peak upstream of a minus-strand peak, same distances) go through the same filters as a
**control**: they arise by chance or from overlapping transcription, not from bidirectional initiation. The MultiQC
table reports both after each filter; only an excess of divergent pairs supports an eRNA-like signal, and on *O.
dioica* embryo CAGE the excess was about 3-fold (156 divergent vs 48 convergent candidates, broad set). Candidates
are hypotheses: unannotated genes and other short RNAs can look the same, so check them in the genome browser
(the session has an eRNA track) and against chromatin data. Use the broad set: the reproducibility filter of the
sharp set removes most weak, eRNA-like signal. `--find_erna false` turns this step off.

### Why StringTie, and what the pipeline corrects

StringTie is only trusted for what TSS-anchored CAGEscan pairs can support: the splice structure near the TSS.

- **Its input is pre-filtered.** StringTie only sees pairs whose READ1 5' end lies in a TSS cluster, so every fragment
  it assembles starts at a capped 5' end.
- **It knows the strand.** CAGE READ1 is sense to the RNA (`--fr`), and STAR adds the `XS` intron-motif tag, so sense and
  antisense transcription are kept apart.
- **It uses both mates of a pair.** READ1 and READ2 are one fragment, which links exons across the unsequenced insert;
  isoforms are paths through the splice graph built from junction-spanning reads.
- **It does not decide 5' ends.** StringTie starts a transcript at the most upstream read start; the pipeline moves
  every 5' end to the dominant CTSS of its TSS cluster (`tss_shift` records by how much).

On simulated data (`tests/data/make_test_data.py`, `tests/data/make_edge_case_data.py`) StringTie recovered exon
skipping near the TSS, alternative first exons, and all four combinations of two independent alternative exons
(including two at 10% each), with 100% intron-chain sensitivity and precision. It also showed four problems; the
first three are corrected by the pipeline:

| StringTie behaviour on CAGEscan data | Correction |
| --- | --- |
| A TSS inside continuously covered sequence (e.g. a second TSS in the same first exon) gets no transcript. | A transcript is **derived** for every TSS cluster lying inside an exon of an anchored transcript, starting at the cluster peak (`anchored "derived"`, `derived_from`; `--derive_transcripts false` to disable). |
| 3' ends sit where coverage runs out (coverage decays with fragment length), and one TSS can get several 3'-truncated copies of the same structure. | Transcripts that are 3' truncations of another transcript of the same TSS are dropped. 3' ends remain approximate. |
| paraclu can split one TSS (e.g. when the extra G aligns to a genomic G for part of the reads), and the transcript may anchor to the minor piece. | Clusters closer than `--tss_merge_distance` (20 bp) are merged, and transcripts anchor to the strongest cluster within `--tss_window`. |
| Isoform abundances assume uniform coverage, which CAGEscan does not have (one 4-isoform test gene: estimated 0.55/0.08/0.20/0.17 vs true 0.40/0.10/0.10/0.40). | Not corrected: use `tss_clusters/<set>/*.tss_cluster_counts.tsv` for expression; transcript TPMs are approximate. |

Splicing further downstream than the library insert size is never observed; that is a limit of CAGEscan data.

Derived transcripts assume that the downstream TSS shares the structure of the enclosing transcript. TSS clusters
inside exons can also come from CAGE artefacts (e.g. "exon painting"), so treat `anchored "derived"` transcripts as
lower confidence.

### READ1 structure, linkers and spliced leaders

Reverse transcriptase usually adds a non-templated G at the cap, so many READ1s start with one extra G. STAR's local
alignment soft-clips it, and the CTSS is taken from the first *aligned* READ1 base. When the genome itself has a G
at that position, the extra G aligns and the CTSS shifts 1 base upstream. This is a known limitation of CAGE: with a
genomic G upstream of the TSS, the dominant CTSS (and so the transcript 5' end) can be 1 base upstream of the true TSS.
Neighbouring clusters created by this shift are merged (`--tss_merge_distance`).

If you are unsure what precedes the cap in your reads, run the pipeline once with defaults and check
`ctss/<sample>.softclip_5p.tsv` and the *CTSS extraction* table in the MultiQC report:

| Most frequent soft clip | Meaning | What to do |
| --- | --- | --- |
| `-` (none) or `G` | READ1 starts at the cap (± extra G) | nothing |
| a constant longer sequence, e.g. `...TATAGGG` | a 5' linker is still present | `--r1_5p_linker TATAGGG` (cutadapt removes it and everything before it) |
| a constant ~20-40 nt sequence present only at some TSSs | spliced leader (SL trans-splicing) | `--sl_sequence <SL>`; SL READ1s are counted per sample and per TSS cluster |
| fixed-length random bases | barcode/UMI-like prefix | `--r1_trim_front <N>` |

Spliced-leader handling is only for organisms with SL trans-splicing (e.g. tunicates, nematodes, flatworms,
trypanosomes) and is off unless `--sl_sequence` is given; without it, no reads are split and no TSS is moved.
With `--sl_sequence`, the TSS of a trans-spliced READ1 is the trans-splice site, and SL reads are recognised in two
ways.

`--sl_sequence` takes one leader or several, for organisms with more than one (e.g. *C. elegans* SL1 and the SL2
family, used mostly by downstream genes of operons): `--sl_sequence SL1=GGTTTAATTACCCAAGTTTGAG,SL2a=GGTTTTAACCCAGTTACTCAAG,...`
(unnamed leaders are called SL, or SL1, SL2, ...). All leaders are looked for in every read; SL clusters pool them,
and `tss_clusters/<set>/all_samples.<set>.tss_cluster_sl_leaders.tsv` gives the SL READ1s per leader and the main
leader of each cluster; MultiQC shows the share of each leader per sample.

1. **Before alignment** (`--sl_split`, on by default). As in
   [splitSpliceLeaderPE](https://github.com/oist/plessy_splitspliceleaderpe), which splits raw read pairs by the
   presence of the SL at the start of READ1 (with TagDust2), the pairs are split into an SL and a non-SL set
   (cutadapt): the SL must start at READ1 base 1 or right after one G (the extra G), mismatches are allowed
   (`--sl_max_error_rate`, default 0.15) but no indels, so that the SL is cut exactly at its end, and the trimmed READ1
   starts at the trans-splice site. SL pairs are renamed `SL.<leader>:<name>` (or `SLG.<leader>:<name>` when the extra G
   was present), the others `no_adapter:<name>`. Both sets are then aligned together, and the read names carry the
   SL status and the leader into the BAM files (`samtools view -e 'qname=~"^SLG?\\."'` selects trans-spliced pairs). A trimmed READ1 holds no
   non-templated bases, so up to 3 bases that STAR soft-clips at its 5' end (sequencing errors or strain variants
   that differ from the assembly) are added back: its CTSS is the trans-splice site even when the first exon bases
   do not match the genome. Trimmed READ1s with a longer 5' clip (extra non-genomic bases after the leader, or an
   indel near the 5' end) cannot be placed and give no TSS (1-3% of READ1s in Osaka CAGE; MultiQC *CTSS extraction*).
   The split FASTQ files and cutadapt logs are in `sl_split/` (FASTQs only with `--save_sl_split_reads`, default on).
2. **After alignment**, in the READ1 5' soft clip, for SL reads not found before alignment (e.g. with only part of the
   leader, or with `--sl_split false`); exact match only.

Why both: with the soft-clip match alone, 11-15% of SL READ1s were missed on *O. dioica* MGI and Illumina CAGE,
because of sequencing errors in the 40-nt leader, and a further ~5% because STAR aligned the leader's last 4-27 bases to
the genome with mismatches. On 500,000 raw pairs per library, the cutadapt step called the same SL reads as
splitSpliceLeaderPE/TagDust2 (99.6% and 99.98% of TagDust2's SL reads; identical trimmed sequences for >99.9%),
in 4 s instead of 290 s, and also recovered clean SL reads that TagDust2 filtered as low complexity. Without indels,
cutadapt cut 100% of the SL reads exactly at the end of the leader (TagDust2 and cutadapt with indels: 99.3-99.5%).
Adding back the 5' soft clip matters: in the Osaka strain, the first base after the acceptor of a highly expressed
trans-spliced gene differs from the O10 assembly in most reads, and STAR's clip of that base put the SL peak one
base downstream of the acceptor. With these rules, 98.5% of the SL cluster peaks of the full Osaka data follow
an AG acceptor (CAGEr consensus clusters of oist/LuscombeU-CAGE_libraries: 97.2%).

**Separate clusters for trans-splice sites** (`--sl_separate_clusters`, on by default with `--sl_sequence`).
Trans-splice sites and TSSs of the same gene are often a few bases apart. Pooled into one cluster, only the stronger
of the two is reported as the peak. The broad set therefore clusters SL READ1s and the other READ1s separately:
`SL…` clusters (trans-splice sites) and `TC…` clusters (TSSs), which may overlap. SL READ1s are assigned to SL
clusters only and the other READ1s to TC clusters only, so the counts of a TC cluster are its non-SL READ1s and its
`tss_cluster_sl_counts` are 0. Compared with the CAGEr consensus clusters of the same Osaka libraries
(oist/LuscombeU-CAGE_libraries, which cluster SL and non-SL libraries separately), pooled clustering put the peak
on the trans-splice site for 830 TSS clusters (3.5% with a pyrimidine-purine initiator at the peak, against 83% for
the separate CAGEr clusters). The sharp (RECLU) set is not split.

The last bases of a spliced leader
often also match the genome, because the splice acceptor (`...AG`, `...CAG`) resembles the leader's 3' end (e.g.
`...AACAG` in *Oikopleura dioica*). STAR then aligns those leader bases, and the read's 5' end lands up to 3 bases
upstream, inside the acceptor. Such reads are recognised (soft clip ending in the leader minus its last k bases, and
the first k aligned bases equal to the leader's last k bases) and their CTSS is moved k bases downstream. On
*O. dioica* embryo CAGE this put 97% of trans-spliced TSS peaks exactly on annotated 5' ends, against 0% (62% at -2,
35% at -3) without the correction.

READ2 of CAGEscan libraries starts in the random primer, whose mismatches to the RNA are tolerated by the reverse
transcriptase. FANTOM5 CAGEscan ([Bertin et al. 2017](https://doi.org/10.1038/sdata.2017.147)) trimmed the first 6
bases of READ2 for random hexamers; use `--r2_trim_front 6` (or the length of your random primer) to do the same.

### Subsampling

By default each sample is randomly subsampled to at most 10 million read pairs (`--subsample_reads 10000000`, seqtk
`sample` with the same seed for READ1 and READ2, so pairs stay together; samples with fewer pairs are used whole).
This bounds run time and disk use for very large libraries (an MGI lane can hold several hundred million pairs per
sample; the BAMs need about 0.4 GB per million pairs). Subsampling is the first step after merging lanes, so FastQC,
the alignment statistics and all counts refer to the subsampled reads; CPM/TPM values stay comparable between samples.

10 M pairs are enough to place the TSSs of moderately and highly expressed genes and to build their transcripts. Use
`--subsample_reads 0` (all reads) when depth matters: weak TSSs, enhancer RNAs, lowly expressed genes, or quantifying
rare isoforms. Change `--subsample_seed` to draw a different random subset.

### Relation to the FANTOM5 CAGEscan workflow

The FANTOM5 CAGEscan libraries ([Bertin et al. 2017](https://doi.org/10.1038/sdata.2017.147)) had no UMIs either.
That workflow aligned the pairs (BWA), kept proper pairs, removed pairs with identical coordinates, and grouped pairs
whose READ1 starts in a CAGE peak into *CAGEscan clusters* (union of the pairs' blocks). This pipeline follows the same
logic (STAR instead of BWA for spliced alignment, paraclu TSS clusters as seeds) and writes the same kind of clusters
to `cagescan_clusters/`, but additionally resolves isoforms with StringTie.

### Duplicates

Without UMIs, PCR duplicates cannot be told apart from independent molecules with the same TSS. With
`--dedup`, read pairs with identical READ1 *and* READ2 positions are marked by `samtools markdup` and ignored
downstream. Because READ2 is random-primed this is a reasonable proxy (FANTOM5 CAGEscan removed such pairs), but
highly expressed TSSs with short fragments also collide by chance, so it is off by default.

### Reference annotation

`--gtf` is optional. When given, it is used for the STAR splice-junction database and to compare the merged
transcripts with the annotation (gffcompare). By default it is **not** given to StringTie as a guide, so that
transcript 5' ends come from the CAGE data; use `--stringtie_guide` to change this.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CAGE_REP1,CAGE_REP1_L002_R1.fastq.gz,CAGE_REP1_L002_R2.fastq.gz
CAGE_REP1,CAGE_REP1_L003_R1.fastq.gz,CAGE_REP1_L003_R2.fastq.gz
CAGE_REP1,CAGE_REP1_L004_R1.fastq.gz,CAGE_REP1_L004_R2.fastq.gz
```

### Full samplesheet

Reads must be paired-end. `fastq_1` must be READ1, the read that starts at the capped 5' end. The optional `group`
column names the condition a sample is a replicate of; it is only used for the sharp (RECLU) TSS clusters.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2,group
LIVER_REP1,LIVER_REP1_R1.fastq.gz,LIVER_REP1_R2.fastq.gz,liver
LIVER_REP2,LIVER_REP2_R1.fastq.gz,LIVER_REP2_R2.fastq.gz,liver
BRAIN_REP1,BRAIN_REP1_R1.fastq.gz,BRAIN_REP1_R2.fastq.gz,brain
BRAIN_REP2,BRAIN_REP2_R1.fastq.gz,BRAIN_REP2_R2.fastq.gz,brain
```

| Column    | Description                                                                                                                                  |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Sample name: letters, digits, `_` and `-` only. Identical for multiple sequencing runs of the same library.                                  |
| `fastq_1` | Full path to the gzipped READ1 FastQ file (cap side). Extension ".fastq.gz" or ".fq.gz".                                                     |
| `fastq_2` | Full path to the gzipped READ2 FastQ file. Extension ".fastq.gz" or ".fq.gz".                                                                |
| `group`   | Optional. Condition the sample is a replicate of (letters, digits, `_`, `-`). Default: all samples form one group.                           |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run luscombeu/cagescan --input ./samplesheet.csv --fasta genome.fa --gtf annotation.gtf --outdir ./results -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run luscombeu/cagescan -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull luscombeu/cagescan
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [luscombeu/cagescan releases page](https://github.com/luscombeu/cagescan/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
