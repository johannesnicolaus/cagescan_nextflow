# luscombeu/cagescan: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Transcripts](#transcripts) - the main result: TSS-anchored transcripts (GTF, GFF3, BED12, FASTA)
- [TSS clusters](#tss-clusters) - consensus TSS clusters, counts and TSS-anchored BAMs
- [CAGEscan clusters](#cagescan-clusters) - FANTOM5-style union of read pairs per TSS cluster
- [CTSS](#ctss) - per-sample CAGE TSSs, bigWigs and READ1 5' soft-clip reports
- [Quantification](#quantification) - per-sample transcript abundance
- [gffcompare](#gffcompare) - comparison with the reference annotation
- [Alignment](#alignment) - STAR alignments and statistics
- [Read QC and trimming](#read-qc-and-trimming) - FastQC, fastp, cutadapt, SortMeRNA
- [StringTie assemblies](#stringtie-assemblies) - raw assemblies before TSS anchoring
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Transcripts

<details markdown="1">
<summary>Output files</summary>

- `transcripts/merged/` - consensus transcript set across all samples (use this one for most analyses)
  - `merged.transcripts.gtf`: TSS-anchored transcripts.
  - `merged.transcripts.gff3`: the same in GFF3 (gffread).
  - `merged.transcripts.bed12`: the same in BED12, for genome browsers.
  - `merged.transcripts.fasta`: spliced transcript sequences (gffread `-w`).
- `transcripts/<sample>/` - the same four files for each sample's own assembly.

</details>

Transcripts are grouped by TSS cluster: `gene_id` is the TSS cluster ID (e.g. `TC000012`) and `transcript_id` is
`<tss_id>.<n>`, numbered by decreasing StringTie coverage, so a TSS with several isoforms has several transcripts.
GTF attributes:

| Attribute | Meaning |
| --- | --- |
| `tss_id` | TSS cluster the transcript starts in |
| `tss_peak` | dominant CTSS of the cluster (`chrom:position:strand`, 1-based); the transcript 5' end was moved here |
| `tss_count` | pooled READ1 count of the cluster |
| `tss_shift` | bases added (positive) or removed (negative) at the 5' end when snapping to the peak |
| `anchored` | `yes` (assembled transcript starting in the cluster), `derived` (see below), or `no` for unanchored transcripts kept with `--keep_unanchored` (IDs prefixed `U_`) |
| `derived_from` | for derived transcripts: the transcript whose structure was copied |
| `cov`, `FPKM`, `TPM` | StringTie values from the assembly (merged set: from `stringtie --merge`) |

Transcripts that are 3' truncations of another transcript of the same TSS are removed, because CAGEscan 3' ends only
mark where read coverage runs out. A TSS cluster that lies inside an exon of another TSS's transcript, and so has no
assembled transcript (StringTie cannot start a transcript inside continuously covered sequence), gets **derived**
transcripts: copies of the enclosing transcript(s) starting at the cluster peak. See
[usage](usage.md#why-stringtie-and-what-the-pipeline-corrects).

The BED12 follows the original CAGEscan colour scheme (green: plus strand, purple: minus strand). The 5'-most block is
drawn thick when the transcript is TSS-anchored, so `awk '$7 != $8'` keeps only TSS-anchored entries.

These are **5'-anchored transcript fragments**: they extend from the TSS to roughly the end of the sequenced
fragments (the library insert size), not necessarily to the polyadenylation site.

### TSS clusters

<details markdown="1">
<summary>Output files</summary>

- `tss_clusters/`
  - `all_samples.ctss.bed`: CTSS counts pooled across samples (input to paraclu).
  - `all_samples.tss_clusters.bed`: consensus TSS clusters (BED9; `thickStart` is the dominant CTSS, `score` the pooled count capped at 1000).
  - `all_samples.tss_clusters.tsv`: cluster table: `tss_id`, location (1-based), strand, peak (1-based), total count, width.
  - `all_samples.tss_cluster_counts.tsv`: READ1 counts per cluster (rows) and sample (columns).
  - `all_samples.tss_cluster_cpm.tsv`: the same, as counts per million TSS-assigned READ1s.
  - `all_samples.tss_cluster_sl_counts.tsv`: READ1s per cluster whose soft-clipped 5' end matches `--sl_sequence` (only written when SL reads are found). Divide by the counts matrix to get the trans-spliced fraction of each TSS.
  - `tss_anchored_bam/<sample>.tss.bam(.bai)`: read pairs whose READ1 5' end falls in a TSS cluster, tagged with `TC:Z:<tss_id>`. These are the reads used to build transcripts.

</details>

TSS clusters are called with [paraclu](https://gitlab.com/mcfrith/paraclu) on the pooled CTSS
(`--paraclu_min_cluster`, simplified with `paraclu-cut`); sparse tails are trimmed and clusters closer than
`--tss_merge_distance` bases are merged. A READ1 is assigned to a cluster when its 5' end lies in
the cluster extended by `--tss_window` bases on the same strand.

### CAGEscan clusters

<details markdown="1">
<summary>Output files</summary>

- `cagescan_clusters/`
  - `<sample>.cagescan_clusters.bed12`: for each TSS cluster, the union of the aligned blocks of all its read pairs. Score = number of pairs (capped at 1000); thickStart/thickEnd = the TSS cluster.
  - `all_samples.cagescan_clusters.bed12`: meta-clusters combining all samples. Score = number of samples contributing.

</details>

These follow the FANTOM5 CAGEscan conventions ([Bertin et al. 2017](https://doi.org/10.1038/sdata.2017.147),
[CAGEscan-Clustering](https://github.com/nicolas-bertin/CAGEscan-Clustering)). They are a model-free summary of where
the pairs of each TSS land, so they keep low-coverage evidence that StringTie may drop. They are **not
isoform-resolved**: alternative exons are merged into one model, and gaps between READ1 and READ2 that no read covers
look like introns. Use [Transcripts](#transcripts) for isoform structures.

### CTSS

<details markdown="1">
<summary>Output files</summary>

- `ctss/`
  - `<sample>.ctss.bed`: CAGE TSS counts (`chrom, start, end, ., count, strand`; 0-based, one base each).
  - `<sample>.softclip_5p.tsv`: most frequent sequences soft-clipped before the first aligned READ1 base (`-` means no clip). Use this to identify the extra G, leftover linkers or spliced leaders.
  - `bigwig/<sample>.ctss.{plus,minus}.bigWig`: per-strand CTSS signal.

</details>

### Quantification

<details markdown="1">
<summary>Output files</summary>

- `quantification/`
  - `all_samples.transcript_tpm.tsv`, `all_samples.transcript_cov.tsv`: abundance of the merged transcripts per sample (StringTie `-e` on the TSS-anchored BAMs).
  - `<sample>/`: StringTie `-e` output for each sample (`*.quant.transcripts.gtf`, gene (= TSS cluster) abundance, ballgown tables, coverage GTF).

</details>

For TSS-level expression, prefer `tss_clusters/all_samples.tss_cluster_counts.tsv` (one count per capped 5' end).

### gffcompare

<details markdown="1">
<summary>Output files</summary>

- `gffcompare/` (only with `--gtf`)
  - `merged.stats`, `merged.annotated.gtf`, `merged.tracking`, `merged.loci`, `*.tmap`, `*.refmap`: [gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) comparison of the merged transcripts with the reference. Class codes in `*.tmap` tell known (`=`, `c`) from novel isoforms (`j`, `k`, `o`, `u`, ...). Because transcripts are 5' fragments, many match the reference as `c` (contained).

</details>

### Alignment

<details markdown="1">
<summary>Output files</summary>

- `alignment/`
  - `<sample>.Aligned.sortedByCoord.out.bam(.bai)`: STAR alignments (or `<sample>.markdup.bam` with `--dedup`).
  - `star/log/`: STAR logs and splice junctions (`*.SJ.out.tab`).
  - `samtools_stats/`: `samtools stats`, `flagstat` and `idxstats` output.

</details>

### Read QC and trimming

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`: FastQC reports of the raw reads.
- `trimming/fastp/`: fastp reports and logs.
- `trimming/cutadapt/`: cutadapt logs (only with `--r1_5p_linker`).
- `sortmerna/`: SortMeRNA logs (only with `--ribo_database_manifest`).

</details>

### StringTie assemblies

<details markdown="1">
<summary>Output files</summary>

- `stringtie/assembly/`
  - `<sample>.transcripts.gtf`: per-sample StringTie assembly of the TSS-anchored pairs, before anchoring.
  - `merged.gtf`: `stringtie --merge` of the per-sample assemblies, before anchoring.

</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools (FastQC, fastp, cutadapt, SortMeRNA, STAR, samtools) and three pipeline-specific tables: *CTSS extraction* (READ1s used and their 5' soft clips), *TSS cluster assignment* and *Transcript anchoring*. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
