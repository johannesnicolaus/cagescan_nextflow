# luscombeu/cagescan: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of luscombeu/cagescan, created with the [nf-core](https://nf-co.re/) template.

### `Added`

Reads and alignment
- Paired-end CAGE without UMIs (CAGEscan) and single-end CAGE (samplesheet without `fastq_2`; transcripts skipped).
- Random subsampling per sample (`--subsample_reads`, default 10,000,000; seqtk), optional READ1 5' linker removal
  (cutadapt), fastp, optional rRNA removal (SortMeRNA).
- Spliced leaders (`--sl_sequence`, one or several named leaders): reads split into SL and non-SL sets before
  alignment, the leader trimmed exactly (cutadapt; leader at READ1 base 1 or after one G, mismatches up to
  `--sl_max_error_rate`, no indels), reads renamed `SL.<leader>:` / `SLG.<leader>:`; following splitSpliceLeaderPE.
- STAR alignment with READ1 aligned from its first base (`--read1_end_to_end`) and removal of a non-templated first G
  that differs from the genome (as CAGEr's removeFirstG); `--star_extra_args`.
- Optional pair-level duplicate marking (`--dedup`, samtools markdup; paired-end only).

TSSs and trans-splice sites
- CTSS from READ1 5' ends; SL-trimmed READ1s start at the trans-splice site (up to 3 soft-clipped 5' bases added back,
  longer clips excluded); with local alignment, SL recognition in the 5' soft clip and the shift to the acceptor when
  the leader's last bases align.
- Broad TSS clusters with CAGEfightR (default: per-sample TPM relative to all READ1s, CTSSs in >= 2 samples,
  slice-reduce with `--tss_merge_distance` 20, >= 5 TPM in >= 1 sample) or paraclu (`--cluster_method paraclu`);
  trans-splice sites (`SL…`) clustered separately from TSSs (`TC…`, `--sl_separate_clusters`).
- Sharp TSS clusters following RECLU (per-replicate paraclu hierarchy, TPM-per-base filter, hierarchical
  (log-)stability, reciprocal-overlap pairing and IDR); optional `group` samplesheet column for replicates.
- TSS cluster counts, CPM, spliced-leader counts and (several leaders) per-leader counts; TSS-anchored BAMs.
- MultiQC: CTSS extraction, PyPu initiator at TSSs, AG acceptor at trans-splice sites, spliced-leader shares, TSS
  assignment.

Transcripts (paired-end)
- Read pairs grouped by TSS cluster (proper pairs), FANTOM5-style CAGEscan clusters and meta-clusters, StringTie
  assembly per sample and merged, anchoring of 5' ends to cluster peaks, derived transcripts for TSSs inside exons,
  removal of 3'-truncated copies; GTF, GFF3, BED12, FASTA (gffread); quantification (StringTie `-e`); gffcompare.

Other
- eRNA candidates: divergent, distal, balanced (and, paired-end, unspliced) pairs of non-SL TSS clusters, with
  convergent pairs as a control.
- Stranded coverage bigWigs, pooled TSS-signal bigWigs and a self-contained IGV session.
- Reference GTFs: records with an empty `transcript_id` (gene rows of NCBI-style GTFs, which gffcompare rejects) are
  dropped before use.
- Simulated test data (`tests/data/make_test_data.py`) with known TSSs, isoforms and a trans-spliced gene; tests for
  paired-end, single-end, paraclu, NCBI-style GTFs and every local module.

### `Fixed`

- Patched nf-core paraclu module: CTSS are sorted by strand (plus and minus strand CTSS were interleaved and
  fragmented paraclu's per-strand runs).
- Numeric and boolean options given on the command line are typed (`params` block in `main.nf`): with Nextflow 26.04
  they were strings (`--sl_split false` was truthy). Requires Nextflow >= 26.04.

### `Dependencies`

### `Deprecated`
