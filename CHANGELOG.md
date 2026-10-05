# luscombeu/cagescan: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of luscombeu/cagescan, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Reference-guided CAGEscan workflow for paired-end CAGE without UMIs: STAR alignment, CTSS extraction with READ1 5' soft-clip report (extra G / linker / spliced leader), consensus paraclu TSS clusters, grouping of read pairs by TSS cluster, StringTie isoform assembly per sample and merged, TSS anchoring, GTF / GFF3 / BED12 / FASTA output, TSS cluster and transcript count matrices, gffcompare, MultiQC.
- FANTOM5-style CAGEscan clusters per sample and meta-clusters across samples; `--r2_trim_front` for random-primer bases on READ2; proper-pair requirement for TSS assignment (following Bertin et al. 2017).
- TSS clusters closer than `--tss_merge_distance` are merged; transcripts anchor to the strongest cluster in the window; 3'-truncated duplicates are dropped; derived transcripts for TSS clusters inside exons of other transcripts (`--derive_transcripts`).
- Sharp TSS clusters following RECLU (per-replicate paraclu hierarchy, TPM-per-base filter, hierarchical (log-)stability, reciprocal-overlap pairing and IDR) next to the broad pooled clusters; everything downstream runs for both sets. Optional `group` samplesheet column for replicates.
- Fixed the nf-core paraclu module sorting CTSS without strand (patched; plus and minus strand CTSS interleaved and fragmented paraclu's per-strand runs).
- Simulated test data set (`tests/data/make_test_data.py`) with known TSSs and isoforms.

### `Fixed`

### `Dependencies`

### `Deprecated`
