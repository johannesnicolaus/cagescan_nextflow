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
- CTSS of spliced-leader READ1s whose leader end aligned to the splice acceptor are moved to the trans-splice site (up to 3 bases downstream).
- Fixed the nf-core paraclu module sorting CTSS without strand (patched; plus and minus strand CTSS interleaved and fragmented paraclu's per-strand runs).
- Stranded coverage bigWigs per sample, pooled TSS-signal bigWigs, and an `igv_session.xml` for genome browsing.
- eRNA candidates: divergent, distal, balanced, unspliced pairs of non-SL TSS clusters per cluster set, with convergent pairs as a control (`erna/`, MultiQC table, IGV track).
- Read pairs split into SL and non-SL sets before alignment, by the spliced leader at the READ1 start, allowing mismatches, with the SL trimmed (cutadapt; `--sl_split`, `--sl_max_error_rate`, `--save_sl_split_reads`), following splitSpliceLeaderPE; recovers 15-20% more SL reads on *O. dioica* CAGE than the soft-clip match.
- READ1 aligned to its 5' end by STAR (`--read1_end_to_end`, `--alignEndsType Extend5pOfRead1`) and a non-templated first G that differs from the genome removed (as CAGEr's removeFirstG); `--star_extra_args`.
- Broad TSS clusters called with CAGEfightR by default (TPM, CTSS support in >= 2 samples, slice-reduce clustering, >= 5 TPM in >= 1 sample as CAGEr's consensus clusters; `--cluster_method`, `--ctss_min_samples`, `--tc_min_tpm`, `--tc_min_samples`); paraclu remains available.
- Several spliced leaders (`--sl_sequence NAME=SEQ,NAME=SEQ`, e.g. C. elegans SL1 and SL2): all looked for in every read; SL READ1s per leader per cluster (`tss_cluster_sl_leaders.tsv`) and per sample (MultiQC). SL helpers shared in `bin/cagescan_sl.py`.
- MultiQC: initiator dinucleotide (-1/+1) of TSS READ1s and the two bases before trans-splice sites, per sample (`CTSS_MOTIFS`).
- SL split: exact leader boundary (two anchored patterns, with and without the extra G; mismatches but no indels), and up to 3 bases soft-clipped at the 5' end of SL-trimmed READ1s are added back to their CTSS (first exon bases that differ from the assembly no longer move trans-splice sites downstream); SL-trimmed READ1s with longer 5' clips give no TSS.
- Spliced-leader READ1s (trans-splice sites) and the other READ1s are clustered separately in the broad set (`SL…` and `TC…` clusters; `--sl_separate_clusters`), so that a TSS next to a trans-splice site keeps its own peak.
- Random subsampling of read pairs per sample with seqtk, `--subsample_reads` (default 10,000,000; 0 = all reads) and `--subsample_seed`.
- Simulated test data set (`tests/data/make_test_data.py`) with known TSSs and isoforms.

### `Fixed`

- Integer and number options given on the command line (e.g. `--min_mapq 5`) failed parameter validation with Nextflow 26.04, which passes them as strings: the boolean and numeric parameters are now declared with types in a `params` block in `main.nf`, so Nextflow converts command-line values (`--sl_split false` was the truthy string `'false'`). Requires Nextflow >= 26.04.

### `Dependencies`

### `Deprecated`
