# luscombeu/cagescan


[![GitHub Actions CI Status](https://github.com/luscombeu/cagescan/actions/workflows/nf-test.yml/badge.svg)](https://github.com/luscombeu/cagescan/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/luscombeu/cagescan/actions/workflows/linting.yml/badge.svg)](https://github.com/luscombeu/cagescan/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/luscombeu/cagescan)

## Introduction

**luscombeu/cagescan** reconstructs 5'-anchored transcripts from **paired-end CAGE reads without UMIs**
(CAGEscan-style libraries: READ1 starts at the capped 5' end, READ2 is random-primed downstream).

It is a reference-guided re-implementation of the idea behind the original
[CAGEscan pipeline](https://gitlab.com/mcfrith/cagescan-pipeline). That pipeline groups reads into molecules by UMI
and assembles each molecule de novo. Without UMIs, the closest equivalent of a molecule is a **TSS cluster**, so this
pipeline maps the read pairs first, groups them by the TSS cluster that READ1 starts in, and builds one or more
isoforms per TSS from the spliced pair alignments.

1. Merge lanes ([`cat`](https://www.gnu.org/software/coreutils/)) and read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Optional READ1 5' linker removal ([`cutadapt`](https://cutadapt.readthedocs.io/)), adapter/quality trimming ([`fastp`](https://github.com/OpenGene/fastp)) and optional rRNA removal ([`SortMeRNA`](https://github.com/sortmerna/sortmerna))
3. Spliced paired alignment ([`STAR`](https://github.com/alexdobin/STAR)); local alignment soft-clips the non-templated 5' G or a spliced leader on READ1
4. Optional pair-level duplicate marking ([`samtools markdup`](http://www.htslib.org/)) and alignment QC ([`samtools`](http://www.htslib.org/))
5. CTSS extraction from READ1 5' ends, with a report of what was soft-clipped before the cap (extra G, linker, spliced leader)
6. Two TSS cluster sets: **broad** clusters from CTSS pooled across samples ([`paraclu`](https://gitlab.com/mcfrith/paraclu)), and **sharp** clusters reproducible across replicates, following [RECLU](https://doi.org/10.1186/1471-2164-15-269) (paraclu hierarchy + [IDR](https://github.com/nboley/idr)); CTSS bigWigs ([`bedGraphToBigWig`](https://genome.ucsc.edu/goldenPath/help/bigWig.html)). All later steps run for both sets, so you can choose afterwards
7. Grouping of read pairs by TSS cluster (READ1 5' end inside a cluster), per-sample TSS cluster counts, and FANTOM5-style CAGEscan clusters / meta-clusters ([Bertin et al. 2017](https://doi.org/10.1038/sdata.2017.147))
8. Isoform assembly from the TSS-anchored pairs per sample, merged across samples ([`StringTie`](https://ccb.jhu.edu/software/stringtie/))
9. Anchoring: every transcript's 5' end is moved to its TSS cluster's dominant CTSS; transcripts are named `<tss_id>.<n>`
10. Transcript GTF, GFF3, BED12 and FASTA ([`gffread`](https://github.com/gpertea/gffread)), per-sample quantification (StringTie `-e`) and comparison with a reference annotation ([`gffcompare`](https://github.com/gpertea/gffcompare))
11. Report ([`MultiQC`](http://multiqc.info/))

> [!NOTE]
> CAGEscan transcripts span from the TSS to the end of the sequenced fragment (roughly the insert size, a few hundred
> bases to ~1 kb). They are **5'-anchored transcript fragments**, not full-length transcripts.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

Prepare a samplesheet with one row per pair of FASTQ files (rows with the same `sample` are treated as lanes and concatenated):

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,group
CAGE_REP1,CAGE_REP1_R1.fastq.gz,CAGE_REP1_R2.fastq.gz,ctrl
CAGE_REP2,CAGE_REP2_R1.fastq.gz,CAGE_REP2_R2.fastq.gz,ctrl
```

`fastq_1` must be READ1, i.e. the read that starts at the capped 5' end. `group` (optional) names the condition each
sample is a replicate of; sharp TSS clusters need at least two replicates in a group.

Now, you can run the pipeline using:

```bash
nextflow run luscombeu/cagescan \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --fasta genome.fa \
   --gtf annotation.gtf \
   --outdir <OUTDIR>
```

`--gtf` is optional. If you do not know the READ1 structure, run once without `--r1_5p_linker` and look at
`ctss/<sample>.softclip_5p.tsv` (and the *CTSS extraction* table in the MultiQC report): a single `G` is the
expected non-templated base, a recurring longer sequence is a linker (`--r1_5p_linker`) or a spliced leader
(`--sl_sequence`). See [usage](docs/usage.md) and [output](docs/output.md) for details.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Credits

luscombeu/cagescan was originally written by Johannes Nicolaus Wibisana.

We thank the following people for their extensive assistance in the development of this pipeline:

- Martin C. Frith, for the original [CAGEscan pipeline](https://gitlab.com/mcfrith/cagescan-pipeline) whose BED conventions are reused here.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use luscombeu/cagescan for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->


An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
