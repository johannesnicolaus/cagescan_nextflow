/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CAT_FASTQ              } from '../modules/nf-core/cat/fastq/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { CUTADAPT               } from '../modules/nf-core/cutadapt/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { SORTMERNA              } from '../modules/nf-core/sortmerna/main'
include { STAR_ALIGN             } from '../modules/nf-core/star/align/main'
include { SAMTOOLS_COLLATE       } from '../modules/nf-core/samtools/collate/main'
include { SAMTOOLS_FIXMATE       } from '../modules/nf-core/samtools/fixmate/main'
include { SAMTOOLS_SORT          } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_MARKDUP       } from '../modules/nf-core/samtools/markdup/main'
include { SAMTOOLS_INDEX         } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS         } from '../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT      } from '../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS      } from '../modules/nf-core/samtools/idxstats/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { PREPARE_GENOME         } from '../subworkflows/local/prepare_genome'
include { TSS_CLUSTERS           } from '../subworkflows/local/tss_clusters'
include { BUILD_TRANSCRIPTS      } from '../subworkflows/local/build_transcripts'
include { COVERAGE_TRACKS        } from '../subworkflows/local/coverage_tracks'
include { FIND_ERNA              } from '../modules/local/find_erna/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_cagescan_pipeline'
include { igvSessionXml          } from '../subworkflows/local/utils_nfcore_cagescan_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CAGESCAN {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // SUBWORKFLOW: Uncompress and index the genome
    //
    PREPARE_GENOME(params.fasta, params.gtf, params.star_index)

    //
    // MODULE: Concatenate FastQ files from the same sample (multiple lanes)
    //
    def ch_fastq = ch_samplesheet
        .map { meta, fastqs -> [meta + [group: meta.group ?: 'all'], fastqs] }
        .branch { _meta, fastqs ->
            single  : fastqs.size() == 2
            multiple: true
        }
    CAT_FASTQ(ch_fastq.multiple)
    def ch_reads = CAT_FASTQ.out.reads.mix(ch_fastq.single)

    //
    // MODULE: Run FastQC on raw reads
    //
    if (!params.skip_fastqc) {
        FASTQC(ch_reads)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { _meta, zip -> zip })
    }

    //
    // MODULE: Remove a 5' linker from READ1 (e.g. nanoCAGE TATAGGG)
    //
    if (params.r1_5p_linker) {
        CUTADAPT(ch_reads)
        ch_reads = CUTADAPT.out.reads
        ch_multiqc_files = ch_multiqc_files.mix(CUTADAPT.out.log.map { _meta, log -> log })
    }

    //
    // MODULE: Adapter and quality trimming
    //
    if (!params.skip_fastp) {
        FASTP(ch_reads.map { meta, reads -> [meta, reads, []] }, false, false, false)
        ch_reads = FASTP.out.reads
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map { _meta, json -> json })
    }

    //
    // MODULE: Optional rRNA removal
    //
    if (params.ribo_database_manifest) {
        def rrna_fastas = file(params.ribo_database_manifest, checkIfExists: true)
            .readLines()
            .findAll { line -> line.trim() }
            .collect { line -> file(line.trim(), checkIfExists: true) }
        SORTMERNA(ch_reads, channel.value([[id: 'rrna'], rrna_fastas]), channel.value([[:], []]))
        ch_reads = SORTMERNA.out.reads
        ch_multiqc_files = ch_multiqc_files.mix(SORTMERNA.out.log.map { _meta, log -> log })
    }

    //
    // MODULE: Spliced alignment of read pairs
    //
    STAR_ALIGN(ch_reads, PREPARE_GENOME.out.star_index, PREPARE_GENOME.out.gtf, !params.gtf)
    ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.map { _meta, log -> log })
    def ch_bam = STAR_ALIGN.out.bam_sorted_aligned

    //
    // MODULES: Optional pair-level duplicate marking (proxy for the missing UMIs)
    //
    def ch_no_fasta = channel.value([[:], [], []])
    if (params.dedup) {
        SAMTOOLS_COLLATE(ch_bam, ch_no_fasta)
        SAMTOOLS_FIXMATE(SAMTOOLS_COLLATE.out.bam, ch_no_fasta)
        SAMTOOLS_SORT(SAMTOOLS_FIXMATE.out.bam, ch_no_fasta, '')
        SAMTOOLS_MARKDUP(SAMTOOLS_SORT.out.bam, ch_no_fasta)
        ch_bam = SAMTOOLS_MARKDUP.out.bam
    }

    SAMTOOLS_INDEX(ch_bam)
    def ch_bam_bai = ch_bam.join(SAMTOOLS_INDEX.out.index, failOnDuplicate: true, failOnMismatch: true)

    SAMTOOLS_STATS(
        ch_bam_bai,
        PREPARE_GENOME.out.fasta.combine(PREPARE_GENOME.out.fai).map { meta, fasta, _fai_meta, fai -> [meta, fasta, fai] }.first()
    )
    SAMTOOLS_FLAGSTAT(ch_bam_bai)
    SAMTOOLS_IDXSTATS(ch_bam_bai)
    ch_multiqc_files = ch_multiqc_files
        .mix(SAMTOOLS_STATS.out.stats.map { _meta, stats -> stats })
        .mix(SAMTOOLS_FLAGSTAT.out.flagstat.map { _meta, flagstat -> flagstat })
        .mix(SAMTOOLS_IDXSTATS.out.idxstats.map { _meta, idxstats -> idxstats })

    //
    // SUBWORKFLOW: CTSS, consensus TSS clusters, group pairs by TSS
    //
    TSS_CLUSTERS(
        ch_bam_bai,
        PREPARE_GENOME.out.sizes,
        params.paraclu_min_cluster,
        !params.skip_sharp_clusters,
        params.reclu_paraclu_min
    )
    ch_multiqc_files = ch_multiqc_files.mix(TSS_CLUSTERS.out.multiqc_files)

    //
    // SUBWORKFLOW: Build, anchor, convert and quantify transcripts
    //
    BUILD_TRANSCRIPTS(
        TSS_CLUSTERS.out.tss_bam,
        TSS_CLUSTERS.out.clusters_tsv,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.gtf,
        (params.stringtie_guide && params.gtf) as boolean,
        params.gtf as boolean
    )
    ch_multiqc_files = ch_multiqc_files.mix(BUILD_TRANSCRIPTS.out.multiqc_files)

    //
    // MODULE: Enhancer-RNA-like divergent TSS pairs, per cluster set
    //
    if (params.find_erna) {
        def ch_erna = TSS_CLUSTERS.out.clusters_tsv
            .join(TSS_CLUSTERS.out.counts.map { meta, files -> [meta.cluster_set, files] })
            .join(BUILD_TRANSCRIPTS.out.gtf.filter { meta, _gtf -> meta.id == 'merged' }.map { meta, gtf -> [meta.cluster_set, gtf] })
            .map { set, tsv, files, gtf -> [[id: 'all_samples', cluster_set: set], tsv, files, gtf] }
        FIND_ERNA(ch_erna, PREPARE_GENOME.out.gtf.map { _meta, gtf -> gtf })
        ch_multiqc_files = ch_multiqc_files.mix(FIND_ERNA.out.stats.map { _meta, stats -> stats })
    }

    //
    // SUBWORKFLOW: Stranded read coverage tracks
    //
    COVERAGE_TRACKS(ch_bam_bai, TSS_CLUSTERS.out.ctss_stats, PREPARE_GENOME.out.sizes)

    //
    // IGV session with genome, gene models, transcripts, TSS clusters and coverage
    //
    ch_bam_bai.map { meta, _bam, _bai -> meta.id }.toSortedList().map { ids -> ids.join(',') }
        .combine(TSS_CLUSTERS.out.clusters_tsv.map { set, _tsv -> set }.toSortedList().map { sets -> sets.join(',') })
        .map { samples, sets -> igvSessionXml(samples.tokenize(','), sets.tokenize(','), params.fasta, params.gtf, params.find_erna) }
        .collectFile(name: 'igv_session.xml', storeDir: "${outdir}")

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'cagescan_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'cagescan'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
