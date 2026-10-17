#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    luscombeu/cagescan
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/luscombeu/cagescan
----------------------------------------------------------------------------------------
*/

//
// Typed declarations of the boolean and numeric parameters (defaults: nextflow.config). From Nextflow 26.04
// command-line values arrive as strings unless declared with a type here ('--flag false' would be truthy).
//
params {
    subsample_reads          : Integer = 10000000
    subsample_seed           : Integer = 100
    linker_discard_untrimmed : Boolean = true
    r1_trim_front            : Integer = 0
    r2_trim_front            : Integer = 0
    skip_fastqc              : Boolean = false
    skip_fastp               : Boolean = false
    min_mapq                 : Integer = 10
    dedup                    : Boolean = false
    sl_split                 : Boolean = true
    save_sl_split_reads      : Boolean = true
    sl_max_error_rate        : Float   = 0.15
    paraclu_min_cluster      : Integer = 30
    tss_window               : Integer = 50
    tss_merge_distance       : Integer = 20
    skip_sharp_clusters      : Boolean = false
    reclu_paraclu_min        : Integer = 2
    reclu_min_tpm_per_base   : Float   = 0.1
    reclu_min_overlap        : Float   = 0.9
    reclu_max_idr            : Float   = 0.1
    reclu_max_length         : Integer = 200
    reclu_min_pairs          : Integer = 50
    min_transcript_length    : Integer = 100
    stringtie_guide          : Boolean = false
    keep_unanchored          : Boolean = false
    derive_transcripts       : Boolean = true
    version                  : Boolean = false
    plaintext_email          : Boolean = false
    monochrome_logs          : Boolean = false
    validate_params          : Boolean = true
    help_full                : Boolean = false
    show_hidden              : Boolean = false
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CAGESCAN  } from './workflows/cagescan'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_cagescan_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_cagescan_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow LUSCOMBEU_CAGESCAN {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    CAGESCAN (
        samplesheet,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )
    emit:
    multiqc_report = CAGESCAN.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    LUSCOMBEU_CAGESCAN (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        LUSCOMBEU_CAGESCAN.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
