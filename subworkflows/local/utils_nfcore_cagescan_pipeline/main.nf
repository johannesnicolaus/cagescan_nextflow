//
// Subworkflow with functionality specific to the luscombeu/cagescan pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        true // cast command-line values to the schema types (Nextflow >= 26.04 passes them as strings)
    )

    validateSplicedLeaders(params.sl_sequence)

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //

    channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Spliced leaders 'SEQ' or 'NAME=SEQ,...': names (given, or SL / SL1, SL2, ...) must be unique
//
def validateSplicedLeaders(spec) {
    if (!spec) {
        return
    }
    def entries = spec.tokenize(',')*.trim().findAll { entry -> entry }
    def names = entries.withIndex().collect { entry, i ->
        entry.contains('=') ? entry.tokenize('=')[0] : (entries.size() == 1 ? 'SL' : "SL${i + 1}")
    }
    def duplicates = names.findAll { name -> names.count(name) > 1 }.unique()
    if (duplicates) {
        error("--sl_sequence: spliced-leader names must be unique (duplicated: ${duplicates.join(', ')}); name every leader, e.g. SL1=...,SL2=...")
    }
}

//
// Single-end CAGE: every samplesheet row without fastq_2. Mixing single-end and paired-end samples is not supported.
//
def isSingleEnd(input) {
    def rows = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
    def kinds = rows.collect { row -> row[2] ? 'paired' : 'single' }.unique()
    if (kinds.size() > 1) {
        error("Please check input samplesheet -> all samples must be single-end (no fastq_2) or all paired-end")
    }
    return kinds == ['single']
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Generate methods description for MultiQC
//
//
// Tools used in this run, for the methods text: [condition, citation in the text, bibliography entry]
//
def usedTools() {
    def single_end = params.input ? isSingleEnd(params.input) : false
    def sl_split = params.sl_sequence && params.sl_split
    def end_to_end = params.read1_end_to_end && !(params.sl_sequence && !params.sl_split)
    return [
        [!params.skip_fastqc, "FastQC (Andrews 2010)",
            "Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/"],
        [params.subsample_reads as boolean, "seqtk (Li)",
            "Li H. seqtk: toolkit for processing sequences in FASTA/Q formats. https://github.com/lh3/seqtk"],
        [params.r1_5p_linker || sl_split, "cutadapt (Martin 2011)",
            "Martin M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10-12. doi: 10.14806/ej.17.1.200"],
        [sl_split, "spliced-leader splitting as in splitSpliceLeaderPE (Plessy)",
            "Plessy C. splitSpliceLeaderPE: split paired-end reads by spliced leader. https://github.com/oist/plessy_splitspliceleaderpe"],
        [!params.skip_fastp, "fastp (Chen et al. 2018)",
            "Chen S, Zhou Y, Chen Y, Gu J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890. doi: 10.1093/bioinformatics/bty560"],
        [params.ribo_database_manifest as boolean, "SortMeRNA (Kopylova et al. 2012)",
            "Kopylova E, Noé L, Touzet H. (2012). SortMeRNA: fast and accurate filtering of ribosomal RNAs in metatranscriptomic data. Bioinformatics, 28(24), 3211-3217. doi: 10.1093/bioinformatics/bts611"],
        [true, "STAR (Dobin et al. 2013)",
            "Dobin A, Davis CA, Schlesinger F, et al. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 29(1), 15-21. doi: 10.1093/bioinformatics/bts635"],
        [true, "SAMtools (Danecek et al. 2021)",
            "Danecek P, Bonfield JK, Liddle J, et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. doi: 10.1093/gigascience/giab008"],
        [true, "pysam",
            "Heger A, Jacobs K, et al. pysam: a Python module for reading and manipulating SAM/BAM files. https://github.com/pysam-developers/pysam"],
        [end_to_end, "first-G removal as in CAGEr (Haberle et al. 2015)",
            "Haberle V, Forrest ARR, Hayashizaki Y, Carninci P, Lenhard B. (2015). CAGEr: precise TSS data retrieval and high-resolution promoterome mining for integrative analyses. Nucleic Acids Research, 43(8), e51. doi: 10.1093/nar/gkv054"],
        [params.cluster_method == 'cagefightr', "CAGEfightR (Thodberg et al. 2019)",
            "Thodberg M, Thieffry A, Vitting-Seerup K, Andersson R, Sandelin A. (2019). CAGEfightR: analysis of 5'-end data using R/Bioconductor. BMC Bioinformatics, 20, 487. doi: 10.1186/s12859-019-3029-5"],
        [params.cluster_method == 'paraclu' || !params.skip_sharp_clusters, "paraclu (Frith et al. 2008)",
            "Frith MC, Valen E, Krogh A, Hayashizaki Y, Carninci P, Sandelin A. (2008). A code for transcription initiation in mammalian genomes. Genome Research, 18(1), 1-12. doi: 10.1101/gr.6831208"],
        [!params.skip_sharp_clusters, "RECLU (Ohmiya et al. 2014)",
            "Ohmiya H, Vitezic M, Frith MC, et al. (2014). RECLU: a pipeline to discover reproducible transcriptional start sites and their alternative regulation using capped analysis of gene expression (CAGE). BMC Genomics, 15, 269. doi: 10.1186/1471-2164-15-269"],
        [!params.skip_sharp_clusters, "IDR (Li et al. 2011)",
            "Li Q, Brown JB, Huang H, Bickel PJ. (2011). Measuring reproducibility of high-throughput experiments. Annals of Applied Statistics, 5(3), 1752-1779. doi: 10.1214/11-AOAS466"],
        [!single_end, "StringTie (Pertea et al. 2015)",
            "Pertea M, Pertea GM, Antonescu CM, et al. (2015). StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nature Biotechnology, 33(3), 290-295. doi: 10.1038/nbt.3122"],
        [!single_end, (params.gtf ? "GffRead and GffCompare" : "GffRead") + " (Pertea and Pertea 2020)",
            "Pertea G, Pertea M. (2020). GFF Utilities: GffRead and GffCompare. F1000Research, 9, 304. doi: 10.12688/f1000research.23297.2"],
        [true, "BEDTools (Quinlan and Hall 2010)",
            "Quinlan AR, Hall IM. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), 841-842. doi: 10.1093/bioinformatics/btq033"],
        [true, "bedGraphToBigWig (Kent et al. 2010)",
            "Kent WJ, Zweig AS, Barber G, Hinrichs AS, Karolchik D. (2010). BigWig and BigBed: enabling browsing of large distributed datasets. Bioinformatics, 26(17), 2204-2207. doi: 10.1093/bioinformatics/btq351"],
        [true, "IGV (Robinson et al. 2011)",
            "Robinson JT, Thorvaldsdóttir H, Winckler W, et al. (2011). Integrative genomics viewer. Nature Biotechnology, 29(1), 24-26. doi: 10.1038/nbt.1754"],
        [!single_end, "the CAGEscan approach (Plessy et al. 2010; Bertin et al. 2017)",
            "Plessy C, Bertin N, Takahashi H, et al. (2010). Linking promoters to functional transcripts in small samples with nanoCAGE and CAGEscan. Nature Methods, 7(7), 528-534. doi: 10.1038/nmeth.1470"],
        [!single_end, "",
            "Bertin N, Mendez M, Hasegawa A, et al. (2017). Linking FANTOM5 CAGE peaks to annotations with CAGEscan. Scientific Data, 4, 170147. doi: 10.1038/sdata.2017.147"],
        [params.find_erna as boolean, "eRNA-like divergent TSS pairs as defined by Andersson et al. (2014)",
            "Andersson R, Gebhard C, Miguel-Escalada I, et al. (2014). An atlas of active enhancers across human cell types and tissues. Nature, 507, 455-461. doi: 10.1038/nature12787"],
        [true, "MultiQC (Ewels et al. 2016)",
            "Ewels P, Magnusson M, Lundin S, Käller M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047-3048. doi: 10.1093/bioinformatics/btw354"],
    ].findAll { used, _text, _reference -> used }
}

def toolCitationText() {
    def cited = usedTools().collect { _used, text, _reference -> text }.findAll { text -> text }
    return "Tools used in the workflow included: " + cited.join(', ') + "."
}

def toolBibliographyText() {
    return usedTools().collect { _used, _text, reference -> "<li>${reference}</li>" }.join(' ')
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
