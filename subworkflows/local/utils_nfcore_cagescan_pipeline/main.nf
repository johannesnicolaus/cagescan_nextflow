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
        false
    )

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
def toolCitationText() {
    def citation_text = [
            "Tools used in the workflow included:",
            params.skip_fastqc ? "" : "FastQC (Andrews 2010),",
            params.r1_5p_linker ? "cutadapt (Martin 2011)," : "",
            params.skip_fastp ? "" : "fastp (Chen et al. 2018),",
            params.ribo_database_manifest ? "SortMeRNA (Kopylova et al. 2012)," : "",
            "STAR (Dobin et al. 2013),",
            "SAMtools (Danecek et al. 2021),",
            "pysam,",
            "paraclu (Frith et al. 2008),",
            params.skip_sharp_clusters ? "" : "RECLU (Ohmiya et al. 2014) with IDR (Li et al. 2011),",
            "StringTie (Pertea et al. 2015),",
            "GffRead and GffCompare (Pertea and Pertea 2020),",
            "bedGraphToBigWig (Kent et al. 2010),",
            "following the CAGEscan approach (Plessy et al. 2010; Bertin et al. 2017),",
            params.find_erna ? "eRNA-like divergent TSS pairs as defined by Andersson et al. (2014)," : "",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].findAll { text -> text }.join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
            params.skip_fastqc ? "" : "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            params.r1_5p_linker ? "<li>Martin M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10-12. doi: 10.14806/ej.17.1.200</li>" : "",
            params.skip_fastp ? "" : "<li>Chen S, Zhou Y, Chen Y, Gu J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890. doi: 10.1093/bioinformatics/bty560</li>",
            params.ribo_database_manifest ? "<li>Kopylova E, Noé L, Touzet H. (2012). SortMeRNA: fast and accurate filtering of ribosomal RNAs in metatranscriptomic data. Bioinformatics, 28(24), 3211-3217. doi: 10.1093/bioinformatics/bts611</li>" : "",
            "<li>Dobin A, Davis CA, Schlesinger F, et al. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 29(1), 15-21. doi: 10.1093/bioinformatics/bts635</li>",
            "<li>Danecek P, Bonfield JK, Liddle J, et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. doi: 10.1093/gigascience/giab008</li>",
            "<li>Frith MC, Valen E, Krogh A, Hayashizaki Y, Carninci P, Sandelin A. (2008). A code for transcription initiation in mammalian genomes. Genome Research, 18(1), 1-12. doi: 10.1101/gr.6831208</li>",
            params.skip_sharp_clusters ? "" : "<li>Ohmiya H, Vitezic M, Frith MC, et al. (2014). RECLU: a pipeline to discover reproducible transcriptional start sites and their alternative regulation using capped analysis of gene expression (CAGE). BMC Genomics, 15, 269. doi: 10.1186/1471-2164-15-269</li>",
            params.skip_sharp_clusters ? "" : "<li>Li Q, Brown JB, Huang H, Bickel PJ. (2011). Measuring reproducibility of high-throughput experiments. Annals of Applied Statistics, 5(3), 1752-1779. doi: 10.1214/11-AOAS466</li>",
            "<li>Pertea M, Pertea GM, Antonescu CM, et al. (2015). StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nature Biotechnology, 33(3), 290-295. doi: 10.1038/nbt.3122</li>",
            "<li>Pertea G, Pertea M. (2020). GFF Utilities: GffRead and GffCompare. F1000Research, 9, 304. doi: 10.12688/f1000research.23297.2</li>",
            "<li>Kent WJ, Zweig AS, Barber G, Hinrichs AS, Karolchik D. (2010). BigWig and BigBed: enabling browsing of large distributed datasets. Bioinformatics, 26(17), 2204-2207. doi: 10.1093/bioinformatics/btq351</li>",
            "<li>Plessy C, Bertin N, Takahashi H, et al. (2010). Linking promoters to functional transcripts in small samples with nanoCAGE and CAGEscan. Nature Methods, 7(7), 528-534. doi: 10.1038/nmeth.1470</li>",
            "<li>Bertin N, Mendez M, Hasegawa A, et al. (2017). Linking FANTOM5 CAGE peaks to annotations with CAGEscan. Scientific Data, 4, 170147. doi: 10.1038/sdata.2017.147</li>",
            params.find_erna ? "<li>Andersson R, Gebhard C, Miguel-Escalada I, et al. (2014). An atlas of active enhancers across human cell types and tissues. Nature, 507, 455-461. doi: 10.1038/nature12787</li>" : "",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].findAll { text -> text }.join(' ').trim()

    return reference_text
}

//
// IGV session (relative paths inside --outdir; genome and annotation by absolute path)
//
def igvSessionXml(samples, sets, fasta, gtf, erna) {
    def tracks = []
    if (gtf) {
        tracks << [file(gtf).toString(), 'Gene models (--gtf)', 'displayMode="EXPANDED"']
    }
    sets.each { set ->
        tracks << ["transcripts/${set}/merged/merged.${set}.transcripts.bed12", "Transcripts (${set} TSS clusters)", 'displayMode="EXPANDED"']
        tracks << ["tss_clusters/${set}/all_samples.${set}.tss_clusters.bed", "TSS clusters (${set})", 'displayMode="COLLAPSED"']
        if (erna) {
            tracks << ["erna/${set}/all_samples.${set}.erna_candidates.bed", "eRNA-like divergent pairs (${set})", 'displayMode="EXPANDED" color="120,0,160"']
        }
    }
    tracks << ['ctss/bigwig/all_samples.ctss.plus.bigWig', 'TSS signal + (all samples)', 'color="200,0,0" autoScale="true" autoscaleGroup="ctss"']
    tracks << ['ctss/bigwig/all_samples.ctss.minus.bigWig', 'TSS signal - (all samples)', 'color="0,0,200" autoScale="true" autoscaleGroup="ctss"']
    samples.each { sample ->
        tracks << ["coverage/${sample}.coverage.plus.bigWig", "${sample} coverage +", "color=\"200,0,0\" autoScale=\"true\" autoscaleGroup=\"${sample}\""]
        tracks << ["coverage/${sample}.coverage.minus.bigWig", "${sample} coverage -", "color=\"0,0,200\" autoScale=\"true\" autoscaleGroup=\"${sample}\""]
    }
    def xml = new StringBuilder()
    xml << '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
    xml << "<Session genome=\"${file(fasta).toString()}\" locus=\"All\" relativePath=\"true\" version=\"8\">\n"
    xml << '    <Resources>\n'
    tracks.each { t -> xml << "        <Resource path=\"${t[0]}\"/>\n" }
    xml << '    </Resources>\n'
    xml << '    <Panel name="DataPanel">\n'
    tracks.each { t -> xml << "        <Track id=\"${t[0]}\" name=\"${t[1]}\" ${t[2]}/>\n" }
    xml << '    </Panel>\n'
    xml << '</Session>\n'
    return xml.toString()
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
