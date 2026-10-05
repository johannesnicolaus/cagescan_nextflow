process ANCHOR_TRANSCRIPTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0' :
        'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0' }"

    input:
    tuple val(meta), path(gtf, stageAs: 'input/*')
    path(clusters)  // tss_clusters.tsv

    output:
    tuple val(meta), path("*.transcripts.gtf")    , emit: gtf
    tuple val(meta), path("*.transcripts.bed12")  , emit: bed
    tuple val(meta), path("*.anchor_stats_mqc.tsv"), emit: stats
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    anchor_transcripts.py \\
        $gtf \\
        $clusters \\
        --prefix $prefix \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.transcripts.gtf ${prefix}.transcripts.bed12 ${prefix}.anchor_stats_mqc.tsv
    """
}
