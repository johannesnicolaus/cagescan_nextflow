process ASSIGN_PAIRS_TO_TSS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0' :
        'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(clusters)

    output:
    tuple val(meta), path("*.tss.bam"), path("*.tss.bam.bai"), emit: bam
    tuple val(meta), path("*.tss_counts.tsv")                , emit: counts
    tuple val(meta), path("*.cagescan_clusters.bed12")       , emit: cagescan_clusters
    tuple val(meta), path("*.tss_assign_stats_mqc.tsv")      , emit: stats
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), topic: versions, emit: versions_pysam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    assign_pairs_to_tss.py \\
        $bam \\
        $clusters \\
        --prefix $prefix \\
        --threads $task.cpus \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tss.bam ${prefix}.tss.bam.bai ${prefix}.tss_counts.tsv ${prefix}.cagescan_clusters.bed12 ${prefix}.tss_assign_stats_mqc.tsv
    """
}
