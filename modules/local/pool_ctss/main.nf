process POOL_CTSS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0' :
        'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0' }"

    input:
    tuple val(meta), path(ctss, stageAs: 'input/*')

    output:
    tuple val(meta), path("*.ctss.bed")                   , emit: ctss
    tuple val(meta), path("*.ctss.{plus,minus}.bedgraph"), emit: bedgraph
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tss_clusters.py pool \\
        --prefix $prefix \\
        $ctss
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ctss.bed ${prefix}.ctss.plus.bedgraph ${prefix}.ctss.minus.bedgraph
    """
}
