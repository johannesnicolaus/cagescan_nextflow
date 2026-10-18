process COMBINE_TSS_CLUSTERS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0' :
        'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0' }"

    input:
    tuple val(meta), path(clusters), path(sl_clusters) // sl_clusters optional ([])

    output:
    tuple val(meta), path("*.tss_clusters.bed"), emit: bed
    tuple val(meta), path("*.tss_clusters.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tss_clusters.py combine \\
        --prefix $prefix \\
        ${sl_clusters ? "--sl-clusters $sl_clusters" : ''} \\
        $args \\
        $clusters
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tss_clusters.bed
    echo -e "tss_id\\tlocation\\tstrand\\tpeak\\ttotal_count\\twidth" > ${prefix}.tss_clusters.tsv
    """
}
