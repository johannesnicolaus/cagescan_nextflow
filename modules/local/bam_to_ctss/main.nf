process BAM_TO_CTSS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0' :
        'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.ctss.bed")              , emit: ctss
    tuple val(meta), path("*.ctss.{plus,minus}.bedgraph"), emit: bedgraph
    tuple val(meta), path("*.softclip_5p.tsv")       , emit: softclip
    tuple val(meta), path("*.ctss_stats_mqc.tsv")    , emit: stats
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), topic: versions, emit: versions_pysam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bam_to_ctss.py \\
        $bam \\
        --prefix $prefix \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ctss.bed ${prefix}.ctss.plus.bedgraph ${prefix}.ctss.minus.bedgraph
    touch ${prefix}.softclip_5p.tsv ${prefix}.ctss_stats_mqc.tsv
    """
}
