process CTSS_MOTIFS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0' :
        'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0' }"

    input:
    tuple val(meta), path(ctss), path(sl_ctss) // ctss: non-SL CTSS when sl_ctss is given; sl_ctss optional ([])
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.tss_initiator_mqc.tsv")    , emit: initiator
    tuple val(meta), path("*.sl_acceptor_mqc.tsv")      , emit: acceptor, optional: true
    tuple val(meta), path("*.ctss_motif_stats_mqc.tsv") , emit: stats
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), topic: versions, emit: versions_pysam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ctss_motifs.py \\
        --fasta $fasta \\
        --prefix $prefix \\
        --ctss $ctss \\
        ${sl_ctss ? "--sl-ctss $sl_ctss" : ''} \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tss_initiator_mqc.tsv ${prefix}.ctss_motif_stats_mqc.tsv
    ${sl_ctss ? "touch ${prefix}.sl_acceptor_mqc.tsv" : ''}
    """
}
