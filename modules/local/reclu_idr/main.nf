process RECLU_IDR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/idr:2.0.4.2--py39hcbe4a3b_5' :
        'quay.io/biocontainers/idr:2.0.4.2--py39hcbe4a3b_5' }"

    input:
    tuple val(meta), path(rep1, stageAs: 'rep1/*'), path(rep2, stageAs: 'rep2/*')

    output:
    tuple val(meta), path("*.reclu_idr.tsv"), emit: tsv
    tuple val("${task.process}"), val('idr'), eval("python -c 'import idr; print(idr.__version__)'"), topic: versions, emit: versions_idr

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    reclu.py idr \\
        --prefix $prefix \\
        $args \\
        $rep1 \\
        $rep2
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.reclu_idr.tsv
    """
}
