process RECLU_HIERARCHY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/idr:2.0.4.2--py39hcbe4a3b_5' :
        'quay.io/biocontainers/idr:2.0.4.2--py39hcbe4a3b_5' }"

    input:
    tuple val(meta), path(paraclu), path(ctss)

    output:
    tuple val(meta), path("*.reclu.bed"), emit: bed
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    reclu.py hierarchy \\
        --prefix $prefix \\
        --ctss $ctss \\
        $args \\
        $paraclu
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.reclu.bed
    """
}
