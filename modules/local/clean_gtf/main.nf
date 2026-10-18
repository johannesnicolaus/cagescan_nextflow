process CLEAN_GTF {
    tag "$gtf"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0' :
        'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.clean.gtf"), emit: gtf
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    clean_gtf.py $gtf ${gtf.baseName}.clean.gtf
    """

    stub:
    """
    touch ${gtf.baseName}.clean.gtf
    """
}
