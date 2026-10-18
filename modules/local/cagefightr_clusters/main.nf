process CAGEFIGHTR_CLUSTERS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-cagefightr:1.26.0--r44hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-cagefightr:1.26.0--r44hdfd78af_0' }"

    input:
    tuple val(meta), path(ctss, stageAs: 'ctss/*') // per-sample CTSS beds (one class: all, SL or non-SL READ1s)
    tuple val(meta2), path(fai)
    path(library_sizes) // TSV: sample, READ1s (TPM denominators), or []

    output:
    tuple val(meta), path("*.cagefightr.bed"), emit: bed
    tuple val("${task.process}"), val('CAGEfightR'), eval("Rscript -e 'cat(as.character(packageVersion(\"CAGEfightR\")))'"), topic: versions, emit: versions_cagefightr

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cagefightr_clusters.R \\
        --prefix $prefix \\
        --fai $fai \\
        ${library_sizes ? "--library-sizes $library_sizes" : ''} \\
        $args \\
        $ctss
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cagefightr.bed
    """
}
