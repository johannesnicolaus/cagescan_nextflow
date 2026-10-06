process PARACLU_RAW {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/paraclu:10--h9a82719_1' :
        'quay.io/biocontainers/paraclu:10--h9a82719_1' }"

    input:
    tuple val(meta), path(ctss)  // chrom, pos0, pos0+1, ., count, strand
    val(min_value)

    output:
    tuple val(meta), path("*.paraclu.txt"), emit: clusters
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('paraclu'), val("10"), topic: versions, emit: versions_paraclu

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Full cluster hierarchy: paraclu output is deliberately not simplified with paraclu-cut.
    """
    awk -F '\\t' -v OFS='\\t' '{print \$1, \$6, \$2, \$5}' $ctss \\
        | sort -k1,1 -k2,2 -k3,3n \\
        | paraclu $min_value - \\
        > ${prefix}.paraclu.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.paraclu.txt
    """
}
