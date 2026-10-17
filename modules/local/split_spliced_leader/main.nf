process SPLIT_SPLICED_LEADER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:5.2--py311haab0aaa_0' :
        'quay.io/biocontainers/cutadapt:5.2--py311haab0aaa_0' }"

    input:
    tuple val(meta), path(reads)
    val sl_sequence

    output:
    tuple val(meta), path("*.SL_{1,2}.fastq.gz")  , emit: sl
    tuple val(meta), path("*.noSL_{1,2}.fastq.gz"), emit: nosl
    tuple val(meta), path("*.sl_split.log")       , emit: log
    tuple val("${task.process}"), val('cutadapt'), eval('cutadapt --version'), topic: versions, emit: versions_cutadapt

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The SL must start at READ1 base 1 or 2 (room for the extra G). It is trimmed, and the pairs are written to
    // separate files and renamed SL:<name> / no_adapter:<name>, so that the SL status survives alignment.
    // No {comment} in the new name: STAR rejects header lines that end in a space.
    """
    cutadapt \\
        --cores $task.cpus \\
        -g 'SL=^N${sl_sequence}' \\
        --rename '{r1.adapter_name}:{id}' \\
        $args \\
        -o ${prefix}.SL_1.fastq.gz \\
        -p ${prefix}.SL_2.fastq.gz \\
        --untrimmed-output ${prefix}.noSL_1.fastq.gz \\
        --untrimmed-paired-output ${prefix}.noSL_2.fastq.gz \\
        $reads \\
        > ${prefix}.sl_split.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for f in SL_1 SL_2 noSL_1 noSL_2; do echo '' | gzip > ${prefix}.\$f.fastq.gz; done
    touch ${prefix}.sl_split.log
    """
}
