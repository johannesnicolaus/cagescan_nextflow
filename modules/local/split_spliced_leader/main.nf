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
    // Leaders: 'SEQ' or 'NAME=SEQ,NAME=SEQ,...'; unnamed ones are SL (one leader) or SL1, SL2, ...
    def entries = sl_sequence.tokenize(',')*.trim().findAll { entry -> entry }
    def leaders = entries.withIndex().collect { entry, i ->
        def parts = entry.tokenize('=')
        parts.size() > 1 ? [parts[0], parts[1].toUpperCase()] : [entries.size() == 1 ? 'SL' : "SL${i + 1}", parts[0].toUpperCase()]
    }
    def fasta = leaders.collect { name, seq -> ">SL.${name}\\n${seq}\\n>SLG.${name}\\nG${seq}" }.join('\\n')
    // A leader must start at READ1 base 1 or after one G (the extra G); mismatches are allowed but no indels, so that
    // it is cut exactly at its end. It is trimmed, and the pairs are written to separate files and renamed
    // SL.<leader>:<name> or SLG.<leader>:<name> (leader found) / no_adapter:<name>, so that the SL status and the
    // leader survive alignment. No {comment} in the new name: STAR rejects header lines that end in a space.
    """
    printf '${fasta}\\n' > leaders.fa
    cutadapt \\
        --cores $task.cpus \\
        -g ^file:leaders.fa \\
        --no-indels \\
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
