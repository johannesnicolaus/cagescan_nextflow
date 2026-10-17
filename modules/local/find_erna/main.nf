process FIND_ERNA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0' :
        'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0' }"

    input:
    tuple val(meta), path(clusters), path(count_files, stageAs: 'counts/*'), path(transcripts)
    path(annotation)  // reference GTF, or [] (no distal filter)

    output:
    tuple val(meta), path("*.erna_candidates.bed"), emit: bed
    tuple val(meta), path("*.erna_candidates.tsv"), emit: tsv
    tuple val(meta), path("*.erna_stats_mqc.tsv") , emit: stats
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def files      = [count_files].flatten()
    def counts     = files.find { f -> f.name.endsWith('.tss_cluster_counts.tsv') }
    def sl_counts  = files.find { f -> f.name.endsWith('.tss_cluster_sl_counts.tsv') }
    if (!counts) {
        error "FIND_ERNA: no *.tss_cluster_counts.tsv among ${files*.name}"
    }
    def sl_arg     = sl_counts ? "--sl-counts ${sl_counts}" : ''
    def annot_arg  = annotation ? "--annotation ${annotation}" : ''
    """
    find_erna.py \\
        --prefix $prefix \\
        --clusters $clusters \\
        --counts $counts \\
        $sl_arg \\
        --transcripts $transcripts \\
        $annot_arg \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.erna_candidates.bed ${prefix}.erna_candidates.tsv ${prefix}.erna_stats_mqc.tsv
    """
}
