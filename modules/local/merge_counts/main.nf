process MERGE_COUNTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0' :
        'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0' }"

    input:
    tuple val(meta), path(files, stageAs: 'input/*')
    path(clusters)  // tss_clusters.tsv; required for mode 'tss', [] otherwise
    val(mode)       // 'tss' or 'transcripts'

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!(mode in ['tss', 'transcripts'])) {
        error "MERGE_COUNTS: mode must be 'tss' or 'transcripts', got '${mode}'"
    }
    def clusters_arg = mode == 'tss' ? "--clusters $clusters" : ''
    """
    merge_counts.py $mode \\
        --prefix $prefix \\
        $clusters_arg \\
        $files
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${mode}_counts.tsv
    """
}
