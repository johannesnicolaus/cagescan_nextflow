process IGV_SESSION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0' :
        'quay.io/biocontainers/pysam:0.24.1--py312hf5ad864_0' }"

    input:
    tuple val(meta), path(fasta), path(fai)
    path(gtf)        // reference annotation, or []
    val(samples)     // comma-separated sample names
    val(sets)        // comma-separated TSS cluster sets

    output:
    path("igv_session.xml"), emit: session
    path("igv/*")          , emit: files
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def gtf_arg = gtf ? "--gtf ${gtf}" : ''
    """
    igv_session.py \\
        --fasta $fasta \\
        --fai $fai \\
        $gtf_arg \\
        --samples $samples \\
        --sets $sets \\
        $args
    """

    stub:
    """
    mkdir igv
    touch igv_session.xml igv/genome.fa igv/genome.fa.fai
    """
}
