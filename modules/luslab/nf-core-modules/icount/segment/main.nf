process ICOUNT_SEGMENT {
    tag "$gtf"
    label "low_cores"
    label "low_mem"
    label "regular_queue"

    conda (params.enable_conda ? "bioconda::icount-mini=2.0.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/icount-mini:2.0.3--pyh5e36f6f_0' :
        'quay.io/biocontainers/icount-mini:2.0.3--pyh5e36f6f_0' }"

    input:
    path(gtf)
    path(fai)

    output:
    path("${prefix}.gtf"), emit: gtf
    path("regions.gtf.gz"),          emit: regions
    path "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def filename = "icount_segmentation"
    prefix       = task.ext.prefix ? "${filename}${options.suffix}" : "${filename}"
    """
    iCount-Mini segment \\
        $gtf \\
        ${prefix}.gtf \\
        $fai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iCount-Mini: \$(iCount-Mini -v)
    END_VERSIONS
    """
}
