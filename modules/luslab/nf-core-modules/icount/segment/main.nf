include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ICOUNT_SEGMENT {
    tag "$gtf"
    label "low_cores"
    label "low_mem"
    label "regular_queue"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::icount=2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/icount:2.0.0--py_1"
    } else {
        container "quay.io/biocontainers/icount:2.0.0--py_1"
    }

    input:
    path(gtf)
    path(fai)

    output:
    path("*.gtf")       , emit: gtf
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def filename = "icount_segmentation"
    def prefix   = options.suffix ? "${filename}${options.suffix}" : "${filename}"
    """
    iCount segment \\
        $gtf \\
        ${prefix}.gtf \\
        $fai
    echo \$(iCount -v) > ${software}.version.txt
    """
}
