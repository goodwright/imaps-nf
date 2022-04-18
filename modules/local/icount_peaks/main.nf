include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ICOUNT_PEAKS {
    tag "$meta.id"
    label "low_cores"
    label "low_mem"
    label "regular_queue"

    container "quay.io/biocontainers/icount-mini:2.0.3--pyh5e36f6f_0"

    input:
    tuple val(meta), path(bed), path(sigxls)

    output:
    tuple val(meta), path("*peaks.bed.gz"), emit: peaks
    path "*.version.txt"                   , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    iCount-Mini peaks \\
        $bed \\
        $sigxls \\
        ${prefix}.peaks.bed.gz \\
        $options.args
    echo \$(iCount-Mini -v) > ${software}.version.txt
    """
}