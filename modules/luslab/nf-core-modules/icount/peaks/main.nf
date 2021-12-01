include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ICOUNT_PEAKS {
    tag "$meta.id"
    label "low_cores"
    label "low_mem"
    label "regular_queue"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::icount=2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/icount:2.0.0--py_1"
    } else {
        container "quay.io/biocontainers/icount:2.0.0--py_1"
    }

    input:
    tuple val(meta), path(bed)
    path(segmentation)

    output:
    tuple val(meta), path("*.peaks.bed.gz"), emit: peaks
    tuple val(meta), path("*.scores.tsv")  , emit: scores
    path "*.version.txt"                   , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    iCount peaks \\
        $segmentation \\
        $bed \\
        ${prefix}.peaks.bed.gz \\
        --scores ${prefix}.scores.tsv \\
        $options.args
    echo \$(iCount -v) > ${software}.version.txt
    """
}
