// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '10'

process PARACLU_CONVERT {
    tag "$meta.id"
    label "low_cores"
    label "low_mem"
    label "regular_queue"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::paraclu=10" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/paraclu:10--h9a82719_1"
    } else {
        container "quay.io/biocontainers/paraclu:10--h9a82719_1"
    }

    input:
    tuple val(meta), path(peaks)

    output:
    tuple val(meta), path("*.peaks.bed.gz"),    emit: peaks
    path "*.version.txt",                       emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    gzip -d -c $peaks | \
        awk '{OFS = "\t"}{print \$1, \$3-1, \$4, ".", \$6, \$2}' |
        sort -k1,1 -k2,2n | \
        gzip > ${prefix}.peaks.bed.gz

    echo $VERSION > ${software}.version.txt
    """
}
