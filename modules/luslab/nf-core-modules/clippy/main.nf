include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CLIPPY {
    tag "$meta.id"
    label "avg_cores"
    label "high_mem"
    label "regular_queue"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::clippy=1.3.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/clippy:1.3.1--pyhdfd78af_2"
    } else {
        container "quay.io/biocontainers/clippy:1.3.1--pyhdfd78af_2"
    }

    input:
    tuple val(meta), path(crosslinks)
    path(gtf)
    path(fai)

    output:
    tuple val(meta), path("*_broadPeaks.bed.gz"), emit: peaks
    tuple val(meta), path("*[0-9].bed.gz"),       emit: summits
    path "*.version.txt",                         emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    clippy -i $crosslinks \
        -o $prefix \
        -a $gtf \
        -g $fai \
        -t ${task.cpus} \
        $options.args
    gzip -n *_broadPeaks.bed
    gzip -n *[0-9].bed
    echo \$(clippy -v) > ${software}.version.txt
    """
}
