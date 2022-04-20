include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PEKA {
    tag "$meta.id"
    label "low_cores"
    label "low_mem"
    label "regular_queue"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::peka=0.1.6" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/peka:0.1.6--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/peka:0.1.6--pyhdfd78af_0"
    }

    input:
    tuple val(meta), path(peaks)
    tuple val(meta), path(crosslinks)
    path(genome)
    path(fai)
    path(gtf)

    output:
    tuple val(meta), path("*mer_cluster_distribution*"), emit: cluster
    tuple val(meta), path("*mer_distribution*"),         emit: distribution
    tuple val(meta), path("*.pdf"),                      emit: pdf
    path "*.version.txt",                                emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    # If the modification date and time of the fai is before the fasta then
    # there will be an error. Touching the file first avoids that.
    touch $fai
    mkdir tmp
    TMPDIR=\$(pwd)/tmp peka \
        -i $peaks \
        -x $crosslinks \
        -g $genome \
        -gi $fai \
        -r $gtf \
        $options.args
    echo "0.1.6" > ${software}.version.txt
    """
}
