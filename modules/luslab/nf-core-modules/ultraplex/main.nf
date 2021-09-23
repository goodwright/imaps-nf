// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ULTRAPLEX {
    tag "${meta.id}"
    label "max_cores"
    label "max_memory"
    label "regular_queue"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::ultraplex=1.1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ultraplex:1.1.5--py36hc5360cc_0"
    } else {
        container "quay.io/biocontainers/ultraplex:1.1.5--py36hc5360cc_0"
    }

    input:
    tuple val(meta), path(reads)
    path(barcode_file)

    output:
    tuple val(meta), path("*[!no_match].fastq.gz")             , emit: fastq
    tuple val(meta), path("*no_match.fastq.gz"), optional: true, emit: no_match_fastq
    path "*.log"                                               , emit: report
    path "*.version.txt"                                       , emit: version

    script:
    def software = getSoftwareName(task.process)
    args = ""
    if(options.args && options.args != '') {
        ext_args = options.args
        args += ext_args.trim()
    }
    read_list = reads.collect{it.toString()}
    if (read_list.size > 1){
        ultraplex_command = """ultraplex \\
        --inputfastq ${read_list[0]} \\
        --input_2 ${read_list[1]} \\
        --barcodes $barcode_file \\
        --threads ${task.cpus} ${args}"""
    } else {
        ultraplex_command = """ultraplex \\
        --inputfastq ${read_list[0]} \\
        --barcodes $barcode_file \\
        --threads ${task.cpus} ${args}"""
    }
    if (params.verbose){
        println ("[MODULE] ultraplex command: " + ultraplex_command)
    }

    """
    ${ultraplex_command}
    echo \$(ultraplex --version 2>&1) | sed 's/^.*ultraplex //; s/Using.*\$//' > ${software}.version.txt
    """
}
