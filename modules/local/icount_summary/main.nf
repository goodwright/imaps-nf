include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ICOUNT_SUMMARY {
    tag "$meta.id"
    label "low_cores"
    label "low_mem"
    label "regular_queue"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"icount_summary", meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::icount-mini=2.0.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/icount-mini:2.0.3--pyh5e36f6f_0"
    } else {
        container "quay.io/biocontainers/icount-mini:2.0.3--pyh5e36f6f_0"
    }

    input:
    tuple val(meta), path(bed)
    path(segmentation)

    output:
    tuple val(meta), path("*summary_type.tsv"), emit: summary_type
    tuple val(meta), path("*summary_subtype.tsv")  , emit: summary_subtype
    tuple val(meta), path("*summary_gene.tsv")  , emit: summary_gene
    path "*.version.txt"                   , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    iCount-Mini summary \\
        $segmentation \\
        $bed \\
        . \\
        $options.args
    mv summary_type.tsv ${prefix}_summary_type.tsv
    mv summary_subtype.tsv ${prefix}_summary_subtype.tsv
    mv summary_gene.tsv ${prefix}_summary_gene.tsv
    echo \$(iCount-Mini -v) > ${software}.version.txt
    """
}