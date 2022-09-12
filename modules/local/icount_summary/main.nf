process ICOUNT_SUMMARY {
    tag "$meta.id"
    label "low_cores"
    label "low_mem"
    label "regular_queue"

    conda (params.enable_conda ? "bioconda::icount-mini=2.0.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/icount-mini:2.0.3--pyh5e36f6f_0"
    } else {
        container "quay.io/biocontainers/icount-mini:2.0.3--pyh5e36f6f_0"
    }

    input:
    tuple val(meta), path(bed)
    path(regions)

    output:
    tuple val(meta), path("*summary_type.tsv"), emit: summary_type
    tuple val(meta), path("*summary_subtype.tsv")  , emit: summary_subtype
    tuple val(meta), path("*summary_gene.tsv")  , emit: summary_gene
    path "versions.yml"                   , emit: versions

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    iCount-Mini summary \\
        $regions \\
        $bed \\
        . \\
        $args
    mv summary_type.tsv ${prefix}_summary_type.tsv
    mv summary_subtype.tsv ${prefix}_summary_subtype.tsv
    mv summary_gene.tsv ${prefix}_summary_gene.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iCount-Mini: \$(iCount-Mini -v)
    END_VERSIONS
    """
}