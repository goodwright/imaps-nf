process PARACLU_CUT {
    tag "$meta.id"
    label "low_cores"
    label "low_mem"
    label "regular_queue"

    conda (params.enable_conda ? "bioconda::paraclu=10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/paraclu:10--h9a82719_1' :
        'quay.io/biocontainers/paraclu:10--h9a82719_1' }"

    input:
    tuple val(meta), path(sigxls)

    output:
    tuple val(meta), path("*.peaks.tsv.gz"), emit: peaks
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def paraclu_version = '10'
    """
    gzip -d -c $sigxls | \
        paraclu-cut \
        $args | \
        gzip > ${prefix}.peaks.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
        paraclu: $paraclu_version
    END_VERSIONS
    """
}
