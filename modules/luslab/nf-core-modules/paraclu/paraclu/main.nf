process PARACLU_PARACLU {
    tag "$meta.id"
    label "low_cores"
    label "low_mem"
    label "regular_queue"

    conda (params.enable_conda ? "bioconda::paraclu=10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/paraclu:10--h9a82719_1' :
        'quay.io/biocontainers/paraclu:10--h9a82719_1' }"

    input:
    tuple val(meta), path(crosslinks)

    output:
    tuple val(meta), path("*.sigxls.tsv.gz"), emit: sigxls
    path "versions.yml",                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def paraclu_version = '10'
    """
    gzip -d -c $crosslinks | \
        awk '{OFS = "\t"}{print \$1, \$6, \$2+1, \$5}' | \
        sort -k1,1 -k2,2 -k3,3n > paraclu_input.tsv

    paraclu \
        $args \
        paraclu_input.tsv | \
        gzip > ${prefix}.sigxls.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
        paraclu: $paraclu_version
    END_VERSIONS
    """
}
