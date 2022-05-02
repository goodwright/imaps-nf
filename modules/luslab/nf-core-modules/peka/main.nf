process PEKA {
    tag "$meta.id"
    label "low_cores"
    label "low_mem"
    label "regular_queue"

    conda (params.enable_conda ? "bioconda::peka=0.1.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/peka:0.1.6--pyhdfd78af_0' :
        'quay.io/biocontainers/peka:0.1.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(peaks)
    tuple val(meta), path(crosslinks)
    path(genome)
    path(fai)
    path(gtf)

    output:
    tuple val(meta), path("*mer_cluster_distribution*"), emit: cluster,      optional: true
    tuple val(meta), path("*mer_distribution*"),         emit: distribution, optional: true
    tuple val(meta), path("*.pdf"),                      emit: pdf,          optional: true
    path "versions.yml",                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
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
        $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peka: 0.1.6
    END_VERSIONS
    """
}
