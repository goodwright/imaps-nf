process UMICOLLAPSE {
    tag "$meta.id"
    label "process_high"

    container 'elly1502/umicollapse:latest'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bam")             , emit: bam
    tuple val(meta), path("*.log")             , emit: log
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    ./umicollapse \\
        bam \\
        -i $bam \\
        -o ${prefix}.bam \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        umicollapse: NA
    END_VERSIONS
    """
}

