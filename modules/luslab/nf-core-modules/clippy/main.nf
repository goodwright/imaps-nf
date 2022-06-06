process CLIPPY {
    tag "$meta.id"
    label "avg_cores"
    label "avg_mem"
    label "regular_queue"

    conda (params.enable_conda ? "bioconda::clippy=1.4.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clippy:1.4.1--pyhdfd78af_1' :
        'quay.io/biocontainers/clippy:1.4.1--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(crosslinks)
    path(gtf)
    path(fai)

    output:
    tuple val(meta), path("*_Peaks.bed.gz"),           emit: peaks
    tuple val(meta), path("*_Summits.bed.gz"),         emit: summits
    tuple val(meta), path("*_intergenic_regions.gtf"), emit: intergenic_gtf, optional: true
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    clippy -i $crosslinks \
        -o $prefix \
        -a $gtf \
        -g $fai \
        -t ${task.cpus} \
        $args
    gzip -n *_Peaks.bed
    gzip -n *_Summits.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clippy: \$(clippy -v)
    END_VERSIONS
    """
}
