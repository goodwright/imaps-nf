process CLIPPY {
    tag "$meta.id"
    label "avg_cores"
    label "high_mem"
    label "regular_queue"

    // A dependency issue in Dash (https://github.com/plotly/dash/issues/1992)
    // means we have to specify a version for werkzeug. This should be fixed in
    // the bioconda recipe in the future.
    conda (params.enable_conda ? "bioconda::clippy=1.3.3 conda-forge::werkzeug=2.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clippy:1.3.3--pyhdfd78af_0' :
        'quay.io/biocontainers/clippy:1.3.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(crosslinks)
    path(gtf)
    path(fai)

    output:
    tuple val(meta), path("*_Peaks.bed.gz"),   emit: peaks
    tuple val(meta), path("*_Summits.bed.gz"), emit: summits
    path "versions.yml",                       emit: versions

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
