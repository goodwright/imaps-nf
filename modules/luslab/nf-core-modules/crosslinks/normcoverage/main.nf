process CROSSLINKS_NORMCOVERAGE {
    tag "$meta.id"
    label "low_cores"
    label "low_mem"
    label "regular_queue"

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(crosslinks)

    output:
    tuple val(meta), path("$prefix*.bedgraph.gz"), emit: bedgraph
    path "versions.yml",                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    TOTAL=`gunzip -c $crosslinks | awk 'BEGIN {total=0} {total=total+\$5} END {print total}'`

    gzip -d -c $crosslinks | \
        awk -v total=\$TOTAL '{printf "%s\\t%i\\t%i\\t%s\\t%f\\t%s\\n", \$1, \$2, \$3, \$4, 1000000*\$5/total, \$6}' | \
        awk '{OFS = "\t"}{if (\$6 == "+") {print \$1, \$2, \$3, \$5} else {print \$1, \$2, \$3, -\$5}}' | \
        sort -k1,1 -k2,2n | \
        gzip > ${prefix}.norm.bedgraph.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
