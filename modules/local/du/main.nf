process DU {
    label "min_cores"
    label "min_mem"
    label "regular_queue"

    tag "$meta.id"

    container "biocontainers/biocontainers:v1.2.0_cv1"

    input:
    tuple val(meta), path(input_file), path(bai)

    output:
    tuple val(meta), stdout, emit: size

    script:
    """
    echo -n "\$(du -kL $input_file | awk '{print(\$1)}')"
    """
}
