process RENAME {
    tag "Renaming $input_file to $new_name"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    path(input_file)
    val(new_name)

    output:
    path("$new_name"), emit: renamed

    script:
    """
    cp $input_file $new_name
    """
}
