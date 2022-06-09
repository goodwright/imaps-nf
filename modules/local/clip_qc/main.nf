process CLIP_QC {
    tag "Performing QC"
    label "min_cores"
    label "min_mem"
    label "regular_queue"

    conda (params.enable_conda ? "bioconda::pybedtools==0.9.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pybedtools:0.9.0--py39h9a82719_0"
    } else {
        container "quay.io/biocontainers/pybedtools:0.9.0--py39h9a82719_0"
    }

    input:
    path("premap/*")
    path("mapped/*")
    path("xlinks/*")
    path("icount/*")
    path("paraclu/*")
    path("clippy/*")

    output:
    path "*.tsv", emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'clip_qc.py'
}
