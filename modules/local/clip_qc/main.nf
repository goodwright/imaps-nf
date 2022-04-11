include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CLIP_QC {
    tag "Performing QC"
    label "process_low"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

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
    
    script:
    template 'clip_qc.py'
}
