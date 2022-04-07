include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RESOLVE_UNANNOTATED {
    tag "Resolving $filtered_segmentation and $unfiltered_segmentation"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::peka=0.1.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/peka:0.1.4--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/peka:0.1.4--pyhdfd78af_0"
    }

    input:
    path(filtered_segmentation)
    path(unfiltered_segmentation), stageAs: "filtered.regions.gtf.gz"
    path(gtf)
    path(fai)

    output:
    path "*.gtf", emit: annotated_gtf
    
    script:
    def args = task.ext.args ?: ''
    template 'ResolveUnannotated.py'
}
