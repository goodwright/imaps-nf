process RESOLVE_UNANNOTATED {
    tag "Resolving $filtered_segmentation and $unfiltered_segmentation"
    label "process_medium"
    container "quay.io/biocontainers/peka:0.1.4--pyhdfd78af_0"

    input:
    path(filtered_segmentation)
    path(unfiltered_segmentation), stageAs: "filtered.regions.gtf.gz"
    path(gtf)
    path(fai)

    output:
    path "*.gtf", emit: annotated_gtf
    
    script:
    template 'ResolveUnannotated.py'
}
