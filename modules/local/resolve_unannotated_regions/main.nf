process RESOLVE_UNANNOTATED_REGIONS {
    tag "$filtered_regions/$unfiltered_regions"
    label "process_medium"
    container "quay.io/biocontainers/peka:0.1.4--pyhdfd78af_0"

    input:
    path(filtered_regions)
    path(unfiltered_regions), stageAs: "filtered.regions.gtf.gz"
    path(gtf)
    path(fai)

    output:
    path "*.gtf", emit: annotated_gtf
    
    script:
    template 'ResolveUnannotated.py'
}