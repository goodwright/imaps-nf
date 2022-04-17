process FILTER_GTF {
    tag "$gtf"
    label "process_medium"
    container "quay.io/biocontainers/peka:0.1.4--pyhdfd78af_0"
    
    input:
    path(gtf)

    output:
    path "*.gtf", emit: post_filtering_gtf
    
    script:
    template 'FilterGtf.py'
}
