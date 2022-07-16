process FIND_LONGEST_TRANSCRIPT {
    tag "$gtf"
    label "process_medium"
    container "quay.io/biocontainers/peka:0.1.4--pyhdfd78af_0"
    
    input:
        path gtf

    output:
        path "*.longest_transcript.txt",          emit: longest_transcript
        path "*.transcriptome_index.fa.fai",      emit: transcriptome_index
    
    script:
        template "find_longest_transcript.py"
}