#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

// Process definition
process LONGEST_TRANSCRIPT {
    tag "$gtf"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
    
    // Get GNU awk
    conda (params.enable_conda ? "conda-forge::gawk=5.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gawk:5.1.0"
    } else {
        container "quay.io/biocontainers/gawk:5.1.0"
    }

    input:
        path gtf

    output:
        path "longest_transcript.txt",          emit: longest_transcript
        path "transcriptome_index.fa.fai",      emit: transcriptome_index

    script:
        //SHELL
        """
        #Table of CDS lengths for each protein coding transcript.
        cat ${gtf} | \
        gawk -F "\\t" '\$3 == "CDS" { print \$0 }' |
        tr -d "\\"" | \
        grep "transcript_type protein_coding" | \
        gawk -F"\\t" '{ OFS="\\t" } { match(\$0, "gene_id [[:alnum:].]*", a) } { match(\$0, "transcript_id [[:alnum:].]*", b) } { print a[0], b[0], 1 + \$5 - \$4 }' |
        tr -d "gene_id |transcript_id " |
        gawk -F "\\t" 'BEGIN { SUBSEP = OFS = FS } { s[\$1, \$2] += \$3 } END { for (i in s) { print i, s[i] } }' |
        sort -k1,1 -k3,3nr -k2,2 > cds_lengths.txt

        #Table of transcript lengths for each protein coding transcript.
        cat ${gtf} | \
        gawk -F "\\t" '\$3 == "exon" { print \$0 }' |
        tr -d "\\"" | \
        grep "transcript_type protein_coding" | \
        gawk -F"\\t" '{ OFS="\\t" } { match(\$0, "gene_id [[:alnum:].]*", a) } { match(\$0, "transcript_id [[:alnum:].]*", b) } { print a[0], b[0], 1 + \$5 - \$4 }' |
        tr -d "gene_id |transcript_id " |
        gawk -F "\\t" 'BEGIN { SUBSEP = OFS = FS } { s[\$1, \$2] += \$3 } END { for (i in s) { print i, s[i] } }' |
        sort -k1,1 -k3,3nr -k2,2 > tx_lengths.txt

        #Use tables of CDS and transcript lengths to select longest protein coding transcript per gene.
        gawk -F"\\t" 'NR==FNR { e[\$1\$2]=\$3;next} {print \$0 FS e[\$1\$2]}' tx_lengths.txt cds_lengths.txt | \
        sort -k1,1 -k3,3nr -k4,4nr | gawk '\$1 != x { print } { x = \$1 }' | gawk '{ print \$2 }' > longest_transcript.txt

        #Make fake fasta index based on transcript lengths.
        gawk -F"\\t" '{ OFS="\\t" } { print \$2, \$3, sum, 80, 81; sum += \$3 + 1 }' tx_lengths.txt > transcriptome_index.fa.fai
        """
}