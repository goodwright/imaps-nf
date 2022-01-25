#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LONGEST_TRANSCRIPT {
    tag "$gtf"
    label "low_cores"
    label "low_mem"
    label "regular_queue"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
    
    input:
    path(gtf)

    output:
    "$long_tx",     emit: longest_transcript
    "$tr_ind",      emit: transcriptome_index

    shell:
    long_tx = 'longest_transcript.txt'
    tr_ind = 'transcriptome_index.fa.fai'

    //SHELL
    '''
    #Table of CDS lengths for each protein coding transcript.
    cat !{gtf} | \
    awk -F "\t" '$3 == "CDS" { print $0 }' |
    tr -d "\"" | \
    grep "transcript_type protein_coding" | \
    awk -F"\t" '{ OFS="\t" } { match($0, "gene_id [[:alnum:].]*", a) } { match($0, "transcript_id [[:alnum:].]*", b) } { print a[0], b[0], 1 + $5 - $4 }' |
    tr -d "gene_id |transcript_id " |
    awk -F "\t" 'BEGIN { SUBSEP = OFS = FS } { s[$1, $2] += $3 } END { for (i in s) { print i, s[i] } }' |
    sort -k1,1 -k3,3nr -k2,2 > cds_lengths.txt

    #Table of transcript lengths for each protein coding transcript.
    cat !{gtf} | \
    awk -F "\t" '$3 == "exon" { print $0 }' |
    tr -d "\"" | \
    grep "transcript_type protein_coding" | \
    awk -F"\t" '{ OFS="\t" } { match($0, "gene_id [[:alnum:].]*", a) } { match($0, "transcript_id [[:alnum:].]*", b) } { print a[0], b[0], 1 + $5 - $4 }' |
    tr -d "gene_id |transcript_id " |
    awk -F "\t" 'BEGIN { SUBSEP = OFS = FS } { s[$1, $2] += $3 } END { for (i in s) { print i, s[i] } }' |
    sort -k1,1 -k3,3nr -k2,2 > tx_lengths.txt

    #Use tables of CDS and transcript lengths to select longest protein coding transcript per gene.
    join -j 2 -o 1.1,1.2,1.3,2.3 <(sort -t$'\t' cds_lengths.txt) <(sort -t$'\t' tx_lengths.txt) | \
    sort -k1,1 -k3,3nr -k4,4nr | awk '$1 != x { print } { x = $1 }' | awk '{ print $2 }' > !{long_tx}

    #Make fake fasta index based on transcript lengths.
    awk -F"\t" '{ OFS="\t" } { print $2, $3, sum, 80, 81; sum += $3 + 1 }' tx_lengths.txt > !{tr_ind}
    '''
}