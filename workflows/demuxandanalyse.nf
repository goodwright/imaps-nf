#!/usr/bin/env nextflow

/* 
DEMULTIPLEX & ANALYSE
*/

nextflow.enable.dsl=2

include { ULTRAPLEX } from '../modules/luslab/nf-core-modules/ultraplex/main'    addParams( options: [:] )
include { INPUT_CHECK } from '../subworkflows/input_check'    addParams( options: [:] )
include { FASTQC } from '../modules/nf-core/modules/fastqc/main' addParams( options: [:] )
include { TRIMGALORE } from '../modules/nf-core/modules/trimgalore/main'  addParams( options: [:] )
include { BOWTIE_ALIGN } from '../modules/nf-core/modules/bowtie/align/main'    addParams( options: [:] )
include { STAR_ALIGN } from '../modules/nf-core/modules/star/align/main'    addParams( options: [:] )


workflow {
    
// Initialise channels
    ch_software_versions = Channel.empty()
//demultiplexing

    ch_input_meta = Channel.fromPath(params.input)
    ch_input_fasta = file(params.demultiplexed_fastq)

    

    INPUT_CHECK (
        ch_input_meta,
        ch_input_fasta
    )
//fastqc
 //   FASTQC (
  //      ch_input_fasta
  //  )
//trim-galore
//bowtie to small RNA
//unmapped to STAR GENOME
//UMI-TOOLS
//GET CROSSLINKS
//ICOUNT SUMMARY
//ICOUNT RNAMAPS
//GET COVERAGE
//PEKA
//PARACLU
//ICOUNT PEAKS
//CLIPPY
}