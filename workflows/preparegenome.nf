#!/usr/bin/env nextflow

/* 
PREPARE GENOME WORKFLOW
Steps:
1) Index the genome using STAR module
2) Prepare genome segmentation file using iCount
*/

nextflow.enable.dsl=2

include { GUNZIP                      } from '../modules/nf-core/modules/gunzip/main'
include { STAR_GENOMEGENERATE         } from '../modules/nf-core/modules/star/genomegenerate/main'
include { SAMTOOLS_FAIDX              } from '../modules/nf-core/modules/samtools/faidx/main'
include { ICOUNT_SEGMENT              } from '../modules/luslab/nf-core-modules/icount/segment/main'



workflow {

// Check if genome fasta is compressed, if so decompress it
    genome_fasta_file = file(params.fasta)

    if (params.fasta.matches(".*gz")) {
        ch_fasta = GUNZIP( params.fasta ).gunzip
    } else {
        ch_fasta = file(params.fasta)
    }

// STEP 1 - Index the genome using STAR module
    STAR_GENOMEGENERATE(
        ch_fasta,
        params.gtf
    )

//STEP 2 - Prepare genome segmentation file using iCount
//First need to create an index for our genome fasta file
    ch_fai = SAMTOOLS_FAIDX( ch_fasta ).fai

    ICOUNT_SEGMENT(
        params.gtf,
        ch_fai
    )

}