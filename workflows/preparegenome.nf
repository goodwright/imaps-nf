#!/usr/bin/env nextflow

/* 
PREPARE GENOME WORKFLOW
Steps:
1) Index the genome using STAR module
2) Prepare genome segmentation file using iCount
3) Identify coding isoforms with longest CDS per gene.
*/

nextflow.enable.dsl=2

/*
include { GUNZIP                      } from '../modules/nf-core/modules/gunzip/main'
include { STAR_GENOMEGENERATE         } from '../modules/nf-core/modules/star/genomegenerate/main'
include { SAMTOOLS_FAIDX              } from '../modules/nf-core/modules/samtools/faidx/main'
include { ICOUNT_SEGMENT              } from '../modules/luslab/nf-core-modules/icount/segment/main'
*/
include { LONGEST_TRANSCRIPT          } from '../modules/local/find_longest_transcript/main'

workflow {

    /*
    // If genome is compressed, uncompress it
    if (params.fasta.matches(".*gz")) {
        ch_fasta = GUNZIP ( params.fasta ).gunzip
    } else {
        ch_fasta = file(params.fasta)
    }

    // Index the genome using STAR module
    STAR_GENOMEGENERATE (
        ch_fasta,
        params.gtf
    )

    // Create a FAI genome index using samtools
    ch_fai = SAMTOOLS_FAIDX ( ch_fasta ).fai


    // Create genome segmentation file using iCount
    ICOUNT_SEGMENT (
        params.gtf,
        ch_fai
    )
    */

    // Find the longest CDS transcript per gene.
    LONGEST_TRANSCRIPT (
        params.gtf
    )

}