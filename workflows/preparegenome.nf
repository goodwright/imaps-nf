#!/usr/bin/env nextflow

/* 
PREPARE GENOME WORKFLOW
Steps:
1) Index the genome using STAR module
2) Prepare genome segmentation file using iCount
3) Identify coding isoforms with longest CDS per gene.
*/

nextflow.enable.dsl=2

include { GUNZIP as DNA_GUNZIP        } from '../modules/nf-core/modules/gunzip/main'
include { GUNZIP as RNA_GUNZIP        } from '../modules/nf-core/modules/gunzip/main'
include { STAR_GENOMEGENERATE         } from '../modules/nf-core/modules/star/genomegenerate/main'
include { BOWTIE_BUILD                } from '../modules/nf-core/modules/bowtie/build/main'
include { SAMTOOLS_FAIDX              } from '../modules/nf-core/modules/samtools/faidx/main'
include { ICOUNT_SEGMENT              } from '../modules/luslab/nf-core-modules/icount/segment/main'
include { LONGEST_TRANSCRIPT          } from '../modules/local/find_longest_transcript/main'
include { FILTER_GTF                  } from '../modules/local/filter_gtf/main'
include { RESOLVE_UNANNOTATED         } from '../modules/local/resolve_unannotated/main'

workflow {

    // If genome is compressed, uncompress it
    if (params.fasta.matches(".*gz")) {
        ch_fasta = DNA_GUNZIP ( file(params.fasta) ).gunzip
    } else {
        ch_fasta = file(params.fasta)
    }

    // If smrna genome is compressed, uncompress it
    if (params.smrna_fasta.matches(".*gz")) {
        ch_smrna_fasta = RNA_GUNZIP ( file(params.smrna_fasta) ).gunzip
    } else {
        ch_smrna_fasta = file(params.smrna_fasta)
    }


    FILTER_GTF ( params.gtf )

    // Index the genome using STAR module
    STAR_GENOMEGENERATE (
        ch_fasta,
        params.gtf
    )

    BOWTIE_BUILD (
        ch_smrna_fasta
    )

    // Create a FAI genome index using samtools
    ch_fai = SAMTOOLS_FAIDX ( ch_fasta ).fai


    // Create genome segmentation file using iCount
    ICOUNT_SEGMENT (
        FILTER_GTF.out.filtered_gtf,
        ch_fai
    )


    RESOLVE_UNANNOTATED(
        ICOUNT_SEGMENT.out.regions, // segment
        FILTER_GTF.out.filtered_gtf, // gtf
        ch_fai, // fai
    )

    // Find the longest CDS transcript per gene.
    LONGEST_TRANSCRIPT (
        params.gtf
    )

}