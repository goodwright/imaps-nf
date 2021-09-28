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
include { GET_SOFTWARE_VERSIONS       } from '../modules/local/get_software_versions/main'


workflow {
// Initialise channels
    ch_software_versions = Channel.empty()

// Check if genome fasta is compressed, if so decompress it
    genome_fasta_file = file(params.fasta)

    if (params.fasta.matches(".*gz")) {
        ch_fasta = GUNZIP ( params.fasta ).gunzip
    } else {
        ch_fasta = file(params.fasta)
    }
// Collect software versions after a module is used. Then at the end of the pipleine they are gathered together
    ch_software_versions = ch_software_versions.mix(GUNZIP.out.version.ifEmpty(null))

// STEP 1 - Index the genome using STAR module
    STAR_GENOMEGENERATE (
        ch_fasta,
        params.gtf
    )
    ch_software_versions = ch_software_versions.mix(STAR_GENOMEGENERATE.out.version.ifEmpty(null))

//STEP 2 - Prepare genome segmentation file using iCount
//First need to create an index for our genome fasta file
    ch_fai = SAMTOOLS_FAIDX ( ch_fasta ).fai
    ch_software_versions = ch_software_versions.mix(SAMTOOLS_FAIDX.out.version.ifEmpty(null))

    ICOUNT_SEGMENT (
        params.gtf,
        ch_fai
    )
    ch_software_versions = ch_software_versions.mix(ICOUNT_SEGMENT.out.version.ifEmpty(null))

// Collect together all software versions and output them to a file
    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

}