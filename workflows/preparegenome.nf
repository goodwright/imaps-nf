#!/usr/bin/env nextflow

/* 
PREPARE GENOME WORKFLOW
Steps:
1) Index the genome using STAR module
2) Prepare genome segmentation file using iCount
*/

nextflow.enable.dsl=2


include { STAR_GENOMEGENERATE         } from '../modules/nf-core/modules/star/genomegenerate/main'
include { GUNZIP                      } from '../modules/nf-core/modules/gunzip/main'



workflow {

// Check if genome fasta is compressed, if so decompress it
    genome_fasta_file = file(params.fasta)

    if (params.fasta.matches(".*gz")) {
        GUNZIP(
            params.fasta
        )
        ch_fasta = GUNZIP.out.gunzip
    } else {
        Channel
            .fromPath(params.fasta, checkIfExists: true)
            .into { ch_fasta }
    }


    STAR_GENOMEGENERATE(
        ch_fasta,
        params.gtf
    )

}