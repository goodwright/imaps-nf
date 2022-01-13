#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { DEMULTIPLEX } from '../subworkflows/demultiplex'
include { PRIMARY_ANALYSIS } from '../subworkflows/primaryanalysis'

workflow {

    DEMULTIPLEX (
        params.csv,
        params.multiplexed_fastq
    )

    PRIMARY_ANALYSIS (
        DEMULTIPLEX.out, // [meta, reads] pairs
        params.smrna_genome,
        params.star_index,
        params.gtf,
        params.genome_fai,
        params.icount_regions,
        params.icount_segment
    )

}