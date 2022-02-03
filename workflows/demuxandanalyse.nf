#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { DEMULTIPLEX } from '../subworkflows/demultiplex'
include { PRIMARY_ANALYSIS } from '../subworkflows/primaryanalysis'

workflow {

    DEMULTIPLEX (
        params.csv,
        params.multiplexed_fastq
    )

    DEMULTIPLEX.out.map{pair -> params[pair[0].species + "_genome"]}.set { ch_genome }

    PRIMARY_ANALYSIS (
        DEMULTIPLEX.out, // [meta, reads] pairs
        ch_genome
    )

}