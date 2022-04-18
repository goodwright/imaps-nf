#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { DEMULTIPLEX } from '../subworkflows/demultiplex'
include { PRIMARY_ANALYSIS } from '../subworkflows/primaryanalysis'
include { NCRNA_ANALYSIS } from '../subworkflows/ncrna'

workflow {

    DEMULTIPLEX (
        params.annotation,
        params.multiplexed_fastq
    )

    // Get channel for the genome
    DEMULTIPLEX.out.map{pair -> params[pair[0].species + "_genome"]}.set { ch_genome }

    // Run Primary Analysis pipeline where necessary
    DEMULTIPLEX.out.filter{pair -> pair[0].pipeline == "Primary Analysis"}
    .set{ ch_primary_analysis_reads }
    PRIMARY_ANALYSIS (
        ch_primary_analysis_reads, // [meta, reads] pairs
        ch_genome
    )

    // Run ncRNA pipeline where necessary
    DEMULTIPLEX.out.filter{pair -> pair[0].pipeline == "Non-Coding RNA Analysis"}
    .set { ch_ncrna_reads }
    NCRNA_ANALYSIS (
        ch_ncrna_reads, // [meta, reads] pairs
        ch_genome
    )

}