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

    // Run Primary Analysis
    DEMULTIPLEX.out
    .filter{pair -> pair[0].pipeline == "Primary Analysis"}
    .set{ ch_primary_analysis_reads }

    PRIMARY_ANALYSIS (
        ch_primary_analysis_reads, // [meta, reads] pairs
        ch_primary_analysis_reads.map{pair -> params[pair[0].species + "_genome"]}
    )

    // Run Non-Coding RNA Analysis
    DEMULTIPLEX.out
    .filter{pair -> pair[0].pipeline == "Non-Coding RNA Analysis"}
    .set{ ch_ncrna_reads }

    NCRNA_ANALYSIS (
        ch_ncrna_reads, // [meta, reads] pairs
        ch_ncrna_reads.map{pair -> params[pair[0].species + "_genome"]}
    )

}