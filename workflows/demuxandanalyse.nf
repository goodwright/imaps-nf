#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { DEMULTIPLEX } from '../subworkflows/demultiplex'
include { PRIMARY_CLIP_ANALYSIS } from '../subworkflows/primaryclipanalysis'
include { NCRNA_ANALYSIS } from '../subworkflows/ncrna'

workflow {

    DEMULTIPLEX (
        params.annotation,
        params.multiplexed_fastq
    )

    // Run Primary CLIP Analysis
    DEMULTIPLEX.out
    .filter{pair -> pair[0].pipeline == "Primary CLIP Analysis"}
    .set{ ch_primary_clip_analysis_reads }

    PRIMARY_CLIP_ANALYSIS (
        ch_primary_clip_analysis_reads, // [meta, reads] pairs
        ch_primary_clip_analysis_reads.map{pair -> params[pair[0].species + "_genome"]}
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