#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { DEMULTIPLEX } from '../subworkflows/demultiplex'
include { PRIMARY_CLIP_ANALYSIS } from '../subworkflows/primaryclipanalysis'
include { NCRNA_ANALYSIS } from '../subworkflows/ncrna'

workflow {

    DEMULTIPLEX (
        params.annotation,
        params.multiplexed_fastq,
        params.fastqc_single_end
    )

    // Run Primary CLIP Analysis
    DEMULTIPLEX.out.fastq
    .filter{pair -> pair[0].pipeline == "Primary CLIP Analysis"}
    .map{pair -> [pair[0], pair[1], params[pair[0].species + "_genome"]]}
    .set{ ch_primary_clip_analysis_reads }

    PRIMARY_CLIP_ANALYSIS (
        ch_primary_clip_analysis_reads, // [meta, reads, genome_name] triplets
    )

    // Run Non-Coding RNA Analysis
    DEMULTIPLEX.out.fastq
    .filter{pair -> pair[0].pipeline == "Non-Coding RNA Analysis"}
    .map{pair -> [pair[0], pair[1], file(params[pair[0].species + "_genome"])]}
    .set{ ch_ncrna_reads }

    NCRNA_ANALYSIS (
        ch_ncrna_reads, // [meta, reads, genome_name] triplets
    )

}