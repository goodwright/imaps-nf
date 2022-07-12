#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { DEMULTIPLEX } from '../subworkflows/demultiplex'
include { PRIMARY_CLIP_ANALYSIS } from '../subworkflows/primaryclipanalysis'
include { NCRNA_ANALYSIS } from '../subworkflows/ncrna'
include { CLIP_QUALITY_CHECK } from '../subworkflows/clipqualitycheck'

workflow {

    DEMULTIPLEX (
        params.annotation,
        params.multiplexed_fastq
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

    CLIP_QUALITY_CHECK(
        DEMULTIPLEX.out.fastqc_html,
        DEMULTIPLEX.out.fastqc_zip,
        PRIMARY_CLIP_ANALYSIS.out.trimgalore_log,
        PRIMARY_CLIP_ANALYSIS.out.bowtie_align_log,
        PRIMARY_CLIP_ANALYSIS.out.star_align_log_final,
        PRIMARY_CLIP_ANALYSIS.out.umicollapse_log,
        PRIMARY_CLIP_ANALYSIS.out.crosslinks,
        PRIMARY_CLIP_ANALYSIS.out.icount_peaks,
        PRIMARY_CLIP_ANALYSIS.out.paraclu_peaks,
        PRIMARY_CLIP_ANALYSIS.out.clippy_peaks
    )

}
