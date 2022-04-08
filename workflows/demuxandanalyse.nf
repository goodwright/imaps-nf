#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { DEMULTIPLEX } from '../subworkflows/demultiplex'
include { PRIMARY_ANALYSIS } from '../subworkflows/primaryanalysis'
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main.nf'

workflow {

    DEMULTIPLEX (
        params.csv,
        params.multiplexed_fastq
    )

    DEMULTIPLEX.out.ch_reads_with_meta
        .map{pair -> params[pair[0].species + "_genome"]}
        .set { ch_genome }

    PRIMARY_ANALYSIS (
        DEMULTIPLEX.out.ch_reads_with_meta, // [meta, reads] pairs
        ch_genome
    )

    Channel
        .concat(
            DEMULTIPLEX.out.fastqc_html,
            DEMULTIPLEX.out.fastqc_zip,
            PRIMARY_ANALYSIS.out.trimgalore_log,
            PRIMARY_ANALYSIS.out.bowtie_align_log,
            PRIMARY_ANALYSIS.out.star_align_log_final,
            PRIMARY_ANALYSIS.out.clip_qc_log
        )
        .view()

    // MULTIQC ()

}
