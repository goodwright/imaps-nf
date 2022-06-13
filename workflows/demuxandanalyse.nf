#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { DEMULTIPLEX } from '../subworkflows/demultiplex'
include { PRIMARY_CLIP_ANALYSIS } from '../subworkflows/primaryclipanalysis'
include { NCRNA_ANALYSIS } from '../subworkflows/ncrna'
include { MULTIQC } from '../modules/local/multiqc/main.nf'

workflow {

    DEMULTIPLEX (
        params.annotation,
        params.multiplexed_fastq
    )

    // Run Primary CLIP Analysis
    DEMULTIPLEX.out.ch_reads_with_meta
    .filter{pair -> pair[0].pipeline == "Primary CLIP Analysis"}
    .map{pair -> [pair[0], pair[1], params[pair[0].species + "_genome"]]}
    .set{ ch_primary_clip_analysis_reads }

    PRIMARY_CLIP_ANALYSIS (
        ch_primary_clip_analysis_reads, // [meta, reads, genome_name] triplets
    )

    // Run Non-Coding RNA Analysis
    DEMULTIPLEX.out.ch_reads_with_meta
    .filter{pair -> pair[0].pipeline == "Non-Coding RNA Analysis"}
    .map{pair -> [pair[0], pair[1], file(params[pair[0].species + "_genome"])]}
    .set{ ch_ncrna_reads }

    NCRNA_ANALYSIS (
        ch_ncrna_reads, // [meta, reads, genome_name] triplets
    )

    MULTIQC (
        Channel
            .empty()
            .concat(
                DEMULTIPLEX.out.fastqc_html.map{ vec -> vec[1] },
                DEMULTIPLEX.out.fastqc_zip.map{ vec -> vec[1] },
                PRIMARY_CLIP_ANALYSIS.out.trimgalore_log.map{ vec -> vec[1] },
                PRIMARY_CLIP_ANALYSIS.out.bowtie_align_log.map{ vec -> vec[1] },
                PRIMARY_CLIP_ANALYSIS.out.star_align_log_final.map{ vec -> vec[1] },
                PRIMARY_CLIP_ANALYSIS.out.clipqc_log,
                Channel.fromPath("./conf/multiqc_config.yaml")
            )
            .collect()
    )

}
