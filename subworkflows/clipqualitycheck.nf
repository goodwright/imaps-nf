#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CLIPQC } from '../modules/local/clipqc/main.nf'
include { MULTIQC } from '../modules/local/multiqc/main.nf'

def strip_meta = { vec -> vec[1] }

workflow CLIP_QUALITY_CHECK {

    take:
        fastqc_html
        fastqc_zip
        trimgalore_log
        bowtie_align_log
        star_align_log_final
        umicollapse_log
        crosslinks
        icount_peaks
        paraclu_peaks
        clippy_peaks

    main:

    CLIPQC (
        bowtie_align_log.map(strip_meta).collect(),
        star_align_log_final.map(strip_meta).collect(),
        umicollapse_log.map(strip_meta).collect(),
        crosslinks.map(strip_meta).collect(),
        icount_peaks.map(strip_meta).collect(),
        paraclu_peaks.map(strip_meta).collect(),
        clippy_peaks.map(strip_meta).collect()
    )

    MULTIQC (
        Channel
            .empty()
            .concat(
                fastqc_html.map(strip_meta),
                fastqc_zip.map(strip_meta),
                trimgalore_log.map(strip_meta),
                bowtie_align_log.map(strip_meta),
                star_align_log_final.map(strip_meta),
                CLIPQC.out.log,
                Channel.fromPath("$projectDir/../conf/multiqc/clip.yaml")
            )
            .collect()
    )
}
