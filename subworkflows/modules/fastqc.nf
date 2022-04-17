#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { FASTQC } from '../../modules/nf-core/modules/fastqc/main'

workflow {

    meta = [:]
    meta.id = file(params.fastq).name

    reads = [[
        id: params.fastq.split("/")[-1].replace(".gz", "")
            .replace(".fastq", "").replace("ultraplex_demux_", ""),
        single_end: true
    ], file(params.fastq)]

    FASTQC ( reads )
}