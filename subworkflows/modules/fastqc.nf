#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { FASTQC } from "../../modules/nf-core/modules/fastqc/main";

workflow {

    // FASTQC module takes a [meta, reads] pair as its input. This creates that
    // tuple, with the meta object taking the ID from the filename.
    reads = [[
        id: params.fastq.split("/")[-1].replace(".gz", "")
            .replace(".fastq", "").replace("ultraplex_demux_", ""),
        single_end: true
    ], file(params.fastq)]

    FASTQC ( reads )
}