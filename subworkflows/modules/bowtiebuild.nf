#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { BOWTIE_BUILD } from '../../modules/nf-core/modules/bowtie/build/main'

workflow {

    BOWTIE_BUILD ( params.fasta )

}