#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ICOUNT_SEGMENT } from '../../modules/luslab/nf-core-modules/icount/segment/main'

workflow {

    ICOUNT_SEGMENT ( params.gtf, params.fai )
}