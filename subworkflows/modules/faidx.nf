#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { SAMTOOLS_FAIDX } from '../../modules/nf-core/modules/samtools/faidx/main'

workflow {

    SAMTOOLS_FAIDX ( [[:], params.fasta] )

}