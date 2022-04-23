#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DEMUXANDANALYSE } from './workflows/demuxandanalyse'

workflow {
    DEMUXANDANALYSE ()
}