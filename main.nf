#!/usr/bin/env nextflow
/*
========================================================================================
    goodwright/imaps-nf
========================================================================================
    Github : https://github.com/goodwright/imaps-nf
    Website: https://imaps.goodwright.com/
    Slack  : https://join.slack.com/t/imapsgroup/shared_invite/zt-r24y3591-Xbhnym2t38u_urU~I0K0lQ
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta        = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.gtf          = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.gff          = WorkflowMain.getGenomeAttribute(params, 'gff')
params.gene_bed     = WorkflowMain.getGenomeAttribute(params, 'bed12')
params.star_index   = WorkflowMain.getGenomeAttribute(params, 'star')
params.hisat2_index = WorkflowMain.getGenomeAttribute(params, 'hisat2')
params.rsem_index   = WorkflowMain.getGenomeAttribute(params, 'rsem')
params.salmon_index = WorkflowMain.getGenomeAttribute(params, 'salmon')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { DEMUXANDANALYSE } from './workflows/demuxandanalyse'

//
// WORKFLOW: Run main demultiplex and analyse pipeline
//
workflow IMAPS_DEMUXANDANALYSE {
    DEMUXANDANALYSE ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    IMAPS_DEMUXANDANALYSE ()
}

/*
========================================================================================
    THE END
========================================================================================
*/