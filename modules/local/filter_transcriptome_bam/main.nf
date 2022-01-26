#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FILTER_TRANSCRIPTS {
    tag "$meta.id"
    label "process_medium"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0' :
        'quay.io/biocontainers/samtools:1.14--hb421002_0' }"
    
    input:
    tuple val(meta), path(transcriptome_bam)
    path(transcripts)

    output:
    tuple val(meta), path("*.bam"), emit: filtered_bam

    script:
      def prefix    = "${meta.id}"

    //SHELL
    """
    samtools sort $transcriptome_bam > sorted.bam
    samtools index sorted.bam
    samtools view -h sorted.bam `cat $transcripts` > ${prefix}.bam
    """
}