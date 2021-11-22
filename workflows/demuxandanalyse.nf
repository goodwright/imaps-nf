#!/usr/bin/env nextflow

/* 
DEMULTIPLEX & ANALYSE
*/

nextflow.enable.dsl=2

include { ULTRAPLEX } from '../modules/luslab/nf-core-modules/ultraplex/main'    addParams( options: [:] )
include { CSV_TO_BARCODE } from '../modules/local/csv_to_barcode/main'    addParams( options: [:] )
include { INPUT_CHECK } from '../subworkflows/input_check'    addParams( options: [:] )
include { FASTQC } from '../modules/nf-core/modules/fastqc/main' addParams( options: [:] )
include { TRIMGALORE } from '../modules/nf-core/modules/trimgalore/main'  addParams( options: [:] )
include { BOWTIE_ALIGN } from '../modules/nf-core/modules/bowtie/align/main'    addParams( options: [:] )
include { STAR_ALIGN } from '../modules/nf-core/modules/star/align/main'    addParams( options: [:] )


workflow {
    
// Initialise channels
    ch_software_versions = Channel.empty()

//demultiplexing

    ch_input_meta = Channel.fromPath(params.input)
    ch_input_fasta = file(params.multiplexed_fastq)


    def metaid = [:]
    metaid.id           = "test-string"

    CSV_TO_BARCODE (
        ch_input_meta
    )

    ULTRAPLEX (
        [metaid, ch_input_fasta], 
        CSV_TO_BARCODE.out
    )

    INPUT_CHECK (
         ch_input_meta
    )


ULTRAPLEX.out.fastq
    .flatten()
    .filter(~/.*ultraplex.*/) // get rid of dummy meta.id
    .map { it -> ["id":it.toString().replaceAll(/.*ultraplex_demux_/,"").replaceAll(/\.fastq\.gz/,""),"fastq":it] }
    .set {demuxed_reads}

INPUT_CHECK.out.readsMeta
    .flatten()
    .cross(demuxed_reads)
    .set {ch_demuxed_reads}

ch_demuxed_reads
    .map {
        meta, fastq ->
            [ meta, fastq.fastq ] }
    .set {ch_demuxed_readz}

//fastqc
    FASTQC (
        ch_demuxed_readz
    )

//trim-galore
    TRIMGALORE (
        ch_demuxed_readz
    )

//bowtie to small RNA
    BOWTIE_ALIGN (
        TRIMGALORE.out.reads,
        file(params.smrna_genome)
    )
//TRIMGALORE.out.reads
//unmapped to STAR GENOME
//UMI-TOOLS
//GET CROSSLINKS
//ICOUNT SUMMARY
//ICOUNT RNAMAPS
//GET COVERAGE
//PEKA
//PARACLU
//ICOUNT PEAKS
//CLIPPY
}