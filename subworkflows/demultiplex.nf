#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ULTRAPLEX } from '../modules/luslab/nf-core-modules/ultraplex/main'
include { CSV_TO_BARCODE } from '../modules/local/csv_to_barcode/main'
include { XLSX_TO_CSV } from '../modules/local/xlsx_to_csv/main'
include { FASTQC } from '../modules/nf-core/modules/fastqc/main'

workflow {
    DEMULTIPLEX (
        params.annotation,
        params.multiplexed_fastq
    )
}

workflow DEMULTIPLEX {

    take:
        annotation
        multiplexed_fastq
    
    main:

    // Create channel and load the CSV file into it
    if (annotation.matches(".*xlsx")) {
        ch_csv = XLSX_TO_CSV ( annotation ).csv
    } else {
        ch_csv = Channel.fromPath(annotation)
    }

    // Create channel and load the multiplexed reads file into it
    ch_multiplexed_fastq = file(multiplexed_fastq)

    // Get ultraplex barcodes file from CSV file
    CSV_TO_BARCODE ( ch_csv )

    // Create the meta object describing the reads file needed by Ultraplex
    meta = [:]
    meta.id = file(multiplexed_fastq).name

    // Run Ultraplex on the reads and barcodes
    ULTRAPLEX (
        [meta, ch_multiplexed_fastq], 
        CSV_TO_BARCODE.out
    )

    // Ultraplex's fastq output channel produces a tuple containing:
    // (1) the meta object it was passed
    // (2) a list of produced fastq files
    // Create a new channel with just the reads files loaded into it, which are
    // output one by one
    ULTRAPLEX.out.fastq
    .map { it -> it[1]}
    .flatten()
    .set { ch_demultiplexed_reads }

    // Create channel for meta objects and pass meta objects from the CSV file
    // into it one by one.
    ch_csv
    .splitCsv ( header: true, sep:',', strip:true)
    .map { create_fastq_channel(it) }
    .set { ch_reads_meta }

    // For every file that comes out of demultiplexed_reads, map it to a tuple
    // where the first item is an ID (from the filename) and the second is the
    // reads file itself
    ch_demultiplexed_reads.map { it -> [
        "id":it.toString()
        .replaceAll(/.*ultraplex_demux_/,"")
        .replaceAll(/\.fastq\.gz/,""),
        "fastq":it
    ] }.set { ch_reads_with_id }

    // Combine the above two channels to create a channel loaded with tuples of
    // [meta, fastq.gz], which is the form that FASTQC needs them
    ch_reads_meta
    .cross(ch_reads_with_id)
    .map { meta, fastq -> [ meta, fastq.fastq ] }
    .set { ch_reads_with_meta }

    // Run FASTQC on each of the meta-reads pairs
    FASTQC ( ch_reads_with_meta )

    emit:
        fastq       = ch_reads_with_meta
        fastqc_html = FASTQC.out.html
        fastqc_zip  = FASTQC.out.zip

}



def create_fastq_channel(LinkedHashMap row) {
    /** Takes a row from a samples CSV file and creates a meta object which
        describes it.
    */

    def meta = [:]
    meta.id           = row.entrySet().iterator().next().getValue()
    meta.single_end   = true
    meta.species      = row.Species
    meta.pipeline      = row.Pipeline
    return meta
}