#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ULTRAPLEX } from '../modules/luslab/nf-core-modules/ultraplex/main'    addParams( options: [:] )
include { CSV_TO_BARCODE } from '../modules/local/csv_to_barcode/main'    addParams( options: [:] )
include { FASTQC } from '../modules/nf-core/modules/fastqc/main' addParams( options: [:] )

workflow {
    DEMULTIPLEX ( params.csv, params.multiplexed_fastq )
}

workflow DEMULTIPLEX {

    take:
        csv
        multiplexed_fastq
    
    main:
    // Create channel for CSV file
    ch_csv = Channel.fromPath(csv)

    // Get ultraplex barcodes file from CSV file
    CSV_TO_BARCODE ( ch_csv )

    // Create channel for multiplexed reads file
    ch_multiplexed_fastq = file(multiplexed_fastq)

    // Create the meta object required for some of the following processes
    def initialMeta = [:]
    initialMeta.id = "test-string"

    // Run Ultraplex on the reads and barcodes
    ULTRAPLEX (
        [initialMeta, ch_multiplexed_fastq], 
        CSV_TO_BARCODE.out
    )

    // Create a channel which produces a meta object for each row in the CSV
    ch_csv
    .splitCsv ( header:true, sep:',', strip:true)
    .map { create_fastq_channel(it) }
    .set { readsMeta }

    // Create a new channel from the FASTQ output in which every entry is a
    // tuple, with an ID string and a fastq path.
    ULTRAPLEX.out.fastq
    .flatten()
    .map { it -> ["id":it.toString().replaceAll(/.*ultraplex_demux_/,"")
        .replaceAll(/\.fastq\.gz/,""),"fastq":it] }
    .set { ch_demuxed_reads }

    // Create a channel combining the CSV meta objects with the reads
    readsMeta
    .flatten()
    .cross(ch_demuxed_reads)
    .map { meta, fastq -> [ meta, fastq.fastq ] }
    .set {ch_reads_with_meta}

    // Run FASTQC on each of the reads
    FASTQC ( ch_reads_with_meta )

    emit:
      ULTRAPLEX.out.fastq

}



def create_fastq_channel(LinkedHashMap row) {
    /** Takes a row from a samples CSV file and creates a meta object which
        describes it.
    */

    def meta = [:]
    meta.id           = row.entrySet().iterator().next().getValue() // This is janky and means sample id always has to come 1st
    meta.genome       = row.Species
    meta.barcode      = row.FivePrimeBarcode
    meta.single_end   = true
    return meta
}