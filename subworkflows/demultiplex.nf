nextflow.enable.dsl=2

include { XLSX_TO_CSV } from "../modules/local/xlsx_to_csv/main"
include { CSV_TO_BARCODE } from "../modules/local/csv_to_barcode/main"

workflow {

    // Create channel and load the CSV file into it
    if (params.annotation.matches(".*xlsx")) {
        ch_csv = XLSX_TO_CSV ( params.annotation ).csv
    } else {
        ch_csv = Channel.fromPath(file(params.annotation))
    }

    // Create channel and load the multiplexed reads file into it
    ch_multiplexed_fastq = file(params.multiplexed_fastq)

    // Get ultraplex barcodes file from CSV file
    CSV_TO_BARCODE ( ch_csv )

    // Create the meta object describing the reads file needed by Ultraplex
    meta = [id:file(params.multiplexed_fastq).name]

}