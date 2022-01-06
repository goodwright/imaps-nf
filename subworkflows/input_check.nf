//
// Check input samplesheet and get read channels
//

params.options = [:]

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    samplesheet
        .splitCsv ( header:true, sep:',', strip:true)
        .map { create_fastq_channel(it) }
        .set { readsMeta }

    emit:
    readsMeta                                     // channel: [ meta ]
}

// Function to get list of [ meta, fastq_1 ]
def create_fastq_channel(LinkedHashMap row) {
    // Takes a row from a samples CSV file and creates a meta object which
    // describes it.

    def meta = [:]
    meta.id           = row.entrySet().iterator().next().getValue() // This is janky and means sample id always has to come 1st
    meta.genome       = row.Species
    meta.barcode      = row.FivePrimeBarcode
    meta.single_end   = true

    return meta
}