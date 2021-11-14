//
// Check input samplesheet and get read channels
//

params.options = [:]

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    input_fastq

    main:
    samplesheet
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it, input_fastq) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, fastq_1 ]
def create_fastq_channel(LinkedHashMap row, String input_fastq) {
    def meta = [:]
    meta.id           = row['Sample name']
    meta.genome       = row['mapto']
    meta.barcode      = row["5'barcode"]

    def array = []
    array = [ meta, input_fastq ]

    return array
}