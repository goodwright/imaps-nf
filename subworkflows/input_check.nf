//
// Check input samplesheet and get read channels
//

params.options = [:]

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    samplesheet
        .splitCsv ( header:true, sep:',' )
        //.println it.toString()
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, fastq_1 ]
def create_fastq_channel(LinkedHashMap row) {

    def meta = [:]
    meta.id           = row[Sample Name]
    println meta.id
    meta.genome       = row.Species
    meta.barcode      = row["5'barcode"]
    
    fastq_name = "ultraplex_demux_" + meta.id + "_fastq.gz"
    def array = []
    array = [ meta, [fastq_name] ]

    println Arrays.toString(array)
    return array
}