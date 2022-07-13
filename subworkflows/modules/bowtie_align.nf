nextflow.enable.dsl=2

include { BOWTIE_ALIGN } from "../../modules/nf-core/modules/bowtie/align/main"

workflow {

    meta = [id:file(params.fastq).name, single_end: true]

    BOWTIE_ALIGN ( [meta, file(params.fastq)], file(params.index) )

}