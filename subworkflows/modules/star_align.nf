nextflow.enable.dsl=2

include { STAR_ALIGN } from "../../modules/nf-core/modules/star/align/main"

workflow {

    meta = [id:file(params.fastq).name, single_end: true]

    STAR_ALIGN (
        [meta, file(params.fastq)],
        file(params.index),
        file(params.gtf),
        false,
        "",
        ""
    )

}