nextflow.enable.dsl=2

include { UMICOLLAPSE } from "../../modules/local/umicollapse/main"

workflow {

    meta = [id:file(params.bam).name, single_end: true]

    println file(params.bam)
    println file(params.bai)

    UMICOLLAPSE ( [meta, file(params.bam), file(params.bai) ] )

}