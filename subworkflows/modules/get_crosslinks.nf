nextflow.enable.dsl=2

include { GET_CROSSLINKS } from "../../modules/local/get_crosslinks/main"

workflow {

    meta = [id:file(params.bam).name, single_end: true]

    GET_CROSSLINKS ( [meta, file(params.bam), file(params.bai) ], file(params.fai) )

}