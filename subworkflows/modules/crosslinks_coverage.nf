nextflow.enable.dsl=2

include { CROSSLINKS_COVERAGE } from "../../modules/luslab/nf-core-modules/crosslinks/coverage/main"

workflow {

    meta = [id:file(params.crosslinks).name, single_end: true]

    CROSSLINKS_COVERAGE ( [meta, file(params.crosslinks)] )

}