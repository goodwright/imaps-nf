nextflow.enable.dsl=2

include { PARACLU_PARACLU } from "../../modules/luslab/nf-core-modules/paraclu/paraclu/main"

workflow {

    meta = [id:file(params.crosslinks).name, single_end: true]

    PARACLU_PARACLU ( [meta, file(params.crosslinks)] )

}