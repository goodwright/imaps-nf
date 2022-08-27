nextflow.enable.dsl=2

include { PARACLU_CUT } from "../../modules/luslab/nf-core-modules/paraclu/cut/main"

workflow {

    meta = [id:file(params.sigxls).name, single_end: true]

    PARACLU_CUT ( [meta, file(params.sigxls)] )

}