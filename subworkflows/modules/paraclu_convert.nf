nextflow.enable.dsl=2

include { PARACLU_CONVERT } from "../../modules/luslab/nf-core-modules/paraclu/convert/main"

workflow {

    meta = [id:file(params.peaks).name, single_end: true]

    PARACLU_CONVERT ( [meta, file(params.peaks)] )

}