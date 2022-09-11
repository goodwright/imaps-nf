nextflow.enable.dsl=2

include { ICOUNT_RNAMAPS } from "../../modules/local/icount_rnamaps/main"

workflow {

    meta = [id:file(params.bed).name, single_end: true]

    ICOUNT_RNAMAPS ( [meta, file(params.bed)], file(params.regions) )

}