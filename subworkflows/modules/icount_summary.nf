nextflow.enable.dsl=2

include { ICOUNT_SUMMARY } from "../../modules/local/icount_summary/main"

workflow {

    meta = [id:file(params.bed).name, single_end: true]

    ICOUNT_SUMMARY ( [meta, file(params.bed)], file(params.regions) )

}