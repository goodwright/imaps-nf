nextflow.enable.dsl=2

include { ICOUNT_SEGMENT } from "../../modules/luslab/nf-core-modules/icount/segment/main"

workflow {

    ICOUNT_SEGMENT ( file(params.gtf), file(params.fai) )

}