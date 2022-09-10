nextflow.enable.dsl=2

include { ICOUNT_SIGXLS } from "../../modules/luslab/nf-core-modules/icount/sigxls/main"

workflow {

    meta = [id:file(params.bed).name, single_end: true]

    ICOUNT_SIGXLS ( [meta, file(params.bed)], file(params.segmentation) )

}