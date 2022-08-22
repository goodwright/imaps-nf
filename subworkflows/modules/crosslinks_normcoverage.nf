nextflow.enable.dsl=2

include { CROSSLINKS_NORMCOVERAGE } from "../../modules/luslab/nf-core-modules/crosslinks/normcoverage/main"

workflow {

    meta = [id:file(params.crosslinks).name, single_end: true]

    CROSSLINKS_NORMCOVERAGE ( [meta, file(params.crosslinks)] )

}