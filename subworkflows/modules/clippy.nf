nextflow.enable.dsl=2

include { CLIPPY } from "../../modules/luslab/nf-core-modules/clippy/main"

workflow {

    meta = [id:file(params.crosslinks).name, single_end: true]

    CLIPPY (
        [meta, file(params.crosslinks)],
        file(params.gtf),
        file(params.fai),
    )

}