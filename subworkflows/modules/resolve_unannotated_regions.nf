nextflow.enable.dsl=2

include { RESOLVE_UNANNOTATED_REGIONS } from "../../modules/local/resolve_unannotated_regions/main"

workflow {

    RESOLVE_UNANNOTATED_REGIONS (
        file(params.filtered_regions),
        file(params.unfiltered_regions),
        file(params.gtf),
        file(params.fai),
    )

}