nextflow.enable.dsl=2

include { FILTER_GTF } from "../../modules/local/filter_gtf/main"

workflow {

    FILTER_GTF ( file(params.gtf) )

}