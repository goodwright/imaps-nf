nextflow.enable.dsl=2

include { FIND_LONGEST_TRANSCRIPT } from "../../modules/local/find_longest_transcript/main"

workflow {

    FIND_LONGEST_TRANSCRIPT ( file(params.gtf) )

}