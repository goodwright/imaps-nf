nextflow.enable.dsl=2

include { STAR_GENOMEGENERATE } from "../../modules/nf-core/modules/star/genomegenerate/main"

workflow {

    STAR_GENOMEGENERATE ( file(params.fasta), file(params.gtf) )

}