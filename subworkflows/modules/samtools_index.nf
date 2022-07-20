nextflow.enable.dsl=2

include { SAMTOOLS_INDEX } from "../../modules/nf-core/modules/samtools/index/main"

workflow {

    SAMTOOLS_INDEX ( [[id:file(params.bam).name], file(params.bam)] )

}