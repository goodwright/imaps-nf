nextflow.enable.dsl=2

include { FILTER_TRANSCRIPTOME_BAM } from "../../modules/local/filter_transcriptome_bam/main"

workflow {

    meta = [id:file(params.bam).name, single_end: true]

    FILTER_TRANSCRIPTOME_BAM (
        [meta, file(params.bam)],
        file(params.transcripts)
    )

}