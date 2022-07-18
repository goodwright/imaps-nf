process FILTER_TRANSCRIPTOME_BAM {
    tag "$meta.id"
    label "process_medium"

    container "quay.io/biocontainers/samtools:1.14--hb421002_0"
    
    input:
    tuple val(meta), path(transcriptome_bam)
    path(transcripts)

    output:
    tuple val(meta), path("*_filtered_transcriptome.bam"), emit: filtered_bam

    script:
    def prefix = "${meta.id}"

    //SHELL
    """
    samtools sort $transcriptome_bam > sorted.bam
    samtools index sorted.bam
    samtools view -h sorted.bam `cat $transcripts` > filtunsort.bam
    samtools sort filtunsort.bam > ${prefix}_filtered_transcriptome.bam
    """
}