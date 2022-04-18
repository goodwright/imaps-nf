process GET_CROSSLINKS {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/mulled-v2-c8623b4f6522dddd48913bd12dcf405d1d4f0ce1:10e4c359b727e884f6e19ee978f89c44dbaca255-0'

    input:
      tuple val(meta), path(bam), path(bai)
      path(fai)

    output:
      tuple val(meta), path ("*bed.gz"), emit: crosslinkBed

    script:
      def prefix    = "${meta.id}"

      //SHELL
      """
      bedtools bamtobed -i $bam > dedupe.bed
      bedtools shift -m 1 -p -1 -i dedupe.bed -g $fai > shifted.bed
      bedtools genomecov -dz -strand + -5 -i shifted.bed -g $fai | awk '{OFS="\\t"}{print \$1, \$2, \$2+1, ".", \$3, "+"}' > pos.bed
      bedtools genomecov -dz -strand - -5 -i shifted.bed -g $fai | awk '{OFS="\\t"}{print \$1, \$2, \$2+1, ".", \$3, "-"}' > neg.bed
      cat pos.bed neg.bed | sort -k1,1 -k2,2n | pigz > ${prefix}.bed.gz
      """
}