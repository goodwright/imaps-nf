/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

include { ULTRAPLEX } from '../modules/luslab/nf-core-modules/ultraplex'    addParams( options: [:] )
include { FASTQC } from '../modules/nf-core/fastqc' addParams( genome_options: publish_genome_options, index_options: publish_index_options, gffread_options: gffread_options,  star_index_options: star_genomegenerate_options,  hisat2_index_options: hisat2_build_options, rsem_index_options: rsem_preparereference_options, salmon_index_options: salmon_index_options )
include { TRIMGALORE } from '../modules/nf-core/trimgalore'  addParams( calculateexpression_options: rsem_calculateexpression_options, samtools_sort_options: samtools_sort_genome_options, samtools_index_options: samtools_index_genome_options, samtools_stats_options: samtools_index_genome_options, merge_counts_options: modules['rsem_merge_counts'] )
include { BOWTIE } from '../modules/nf-core/bowtie/align'    addParams( genome_options: publish_genome_options, tximport_options: modules['star_salmon_tximport'], salmon_quant_options: modules['star_salmon_quant'], merge_counts_options: modules['star_salmon_merge_counts'] )
include { STAR_ALIGN } from '../modules/nf-core/star/align'    addParams( genome_options: publish_genome_options, tximport_options: modules['salmon_tximport'], salmon_quant_options: salmon_quant_options, merge_counts_options: modules['salmon_merge_counts'] )
