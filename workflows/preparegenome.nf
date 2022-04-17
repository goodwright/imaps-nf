#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { 
    GUNZIP as DNA_GUNZIP
    GUNZIP as RNA_GUNZIP       
                                      } from '../modules/nf-core/modules/gunzip/main'
include { STAR_GENOMEGENERATE         } from '../modules/nf-core/modules/star/genomegenerate/main'
include { BOWTIE_BUILD                } from '../modules/nf-core/modules/bowtie/build/main'
include { SAMTOOLS_FAIDX              } from '../modules/nf-core/modules/samtools/faidx/main'
include {
    ICOUNT_SEGMENT as RAW_ICOUNT_SEGMENT;
    ICOUNT_SEGMENT as FILTERED_ICOUNT_SEGMENT
                                      } from '../modules/luslab/nf-core-modules/icount/segment/main'
include { FIND_LONGEST_TRANSCRIPT          } from '../modules/local/find_longest_transcript/main'
include { FILTER_GTF                  } from '../modules/local/filter_gtf/main'
include {
    RESOLVE_UNANNOTATED;
    RESOLVE_UNANNOTATED as RESOLVE_UNANNOTATED_GENIC_OTHER
                                      } from '../modules/local/resolve_unannotated/main'

workflow {

    // If smRNA genome is compressed, uncompress it
    if (params.smrna_fasta.matches(".*gz")) {
        RNA_GUNZIP ( [[:], file(params.smrna_fasta)] )
        RNA_GUNZIP.out.gunzip.map { it -> it[1]}.set { ch_smrna_fasta }
    } else {
        ch_smrna_fasta = file(params.smrna_fasta)
    }

    // Create a genome index from the smRNA
    BOWTIE_BUILD ( ch_smrna_fasta )

    // Determine the longest CDS transcripts from the input GTF
    FIND_LONGEST_TRANSCRIPT ( params.gtf )

    // Filter the GTF file as required
    FILTER_GTF ( params.gtf )
    ch_post_filter_gtf = FILTER_GTF.out.post_filtering_gtf

    // If genome is compressed, uncompress it
    if (params.fasta.matches(".*gz")) {
        DNA_GUNZIP ( [[:], file(params.fasta)] )
        ch_fasta_with_meta = DNA_GUNZIP.out.gunzip
        ch_fasta_with_meta.map { it -> it[1]}.set { ch_fasta }
    } else {
        ch_fasta = file(params.fasta)
    }

    // Create a FAI genome index using samtools
    ch_fai_with_meta = SAMTOOLS_FAIDX ( ch_fasta_with_meta ).fai
    ch_fai_with_meta.map { it -> it[1]}.set { ch_fai }

    // Index the genome using STAR module
    STAR_GENOMEGENERATE ( ch_fasta, params.gtf )

    // Run iCount-Segment on raw GTF
    RAW_ICOUNT_SEGMENT ( params.gtf, ch_fai )

    // Run iCount-Segment on filtered GTF
    FILTERED_ICOUNT_SEGMENT ( ch_post_filter_gtf, ch_fai )



    // These two processes are likely getting the wrong inputs currently
    RESOLVE_UNANNOTATED (
        FILTERED_ICOUNT_SEGMENT.out.regions, // filtered_segmentation
        RAW_ICOUNT_SEGMENT.out.regions,      // unfiltered_segmentation
        ch_post_filter_gtf,                  // gtf
        ch_fai,                              // fai
    )
    RESOLVE_UNANNOTATED_GENIC_OTHER (
        FILTERED_ICOUNT_SEGMENT.out.regions, // filtered_segmentation
        RAW_ICOUNT_SEGMENT.out.regions,      // unfiltered_segmentation
        ch_post_filter_gtf,                  // gtf
        ch_fai,                              // fai
    )
}