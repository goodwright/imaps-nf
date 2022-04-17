#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { GUNZIP as DNA_GUNZIP        } from '../modules/nf-core/modules/gunzip/main'
include { GUNZIP as RNA_GUNZIP        } from '../modules/nf-core/modules/gunzip/main'
include { STAR_GENOMEGENERATE         } from '../modules/nf-core/modules/star/genomegenerate/main'
include { BOWTIE_BUILD                } from '../modules/nf-core/modules/bowtie/build/main'
include { SAMTOOLS_FAIDX              } from '../modules/nf-core/modules/samtools/faidx/main'
include {
    ICOUNT_SEGMENT as RAW_ICOUNT_SEGMENT;
    ICOUNT_SEGMENT as FILTERED_ICOUNT_SEGMENT
                                      } from '../modules/luslab/nf-core-modules/icount/segment/main'
include { LONGEST_TRANSCRIPT          } from '../modules/local/find_longest_transcript/main'
include { FILTER_GTF                  } from '../modules/local/filter_gtf/main'
include {
    RESOLVE_UNANNOTATED;
    RESOLVE_UNANNOTATED as RESOLVE_UNANNOTATED_GENIC_OTHER
                                      } from '../modules/local/resolve_unannotated/main'

workflow {

    // If genome is compressed, uncompress it
    if (params.fasta.matches(".*gz")) {
        ch_fasta = DNA_GUNZIP ( file(params.fasta) ).gunzip
    } else {
        ch_fasta = file(params.fasta)
    }

    // If smrna genome is compressed, uncompress it
    if (params.smrna_fasta.matches(".*gz")) {
        ch_smrna_fasta = RNA_GUNZIP ( file(params.smrna_fasta) ).gunzip
    } else {
        ch_smrna_fasta = file(params.smrna_fasta)
    }

    // Index the genome using STAR module
    STAR_GENOMEGENERATE (
        ch_fasta,
        params.gtf
    )

    BOWTIE_BUILD (
        ch_smrna_fasta
    )

    // Create a FAI genome index using samtools
    ch_fai = SAMTOOLS_FAIDX ( ch_fasta ).fai

    FILTER_GTF ( params.gtf )

    // iCount segment the raw annotation
    RAW_ICOUNT_SEGMENT (
        params.gtf,
        ch_fai
    )

    // iCount segment the filtered annotation
    FILTERED_ICOUNT_SEGMENT (
        FILTER_GTF.out.post_filtering_gtf,
        ch_fai
    )

    // Resolve unannotated for iCount summary and other processes
    RESOLVE_UNANNOTATED (
        FILTERED_ICOUNT_SEGMENT.out.regions, // filtered_segmentation
        RAW_ICOUNT_SEGMENT.out.regions,      // unfiltered_segmentation
        FILTER_GTF.out.post_filtering_gtf,   // gtf
        ch_fai,                              // fai
    )

    // Resolve unannotated for PEKA
    RESOLVE_UNANNOTATED_GENIC_OTHER (
        FILTERED_ICOUNT_SEGMENT.out.regions, // filtered_segmentation
        RAW_ICOUNT_SEGMENT.out.regions,      // unfiltered_segmentation
        FILTER_GTF.out.post_filtering_gtf,   // gtf
        ch_fai,                              // fai
    )

    // Find the longest CDS transcript per gene.
    LONGEST_TRANSCRIPT (
        params.gtf
    )

}