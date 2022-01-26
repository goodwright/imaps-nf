#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { TRIMGALORE } from '../modules/nf-core/modules/trimgalore/main'  addParams( options: [:] )
include { BOWTIE_ALIGN } from '../modules/nf-core/modules/bowtie/align/main'    addParams( save_unaligned: true, options: [args:"-v 2 -m 1 --norc --best --strata"] )
include { STAR_ALIGN } from '../modules/nf-core/modules/star/align/main'    addParams( options: [args:"--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM"] )
include { UMITOOLS_DEDUP } from '../modules/nf-core/modules/umitools/dedup/main'    addParams( options: [args:"--umi-separator='rbc:'"] )
include { SAMTOOLS_INDEX as STAR_SAMTOOLS_INDEX} from '../modules/nf-core/modules/samtools/index/main'    addParams( options: [:] )
include { SAMTOOLS_INDEX as UMITOOLS_SAMTOOLS_INDEX} from '../modules/nf-core/modules/samtools/index/main'    addParams( options: [:] )
include { GET_CROSSLINKS } from '../modules/local/get_crosslinks/main'    addParams( options: [:] )
include { CROSSLINKS_COVERAGE } from '../modules/luslab/nf-core-modules/crosslinks/coverage/main'    addParams( options: [:] )
include { CROSSLINKS_NORMCOVERAGE } from '../modules/luslab/nf-core-modules/crosslinks/normcoverage/main'    addParams( options: [:] )

include { FILTER_TRANSCRIPTS } from '../modules/local/filter_transcriptome_bam/main'    addParams( options: [args:"--umi-separator='rbc:'"] )
include { UMITOOLS_DEDUP as TOME_UMITOOLS_DEDUP } from '../modules/nf-core/modules/umitools/dedup/main'    addParams( options: [args:"--umi-separator='rbc:'"] )
include { SAMTOOLS_INDEX as TOME_STAR_SAMTOOLS_INDEX } from '../modules/nf-core/modules/samtools/index/main'    addParams( options: [:] )
include { SAMTOOLS_INDEX as TOME_UMITOOLS_SAMTOOLS_INDEX } from '../modules/nf-core/modules/samtools/index/main'    addParams( options: [:] )
include { GET_CROSSLINKS as TOME_GET_CROSSLINKS } from '../modules/local/get_crosslinks/main'    addParams( options: [:] )
include { CROSSLINKS_COVERAGE as TOME_CROSSLINKS_COVERAGE } from '../modules/luslab/nf-core-modules/crosslinks/coverage/main'    addParams( options: [:] )
include { CROSSLINKS_NORMCOVERAGE as TOME_CROSSLINKS_NORMCOVERAGE } from '../modules/luslab/nf-core-modules/crosslinks/normcoverage/main'    addParams( options: [:] )


/*
include { ICOUNT_SIGXLS } from '../modules/luslab/nf-core-modules/icount/sigxls/main'    addParams( options: [:] )
include { ICOUNT_SUMMARY } from '../modules/local/icount_summary/main'    addParams( options: [:] )
include { ICOUNT_RNAMAPS } from '../modules/local/icount_rnamaps/main'    addParams( options: [:] )
include { ICOUNT_PEAKS } from '../modules/local/icount_peaks/main'    addParams( options: [:] )
include { CLIPPY } from '../modules/luslab/nf-core-modules/clippy/main'    addParams( options: [:] )
include { PARACLU_CONVERT } from '../modules/luslab/nf-core-modules/paraclu/convert/main'    addParams( options: [:] )
include { PARACLU_PARACLU } from '../modules/luslab/nf-core-modules/paraclu/paraclu/main'    addParams( options: [:] )
include { PARACLU_CUT } from '../modules/luslab/nf-core-modules/paraclu/cut/main'    addParams( options: [:] )
*/

workflow {
    // If running straight from command line, will need to construct the
    // [meta, reads] pair channel first
    reads = [[id: params.id, single_end: params.single_end], file(params.fastq)]

    // Now just pass that along with the rest of params
    PRIMARY_ANALYSIS (
        reads,
        params.smrna_genome,
        params.star_index,
        params.gtf,
        params.genome_fai,
        params.icount_regions,
        params.icount_segment,
        params.longest_transcripts,
        params.transcriptome_fai
    )
}

workflow PRIMARY_ANALYSIS {

    take:
        reads
        smrna_genome
        star_index
        gtf
        genome_fai
        icount_regions
        icount_segment
        longest_transcripts
        transcriptome_index

    main:

    // Start things off with TrimGalore
    TRIMGALORE ( reads )

    // Run Bowtie Align on each of the reads, with the provided genome index
    // files as reference
    BOWTIE_ALIGN (
        TRIMGALORE.out.reads,
        file(params.smrna_genome)
    )

    // Run STAR Align on the reads which didn't match above
    STAR_ALIGN (
        BOWTIE_ALIGN.out.fastq,
        file(params.star_index),
        file(params.gtf)
    )

    // Preparing crosslinks from genomic mapping

    STAR_SAMTOOLS_INDEX ( STAR_ALIGN.out.bam_sorted )

    ch_umi_input = STAR_ALIGN.out.bam_sorted.combine(STAR_SAMTOOLS_INDEX.out.bai, by: 0)

    //UMI-TOOLS
    UMITOOLS_DEDUP ( ch_umi_input )

    //SAMTOOLS INDEX the deduped BAM
    UMITOOLS_SAMTOOLS_INDEX ( UMITOOLS_DEDUP.out.bam )

    ch_xl_input = UMITOOLS_DEDUP.out.bam.combine(UMITOOLS_SAMTOOLS_INDEX.out.bai, by: 0)

    //Get crosslinks
    GET_CROSSLINKS (
        ch_xl_input,
        file(params.genome_fai)
    )

    // Get coverage and normalized coverage
    CROSSLINKS_COVERAGE ( GET_CROSSLINKS.out.crosslinkBed )
    CROSSLINKS_NORMCOVERAGE ( GET_CROSSLINKS.out.crosslinkBed )


    // Preparing crosslinks from transcriptome maping

    FILTER_TRANSCRIPTS ( 
        STAR_ALIGN.out.bam_transcript, 
        file(params.longest_transcripts)
    )

    TOME_STAR_SAMTOOLS_INDEX ( FILTER_TRANSCRIPTS.out.filtered_bam )

    tome_ch_umi_input = FILTER_TRANSCRIPTS.out.filtered_bam.combine(TOME_STAR_SAMTOOLS_INDEX.out.bai, by: 0)

    //UMI-TOOLS
    TOME_UMITOOLS_DEDUP ( tome_ch_umi_input )

    //SAMTOOLS INDEX the deduped BAM
    TOME_UMITOOLS_SAMTOOLS_INDEX ( TOME_UMITOOLS_DEDUP.out.bam )

    tome_ch_xl_input = TOME_UMITOOLS_DEDUP.out.bam.combine(TOME_UMITOOLS_SAMTOOLS_INDEX.out.bai, by: 0)

    //Get crosslinks
    TOME_GET_CROSSLINKS (
        tome_ch_xl_input,
        file(params.transcriptome_fai)
    )

    // Get coverage and normalized coverage
    TOME_CROSSLINKS_COVERAGE ( TOME_GET_CROSSLINKS.out.crosslinkBed )
    TOME_CROSSLINKS_NORMCOVERAGE ( TOME_GET_CROSSLINKS.out.crosslinkBed )

    // Peak calling, summary statistics, RNA-maps and PEKA
    /*
    //iCount summary
    ICOUNT_SUMMARY (
        GET_CROSSLINKS.out.crosslinkBed,
        file(params.icount_regions)
    )

    //iCount RNA-maps
    ICOUNT_RNAMAPS (
        GET_CROSSLINKS.out.crosslinkBed,
        file(params.icount_regions)
    )

    // Run peak callers - starting with Paraclu
    PARACLU_PARACLU ( GET_CROSSLINKS.out.crosslinkBed )
    PARACLU_CUT ( PARACLU_PARACLU.out.sigxls )
    PARACLU_CONVERT ( PARACLU_CUT.out.peaks )

    //ICOUNT SIGXLS
    ICOUNT_SIGXLS (
        GET_CROSSLINKS.out.crosslinkBed,
        file(params.icount_segment)
    )

    ch_icount_peaks = GET_CROSSLINKS.out.crosslinkBed.combine(ICOUNT_SIGXLS.out.sigxls, by: 0)

    //ICOUNT PEAKS
    ICOUNT_PEAKS ( ch_icount_peaks )

    //CLIPPY
    CLIPPY (
        GET_CROSSLINKS.out.crosslinkBed,
        file(params.gtf),
        file(params.genome_fai)
    )
    */
}