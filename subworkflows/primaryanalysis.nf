#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { TRIMGALORE } from '../modules/nf-core/modules/trimgalore/main'  addParams( options: [:] )
include { BOWTIE_ALIGN } from '../modules/nf-core/modules/bowtie/align/main'    addParams( save_unaligned: true, options: [args:"-v 2 -m 1 --norc --best --strata"] )
include { STAR_ALIGN } from '../modules/nf-core/modules/star/align/main'    addParams( options: [args:"--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate"] )
include { UMITOOLS_DEDUP } from '../modules/nf-core/modules/umitools/dedup/main'    addParams( options: [args:"--umi-separator='rbc:'"] )
include { SAMTOOLS_INDEX as STAR_SAMTOOLS_INDEX} from '../modules/nf-core/modules/samtools/index/main'    addParams( options: [:] )
include { SAMTOOLS_INDEX as UMITOOLS_SAMTOOLS_INDEX} from '../modules/nf-core/modules/samtools/index/main'    addParams( options: [:] )
include { GET_CROSSLINKS } from '../modules/local/get_crosslinks/main'    addParams( options: [:] )
include { CROSSLINKS_COVERAGE } from '../modules/luslab/nf-core-modules/crosslinks/coverage/main'    addParams( options: [:] )
include { CROSSLINKS_NORMCOVERAGE } from '../modules/luslab/nf-core-modules/crosslinks/normcoverage/main'    addParams( options: [:] )
include { ICOUNT_SIGXLS } from '../modules/luslab/nf-core-modules/icount/sigxls/main'    addParams( options: [:] )
include { ICOUNT_SUMMARY } from '../modules/local/icount_summary/main'    addParams( options: [:] )
include { ICOUNT_RNAMAPS } from '../modules/local/icount_rnamaps/main'    addParams( options: [:] )
include { ICOUNT_PEAKS } from '../modules/local/icount_peaks/main'    addParams( options: [:] )
include { CLIPPY } from '../modules/luslab/nf-core-modules/clippy/main'    addParams( options: [:] )
include { PARACLU_CONVERT } from '../modules/luslab/nf-core-modules/paraclu/convert/main'    addParams( options: [:] )
include { PARACLU_PARACLU } from '../modules/luslab/nf-core-modules/paraclu/paraclu/main'    addParams( options: [:] )
include { PARACLU_CUT } from '../modules/luslab/nf-core-modules/paraclu/cut/main'    addParams( options: [:] )

workflow {
    PRIMARY_ANALYSIS (
        params.id,
        params.species,
        params.barcode,
        params.single_end,
        params.fastq,
        params.smrna_genome,
        params.star_index,
        params.gtf,
        params.genome_fai,
        params.icount_regions,
        params.icount_segment,
    )
}

workflow PRIMARY_ANALYSIS {

    take:
        id
        species
        barcode
        single_end
        fastq
        smrna_genome
        star_index
        gtf
        genome_fai
        icount_regions
        icount_segment
        

    main:
    // Create channel for demultiplexed fastq file
    ch_fastq = file(fastq)

    // Create the meta object that will be passed throughout the analysis
    meta = [:]
    meta.id           = id
    meta.genome       = species
    meta.barcode      = barcode
    meta.single_end   = single_end

    // Start things off with TrimGalore
    TRIMGALORE ( [meta, ch_fastq] )

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
}