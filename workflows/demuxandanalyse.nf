#!/usr/bin/env nextflow

/* 
DEMULTIPLEX & ANALYSE
*/

nextflow.enable.dsl=2

include { ULTRAPLEX } from '../modules/luslab/nf-core-modules/ultraplex/main'    addParams( options: [:] )
include { CSV_TO_BARCODE } from '../modules/local/csv_to_barcode/main'    addParams( options: [:] )
include { INPUT_CHECK } from '../subworkflows/input_check'    addParams( options: [:] )
include { FASTQC } from '../modules/nf-core/modules/fastqc/main' addParams( options: [:] )
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

workflow {
    
// Initialise channels
    ch_software_versions = Channel.empty()

//demultiplexing

    ch_input_meta = Channel.fromPath(params.input)
    ch_input_fasta = file(params.multiplexed_fastq)


    def metaid = [:]
    metaid.id           = "test-string"

    CSV_TO_BARCODE (
        ch_input_meta
    )

    ULTRAPLEX (
        [metaid, ch_input_fasta], 
        CSV_TO_BARCODE.out
    )

    INPUT_CHECK (
         ch_input_meta
    )

ULTRAPLEX.out.fastq
    .flatten()
    .filter(~/.*ultraplex.*/) // get rid of dummy meta.id
    .map { it -> ["id":it.toString().replaceAll(/.*ultraplex_demux_/,"").replaceAll(/\.fastq\.gz/,""),"fastq":it] }
    .set {demuxed_reads}

INPUT_CHECK.out.readsMeta
    .flatten()
    .cross(demuxed_reads)
    .set {ch_demuxed_reads}

ch_demuxed_reads
    .map {
        meta, fastq ->
            [ meta, fastq.fastq ] }
    .set {ch_demuxed_readz}

//fastqc
    FASTQC (
        ch_demuxed_readz
    )

//trim-galore
    TRIMGALORE (
        ch_demuxed_readz
    )

//bowtie to small RNA
    BOWTIE_ALIGN (
        TRIMGALORE.out.reads,
        file(params.smrna_genome)
    )

//unmapped reads to STAR GENOME
    STAR_ALIGN (
        BOWTIE_ALIGN.out.fastq,
        file(params.star_index),
        file(params.gtf)
    )

// Index the BAM
    STAR_SAMTOOLS_INDEX (
        STAR_ALIGN.out.bam_sorted
    )

ch_umi_input = STAR_ALIGN.out.bam_sorted.combine(STAR_SAMTOOLS_INDEX.out.bai, by: 0)

//UMI-TOOLS
    UMITOOLS_DEDUP (
        ch_umi_input
    )

//SAMTOOLS INDEX THE DEDUPED BAM
    UMITOOLS_SAMTOOLS_INDEX (
        UMITOOLS_DEDUP.out.bam
    )

ch_xl_input = UMITOOLS_DEDUP.out.bam.combine(UMITOOLS_SAMTOOLS_INDEX.out.bai, by: 0)

//GET CROSSLINKS
    GET_CROSSLINKS (
        ch_xl_input,
        file(params.genome_fai)
    )

//GET COVERAGE AND NORMALIZED COVERAGE
    CROSSLINKS_COVERAGE (
        GET_CROSSLINKS.out.crosslinkBed
    )

    CROSSLINKS_NORMCOVERAGE (
        GET_CROSSLINKS.out.crosslinkBed
    )

//ICOUNT SUMMARY
    ICOUNT_SUMMARY (
        GET_CROSSLINKS.out.crosslinkBed,
        file(params.icount_regions)
    )

//ICOUNT RNAMAPS
    ICOUNT_RNAMAPS (
        GET_CROSSLINKS.out.crosslinkBed,
        file(params.icount_regions)
    )

/**
 * Peak Callers *
 */

//PARACLU

//ICOUNT SIGXLS
    ICOUNT_SIGXLS (
        GET_CROSSLINKS.out.crosslinkBed,
        file(params.icount_segment)
    )

ch_icount_peaks = GET_CROSSLINKS.out.crosslinkBed.combine(ICOUNT_SIGXLS.out.sigxls, by: 0)

//ICOUNT PEAKS
    ICOUNT_PEAKS (
        ch_icount_peaks
    )

//CLIPPY


/**
 * Post-peak calling analysis *
 */

//PEKA

}