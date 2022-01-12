#!/usr/bin/env nextflow

/* 
DEMULTIPLEX & ANALYSE
*/


/* 
RESOLVE PIPELINE FLOW
*/

// transcriptomic TRUE/FALSE
// clippy TRUE/FALSE (followed by peka using clippy peaks)
// icount TRUE/FALSE (followed by peka using icount peaks)
// paraclu TRUE/FALSE (followed by peka using paraclu peaks)
// which-workflow? group, primary, demultiplex-and-analyse (reseq handled above this pipeline)


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
include { CLIPPY } from '../modules/luslab/nf-core-modules/clippy/main'    addParams( options: [:] )
include { PARACLU_CONVERT } from '../modules/luslab/nf-core-modules/paraclu/convert/main'    addParams( options: [:] )
include { PARACLU_PARACLU } from '../modules/luslab/nf-core-modules/paraclu/paraclu/main'    addParams( options: [:] )
include { PARACLU_CUT } from '../modules/luslab/nf-core-modules/paraclu/cut/main'    addParams( options: [:] )


workflow {
    
    // Create channel for CSV file
    ch_csv = Channel.fromPath(params.csv)

    // Get ultraplex barcodes file from CSV file
    CSV_TO_BARCODE ( ch_csv )

    // Create channel for multiplexed reads file
    ch_multiplexed_fastq = file(params.multiplexed_fastq)

    // Create the meta object required for some of the following processes
    def initialMeta = [:]
    initialMeta.id = "test-string"

    // Run Ultraplex on the reads and barcodes
    ULTRAPLEX (
        [initialMeta, ch_multiplexed_fastq], 
        CSV_TO_BARCODE.out
    )

    // Create a channel which produces a meta object for each row in the CSV
    ch_csv
    .splitCsv ( header:true, sep:',', strip:true)
    .map { create_fastq_channel(it) }
    .set { readsMeta }

    // Create a new channel from the FASTQ output in which every entry is a
    // tuple, with an ID string and a fastq path.
    ULTRAPLEX.out.fastq
    .flatten()
    .map { it -> ["id":it.toString().replaceAll(/.*ultraplex_demux_/,"")
        .replaceAll(/\.fastq\.gz/,""),"fastq":it] }
    .set { ch_demuxed_reads }

    // Create a channel combining the CSV meta objects with the reads
    readsMeta
    .flatten()
    .cross(ch_demuxed_reads)
    .map { meta, fastq -> [ meta, fastq.fastq ] }
    .set {ch_reads_with_meta}

    // Run FASTQC on each of the reads
    FASTQC ( ch_reads_with_meta )

    // Run TrimGalore on each of the reads
    TRIMGALORE ( ch_reads_with_meta )

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




 /*    ch_reads_with_meta
    .map {
        meta, fastq ->
            [ meta, fastq.fastq ] }
    .set {ch_demuxed_readz}
    println ch_demuxed_readz.view() */


//demultiplexing

    /* ch_input_fasta = file(params.multiplexed_fastq)


    def metaid = [:]
    metaid.id           = "test-string"


    ULTRAPLEX (
        [metaid, ch_input_fasta], 
        CSV_TO_BARCODE.out
    )

    INPUT_CHECK (
         ch_input_meta
    )

ULTRAPLEX.out.fastq
    .flatten() */
    //.filter(~/.*ultraplex.*/) // get rid of dummy meta.id
    /* .map { it -> ["id":it.toString().replaceAll(/.*ultraplex_demux_/,"").replaceAll(/\.fastq\.gz/,""),"fastq":it] }
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
        Channel.fromPath(params.smrna_genome)
    ) */

/* //unmapped reads to STAR GENOME
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
/**
//PARACLU
    PARACLU_PARACLU (
        GET_CROSSLINKS.out.crosslinkBed
    )

    PARACLU_CUT (
        PARACLU_PARACLU.out.sigxls
    )

    PARACLU_CONVERT (
        PARACLU_CUT.out.peaks
    )

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
    CLIPPY (
        GET_CROSSLINKS.out.crosslinkBed,
        file(params.gtf),
        file(params.genome_fai)
    ) */

/**
 * Post-peak calling analysis *
 */

//PEKA

}



def create_fastq_channel(LinkedHashMap row) {
    /** Takes a row from a samples CSV file and creates a meta object which
        describes it.
    */

    def meta = [:]
    meta.id           = row.entrySet().iterator().next().getValue() // This is janky and means sample id always has to come 1st
    meta.genome       = row.Species
    meta.barcode      = row.FivePrimeBarcode
    meta.single_end   = true
    return meta
}