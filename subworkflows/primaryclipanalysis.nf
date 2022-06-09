#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { TRIMGALORE } from '../modules/nf-core/modules/trimgalore/main'
include { BOWTIE_ALIGN } from '../modules/nf-core/modules/bowtie/align/main'
include { STAR_ALIGN } from '../modules/nf-core/modules/star/align/main'
include { DU } from '../modules/local/du/main'
include { GET_UMI_LENGTH } from '../modules/local/get_umi_length/main'
include { UMITOOLS_DEDUP } from '../modules/local/umitools_dedup/main'
include { SAMTOOLS_INDEX as STAR_SAMTOOLS_INDEX} from '../modules/nf-core/modules/samtools/index/main'
include { SAMTOOLS_INDEX as UMITOOLS_SAMTOOLS_INDEX} from '../modules/nf-core/modules/samtools/index/main'
include { GET_CROSSLINKS } from '../modules/local/get_crosslinks/main'
include { CROSSLINKS_COVERAGE } from '../modules/luslab/nf-core-modules/crosslinks/coverage/main'
include { CROSSLINKS_NORMCOVERAGE } from '../modules/luslab/nf-core-modules/crosslinks/normcoverage/main'

include { FILTER_TRANSCRIPTS } from '../modules/local/filter_transcriptome_bam/main'
include { DU as TOME_DU } from '../modules/local/du/main'
include { GET_UMI_LENGTH as TOME_GET_UMI_LENGTH } from '../modules/local/get_umi_length/main'
include { UMITOOLS_DEDUP as TOME_UMITOOLS_DEDUP } from '../modules/local/umitools_dedup/main'
include { SAMTOOLS_INDEX as TOME_STAR_SAMTOOLS_INDEX } from '../modules/nf-core/modules/samtools/index/main'
include { SAMTOOLS_INDEX as TOME_UMITOOLS_SAMTOOLS_INDEX } from '../modules/nf-core/modules/samtools/index/main'
include { GET_CROSSLINKS as TOME_GET_CROSSLINKS } from '../modules/local/get_crosslinks/main'
include { CROSSLINKS_COVERAGE as TOME_CROSSLINKS_COVERAGE } from '../modules/luslab/nf-core-modules/crosslinks/coverage/main'
include { CROSSLINKS_NORMCOVERAGE as TOME_CROSSLINKS_NORMCOVERAGE } from '../modules/luslab/nf-core-modules/crosslinks/normcoverage/main'

include { ICOUNT_SIGXLS } from '../modules/luslab/nf-core-modules/icount/sigxls/main'
include { ICOUNT_SUMMARY } from '../modules/local/icount_summary/main'
include { ICOUNT_RNAMAPS } from '../modules/local/icount_rnamaps/main'
include { ICOUNT_PEAKS } from '../modules/local/icount_peaks/main'
include { CLIPPY } from '../modules/luslab/nf-core-modules/clippy/main'
include { PARACLU_CONVERT } from '../modules/luslab/nf-core-modules/paraclu/convert/main'
include { PARACLU_PARACLU } from '../modules/luslab/nf-core-modules/paraclu/paraclu/main'
include { PARACLU_CUT } from '../modules/luslab/nf-core-modules/paraclu/cut/main'
include { PEKA } from '../modules/luslab/nf-core-modules/peka/main'

include { CLIP_QC } from '../modules/local/clip_qc/main.nf'

// Closure to annotate UMITools Input
annotate_umitools_input = { it ->
    def meta = it[0].clone()
    if (it[3].toInteger() > params.max_kilobytes &
        it[4].toInteger() > params.max_umi_length) {
        meta["low_memory"] = true
    }
    return [meta, it[1], it[2]]
}

workflow {
    // If running straight from command line, will need to construct the
    // [meta, reads, genome_name] triplet channel first

    genomeParamName = params.keySet().find{k -> k.endsWith("_genome")}
    species = genomeParamName.substring(0, 2)
    reads = [[
        id: params.fastq.split("/")[-1].replace(".gz", "").replace(".fastq", "").replace("ultraplex_demux_", ""),
        species: species,
        single_end: true
    ], file(params.fastq), params[genomeParamName]]

    // Now just pass that along
    PRIMARY_CLIP_ANALYSIS (
        Channel.from([reads])
    )
}

workflow PRIMARY_CLIP_ANALYSIS {

    take:
        // The input channel brings in tuples of
        // [reads_meta, reads_file, genome_location]
        reads

    main:

    // Create channel just for reads and its meta data
    reads.map{triplet -> [triplet[0], triplet[1]]}.set{ ch_reads }

    // Trim the adapters of everything that comes out of ch_reads
    TRIMGALORE ( ch_reads )

    // Create a channel which outputs [reads_meta, bowtie_index] pairs
    reads.map{triplet -> [
        triplet[0], file(triplet[2] + "/BOWTIE_BUILD/bowtie")
    ]}.set{ ch_bowtie }

    // Create a channel which creates
    // [reads_meta, trimmed_reads, bowtie_index] triplets
    TRIMGALORE.out.reads.join(ch_bowtie).set { ch_reads_with_bowtie_index }

    // Multi-map channel into named outputs
    ch_reads_with_bowtie_index.multiMap { triplet ->
        trimgalore: [triplet[0], triplet[1]]
        bowtie: triplet[2]
    }.set { ch_bowtie_input }

    // Align reads to associated genome bowtie index
    BOWTIE_ALIGN ( ch_bowtie_input.trimgalore, ch_bowtie_input.bowtie )

    // Create a channel which outputs [reads_meta, star_index] pairs
    reads.map{triplet -> [
        triplet[0], file(triplet[2] + "/STAR_GENOMEGENERATE/star")
    ]}.set{ ch_star }

    // Create a channel which outputs [reads_meta, star_index] pairs
    reads.map{triplet -> [
        triplet[0], file(triplet[2] + "/FILTER_GTF/*.gtf")
    ]}.set{ ch_gtf }

    // Create a channel which creates
    // [reads_meta, bowtie_output, star_index, gtf] tuples
    BOWTIE_ALIGN.out.fastq.join( ch_star ).join( ch_gtf ).set{ ch_reads_with_star_index }

    // Multi-map channel into named outputs
    ch_reads_with_star_index.multiMap { tuple ->
        bowtie: [tuple[0], tuple[1]]
        star: tuple[2]
        gtf: tuple[3]
    }.set { ch_star_input }

    // Align reads to associated genome STAR index
    STAR_ALIGN (
        ch_star_input.bowtie, ch_star_input.star, ch_star_input.gtf,
        false, "", ""
    )

    // Create a channel which outputs [reads_meta, transcript_txt] pairs
    reads.map{triplet -> [
        triplet[0], file(triplet[2] + "/FIND_LONGEST_TRANSCRIPT/*.txt")
    ]}.set{ ch_longest_transcript }

    // Create a channel which creates
    // [reads_meta, star_transcripts, transcript_txt] triplets
    STAR_ALIGN.out.bam_transcript
    .join( ch_longest_transcript )
    .set{ ch_star_transcripts_with_longest }

    // Multi-map channel into named outputs
    ch_star_transcripts_with_longest.multiMap { triplet ->
        star: [triplet[0], triplet[1]]
        transcripts: triplet[2]
    }.set { ch_filter_input }

    // Filter transcripts
    FILTER_TRANSCRIPTS ( ch_filter_input.star, ch_filter_input.transcripts )

    // Get TOME crosslinks
    TOME_STAR_SAMTOOLS_INDEX ( FILTER_TRANSCRIPTS.out.filtered_bam )
    FILTER_TRANSCRIPTS.out.filtered_bam.join(TOME_STAR_SAMTOOLS_INDEX.out.bai)
        .set{ tome_ch_umi_input }

    // Determine if UMITools needs to be run in "low_memory" mode
    TOME_DU ( tome_ch_umi_input.map{it -> it[0, 1]} )
    TOME_GET_UMI_LENGTH ( tome_ch_umi_input )
    tome_ch_umi_input
        .join( TOME_DU.out.size )
        .join( TOME_GET_UMI_LENGTH.out.length )
        .map( annotate_umitools_input )
        .set{ tome_ch_umi_input_annotated }

    TOME_UMITOOLS_DEDUP ( tome_ch_umi_input_annotated )
    TOME_UMITOOLS_DEDUP.out.bam
        .map{ it -> [it[0].findAll{key, val -> key != "low_memory"}, it[1]] }
        .set{ ch_tome_umitools_bam }
    TOME_UMITOOLS_SAMTOOLS_INDEX ( ch_tome_umitools_bam )
    reads.map{triplet -> [
        triplet[0], file(triplet[2] + "/FIND_LONGEST_TRANSCRIPT/*.fa.fai")
    ]}.set{ ch_longest_transcript_index }
    tome_ch_xl_input = ch_tome_umitools_bam.join(TOME_UMITOOLS_SAMTOOLS_INDEX.out.bai)
    tome_ch_xl_input.join( ch_longest_transcript_index ).set{ tome_with_index }
    tome_with_index.multiMap { tuple ->
        bam: [tuple[0], tuple[1], tuple[2]]
        transcript: tuple[3]
    }.set { ch_tome_input }
    TOME_GET_CROSSLINKS ( ch_tome_input.bam, ch_tome_input.transcript )
    TOME_CROSSLINKS_COVERAGE ( TOME_GET_CROSSLINKS.out.crosslinkBed )
    TOME_CROSSLINKS_NORMCOVERAGE ( TOME_GET_CROSSLINKS.out.crosslinkBed )


    // Get crosslinks
    STAR_SAMTOOLS_INDEX ( STAR_ALIGN.out.bam_sorted )
    ch_umi_input = STAR_ALIGN.out.bam_sorted.combine(STAR_SAMTOOLS_INDEX.out.bai, by: 0)

    // Determine if UMITools needs to be run in "low_memory" mode
    DU ( ch_umi_input.map{it -> it[0, 1]} )
    GET_UMI_LENGTH ( ch_umi_input )
    ch_umi_input
        .join( DU.out.size )
        .join( GET_UMI_LENGTH.out.length )
        .map( annotate_umitools_input )
        .set{ ch_umi_input_annotated }

    UMITOOLS_DEDUP ( ch_umi_input_annotated )

    // Strip out the low_memory key from the meta value so that the later joins
    // actually work
    UMITOOLS_DEDUP.out.bam
        .map{ it -> [it[0].findAll{key, val -> key != "low_memory"}, it[1]] }
        .set{ ch_umitools_bam }
    // Keep a channel for converting between meta with the low_memory key and
    // without, in case in the future you want to keep track of which files were
    // run as low mem
    UMITOOLS_DEDUP.out.bam
        .map{ it -> [it[0].findAll{key, val -> key != "low_memory"}, it[0]] }
        .set{ ch_meta_conversion }

    UMITOOLS_SAMTOOLS_INDEX ( ch_umitools_bam )

    reads.map{triplet -> [
        triplet[0], file(triplet[2] + "/SAMTOOLS_FAIDX/*.fa.fai")
    ]}.set{ ch_genome_fai }

    ch_xl_input = ch_umitools_bam.join(UMITOOLS_SAMTOOLS_INDEX.out.bai)

    ch_xl_input.join( ch_genome_fai ).set{ ch_with_index }

    ch_with_index.multiMap { tuple ->
        bam: [tuple[0], tuple[1], tuple[2]]
        fai: tuple[3]
    }.set { ch_xl_input }
    GET_CROSSLINKS ( ch_xl_input.bam, ch_xl_input.fai )
    CROSSLINKS_COVERAGE ( GET_CROSSLINKS.out.crosslinkBed )
    CROSSLINKS_NORMCOVERAGE ( GET_CROSSLINKS.out.crosslinkBed )

    // CLIPPY Peak Calling
    GET_CROSSLINKS.out.crosslinkBed.join( ch_gtf ).join( ch_genome_fai ).set{ ch_crosslinks }
    ch_crosslinks.multiMap { tuple ->
        crosslinks: [tuple[0], tuple[1]]
        gtf: tuple[2]
        fai: tuple[3]
    }.set { ch_clippy_input }
    CLIPPY ( ch_clippy_input.crosslinks,  ch_clippy_input.gtf,  ch_clippy_input.fai )

    // Paraclu Peak Calling
    PARACLU_PARACLU ( GET_CROSSLINKS.out.crosslinkBed )
    PARACLU_CUT ( PARACLU_PARACLU.out.sigxls )
    PARACLU_CONVERT ( PARACLU_CUT.out.peaks )

    // iCount
    reads.map{triplet -> [
        triplet[0], file(triplet[2] + "/RAW_ICOUNT_SEGMENT/*segmentation*")
    ]}.set{ ch_segmentation }
    reads.map{triplet -> [
        triplet[0], file(triplet[2] + "/RESOLVE_UNANNOTATED/*regions*")
    ]}.set{ ch_regions }
    reads.map{triplet -> [
        triplet[0], file(triplet[2] + "/RESOLVE_UNANNOTATED_GENIC_OTHER/*regions*")
    ]}.set{ ch_regions_genic_other }
    GET_CROSSLINKS.out.crosslinkBed
    .join( ch_regions )
    .set{ ch_crosslinks_with_regions }
    GET_CROSSLINKS.out.crosslinkBed
    .join( ch_segmentation )
    .set{ ch_crosslinks_with_segmentation }
    ch_crosslinks_with_regions.multiMap { triplet ->
        crosslinks: [triplet[0], triplet[1]]
        regions: triplet[2]
    }.set { ch_regions_input }
    ch_crosslinks_with_segmentation.multiMap { triplet ->
        crosslinks: [triplet[0], triplet[1]]
        segmentation: triplet[2]
    }.set { ch_segmentation_input }
    ICOUNT_SUMMARY ( ch_regions_input.crosslinks, ch_regions_input.regions )
    ICOUNT_RNAMAPS ( ch_regions_input.crosslinks, ch_regions_input.regions )
    ICOUNT_SIGXLS ( ch_segmentation_input.crosslinks, ch_segmentation_input.segmentation )
    ch_icount_peaks = GET_CROSSLINKS.out.crosslinkBed.combine(ICOUNT_SIGXLS.out.sigxls, by: 0)
    ICOUNT_PEAKS ( ch_icount_peaks )

    // PEKA
    reads.map{triplet -> [
        triplet[0],
        file(triplet[2] + "/DNA_GUNZIP/*.fa") ?
        file(triplet[2] + "/DNA_GUNZIP/*.fa") : file(triplet[2] + "/inputs/fasta/*.fa")
    ]}.set{ ch_fasta }
    CLIPPY.out.peaks
    .join(GET_CROSSLINKS.out.crosslinkBed)
    .join(ch_fasta)
    .join(ch_genome_fai)
    .join(ch_regions_genic_other)
    .set{ ch_peka_joins }
    ch_peka_joins.multiMap { tuple ->
        peaks: [tuple[0], tuple[1]]
        crosslinks: [tuple[0], tuple[2]]
        fasta: tuple[3]
        fai: tuple[4]
        regions: tuple[5]
    }.set { ch_peka_input }
    PEKA (
        ch_peka_input.peaks,
        ch_peka_input.crosslinks,
        ch_peka_input.fasta,
        ch_peka_input.fai,
        ch_peka_input.regions,
    )

    CLIP_QC (
        BOWTIE_ALIGN.out.log.map{ vec -> vec[1] }.collect(),
        STAR_ALIGN.out.log_final.map{ vec -> vec[1] }.collect(),
        GET_CROSSLINKS.out.crosslinkBed.map{ vec -> vec[1] }.collect(),
        ICOUNT_PEAKS.out.peaks.map{ vec -> vec[1] }.collect(),
        PARACLU_CONVERT.out.peaks.map{ vec -> vec[1] }.collect(),
        CLIPPY.out.peaks.map{ vec -> vec[1] }.collect()
    )

    emit:
        trimgalore_log       = TRIMGALORE.out.log
        bowtie_align_log     = BOWTIE_ALIGN.out.log
        star_align_log_final = STAR_ALIGN.out.log_final
        clip_qc_log          = CLIP_QC.out.log
}