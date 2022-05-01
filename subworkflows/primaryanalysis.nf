#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { TRIMGALORE } from '../modules/nf-core/modules/trimgalore/main'
include { BOWTIE_ALIGN } from '../modules/nf-core/modules/bowtie/align/main'
include { STAR_ALIGN } from '../modules/nf-core/modules/star/align/main'
include { UMITOOLS_DEDUP } from '../modules/nf-core/modules/umitools/dedup/main'
include { SAMTOOLS_INDEX as STAR_SAMTOOLS_INDEX} from '../modules/nf-core/modules/samtools/index/main'
include { SAMTOOLS_INDEX as UMITOOLS_SAMTOOLS_INDEX} from '../modules/nf-core/modules/samtools/index/main'
include { GET_CROSSLINKS } from '../modules/local/get_crosslinks/main'
include { CROSSLINKS_COVERAGE } from '../modules/luslab/nf-core-modules/crosslinks/coverage/main'
include { CROSSLINKS_NORMCOVERAGE } from '../modules/luslab/nf-core-modules/crosslinks/normcoverage/main'

include { FILTER_TRANSCRIPTS } from '../modules/local/filter_transcriptome_bam/main'
include { UMITOOLS_DEDUP as TOME_UMITOOLS_DEDUP } from '../modules/nf-core/modules/umitools/dedup/main'
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

workflow {
    // If running straight from command line, will need to construct the
    // [meta, reads] pair channel first
    reads = [[
        id: params.fastq.split("/")[-1].replace(".gz", "").replace(".fastq", "").replace("ultraplex_demux_", ""),
        single_end: true
    ], file(params.fastq)]

    // What is the genome param called?
    genomeParamName = params.keySet().find{k -> k.endsWith("_genome")}

    // Now just pass that along with the rest of params
    PRIMARY_ANALYSIS (
        reads,
        Channel.from(params[genomeParamName])
    )
}

workflow PRIMARY_ANALYSIS {

    take:
        reads
        genome

    main:

    fasta = genome.map{ folder -> (
        file(folder + "/DNA_GUNZIP/*.fa") ?
        file(folder + "/DNA_GUNZIP/*.fa") : file(folder + "/inputs/fasta/*.fa")
    ) }
    bowtie_index = genome.map{ folder -> file(folder + "/BOWTIE_BUILD/bowtie")}
    star_index = genome.map{ folder -> file(folder + "/STAR_GENOMEGENERATE/star")}
    genome_gtf = genome.map{ folder -> file(folder + "/FILTER_GTF/*.gtf")}
    genome_fai = genome.map{ folder -> file(folder + "/SAMTOOLS_FAIDX/*.fa.fai")}
    longest_transcript = genome.map{ folder -> file(folder + "/FIND_LONGEST_TRANSCRIPT/*.txt")}
    longest_transcript_index = genome.map{ folder -> file(folder + "/FIND_LONGEST_TRANSCRIPT/*.fa.fai")}
    segmentation_gtf = genome.map{ folder -> file(folder + "/RAW_ICOUNT_SEGMENT/*segmentation*")}
    regions_gtf = genome.map{ folder -> file(folder + "/RESOLVE_UNANNOTATED/*regions*")}
    genic_other_regions_gtf = genome.map{ folder -> file(folder + "/RESOLVE_UNANNOTATED_GENIC_OTHER/*regions*")}

    // Start things off with TrimGalore
    TRIMGALORE ( reads )

    // Run Bowtie Align on each of the reads, with the provided genome index
    // files as reference
    BOWTIE_ALIGN ( TRIMGALORE.out.reads, bowtie_index )

    // Run STAR Align on the reads which didn't match above
    STAR_ALIGN ( BOWTIE_ALIGN.out.fastq, star_index, genome_gtf, false, "", "" )


    // Get TOME crosslinks
    FILTER_TRANSCRIPTS ( STAR_ALIGN.out.bam_transcript, longest_transcript )

    TOME_STAR_SAMTOOLS_INDEX ( FILTER_TRANSCRIPTS.out.filtered_bam )
    tome_ch_umi_input = FILTER_TRANSCRIPTS.out.filtered_bam.combine(TOME_STAR_SAMTOOLS_INDEX.out.bai, by: 0)
    TOME_UMITOOLS_DEDUP ( tome_ch_umi_input )

    TOME_UMITOOLS_SAMTOOLS_INDEX ( TOME_UMITOOLS_DEDUP.out.bam )
    tome_ch_xl_input = TOME_UMITOOLS_DEDUP.out.bam.combine(TOME_UMITOOLS_SAMTOOLS_INDEX.out.bai, by: 0)
    TOME_GET_CROSSLINKS ( tome_ch_xl_input, longest_transcript_index )

    TOME_CROSSLINKS_COVERAGE ( TOME_GET_CROSSLINKS.out.crosslinkBed )
    TOME_CROSSLINKS_NORMCOVERAGE ( TOME_GET_CROSSLINKS.out.crosslinkBed )


    // Get crosslinks
    STAR_SAMTOOLS_INDEX ( STAR_ALIGN.out.bam_sorted )
    ch_umi_input = STAR_ALIGN.out.bam_sorted.combine(STAR_SAMTOOLS_INDEX.out.bai, by: 0)
    UMITOOLS_DEDUP ( ch_umi_input )

    UMITOOLS_SAMTOOLS_INDEX ( UMITOOLS_DEDUP.out.bam )
    ch_xl_input = UMITOOLS_DEDUP.out.bam.combine(UMITOOLS_SAMTOOLS_INDEX.out.bai, by: 0)
    GET_CROSSLINKS ( ch_xl_input, genome_fai )

    CROSSLINKS_COVERAGE ( GET_CROSSLINKS.out.crosslinkBed )
    CROSSLINKS_NORMCOVERAGE ( GET_CROSSLINKS.out.crosslinkBed )


    // CLIPPY Peak Calling
    CLIPPY ( GET_CROSSLINKS.out.crosslinkBed, genome_gtf, genome_fai )


    // Paraclu Peak Calling
    PARACLU_PARACLU ( GET_CROSSLINKS.out.crosslinkBed )
    PARACLU_CUT ( PARACLU_PARACLU.out.sigxls )
    PARACLU_CONVERT ( PARACLU_CUT.out.peaks )


    // iCount
    ICOUNT_SUMMARY ( GET_CROSSLINKS.out.crosslinkBed, regions_gtf )
    ICOUNT_RNAMAPS ( GET_CROSSLINKS.out.crosslinkBed, regions_gtf )
    ICOUNT_SIGXLS ( GET_CROSSLINKS.out.crosslinkBed, segmentation_gtf )

    ch_icount_peaks = GET_CROSSLINKS.out.crosslinkBed.combine(ICOUNT_SIGXLS.out.sigxls, by: 0)
    ICOUNT_PEAKS ( ch_icount_peaks )

    // PEKA
    PEKA (
        CLIPPY.out.peaks,
        GET_CROSSLINKS.out.crosslinkBed,
        fasta,
        genome_fai,
        genic_other_regions_gtf
    )
}