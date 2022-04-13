#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { TRIMGALORE } from '../modules/nf-core/modules/trimgalore/main'  addParams( options: [:] )
include { BOWTIE_ALIGN } from '../modules/nf-core/modules/bowtie/align/main'    addParams( save_unaligned: true, options: [args:"-v 2 -m 1 --norc --best --strata"] )
include { STAR_ALIGN } from '../modules/nf-core/modules/star/align/main'    addParams( options: [args:"--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM"] )

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
    NCRNA_ANALYSIS (
        reads,
        Channel.from(params[genomeParamName])
    )
}

workflow NCRNA_ANALYSIS {

    take:
        reads
        genome

    main:

    bowtie_index = genome.map{ folder -> file(folder + "/BOWTIE_BUILD/bowtie")}
    star_index = genome.map{ folder -> file(folder + "/STAR_GENOMEGENERATE/star")}
    genome_gtf = genome.map{ folder -> file(folder + "/FILTER_GTF/*.gtf")}
    genome_fai = genome.map{ folder -> file(folder + "/SAMTOOLS_FAIDX/*.fa.fai")}
    longest_transcript = genome.map{ folder -> file(folder + "/LONGEST_TRANSCRIPT/*.txt")}
    longest_transcript_index = genome.map{ folder -> file(folder + "/LONGEST_TRANSCRIPT/*.fa.fai")}
    segmentation_gtf = genome.map{ folder -> file(folder + "/RAW_ICOUNT_SEGMENT/*segmentation*")}
    regions_gtf = genome.map{ folder -> file(folder + "/RAW_ICOUNT_SEGMENT/*regions*")}

    // Start things off with TrimGalore
    TRIMGALORE ( reads )

    // Run Bowtie Align on each of the reads, with the provided genome index
    // files as reference
    BOWTIE_ALIGN (
        TRIMGALORE.out.reads,
        bowtie_index
    )

    // Run STAR Align on the reads which didn't match above
    STAR_ALIGN (
        BOWTIE_ALIGN.out.fastq,
        star_index,
        genome_gtf
    )    
}