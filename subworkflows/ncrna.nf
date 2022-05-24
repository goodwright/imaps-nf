#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { TRIMGALORE } from '../modules/nf-core/modules/trimgalore/main'
include { BOWTIE_ALIGN } from '../modules/nf-core/modules/bowtie/align/main'
include { STAR_ALIGN } from '../modules/nf-core/modules/star/align/main'

workflow {
    // If running straight from command line, will need to construct the
    // [meta, reads, genome_name] triplet channel first

    genomeParamName = params.keySet().find{k -> k.endsWith("_genome")}
    species = genomeParamName.substring(0, 2)
    reads = [[
        id: params.fastq.split("/")[-1].replace(".gz", "").replace(".fastq", "").replace("ultraplex_demux_", ""),
        species: species,
        single_end: true
    ], file(params.fastq), genomeParamName]

    // Now just pass that along
    NCRNA_ANALYSIS (
        Channel.from([reads])
    )
}

workflow NCRNA_ANALYSIS {

    take:
        // The input channel brings in tuples of
        // [reads_meta, reads_file, genome_location]
        reads

    main:

    // Create channel just for reads and its meta data
    reads.map{triplet -> [triplet[0], triplet[1]]}.set{ ch_reads }

    // Trim the adapters of everything that comes out of ch_reads
    TRIMGALORE ( ch_reads )
}