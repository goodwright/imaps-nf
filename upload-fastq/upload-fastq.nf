#!/usr/bin/env nextflow

params.src = "*.fastq"
params.merge_lanes = true

src = file(params.src)

process uploadFastq {
    publishDir "."

    input:
    file input from src

    output:
    stdout ch
    file '*.html' into html
    file '*.zip' into zip


    """
    fastqc $params.src
    """
}

html.collectFile()
zip.collectFile()