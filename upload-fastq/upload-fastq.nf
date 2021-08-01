#!/usr/bin/env nextflow

params.src = "*.fastq"
params.merge_lanes = true

src = file(params.src)

process uploadFastq {
    publishDir "."

    input:
    file input from src

    output:
    file '*.html' into html
    file '*.zip' into zip


    """
    VERSION=`fastqc -v`
    printf "fastqc_version: \${VERSION}\\n"
    fastqc $params.src 2>&1
    """
}

html.collectFile()
zip.collectFile()