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
    file 'output.txt' into output


    """
    VERSION=`fastqc -v`
    printf "fastqc_version: \${VERSION}\n" >> output.txt 2>&1
    fastqc $params.src >> output.txt 2>&1
    """
}

html.collectFile()
zip.collectFile()
output.collectFile()