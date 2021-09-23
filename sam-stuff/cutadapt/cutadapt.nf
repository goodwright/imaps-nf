#!/usr/bin/env nextflow

params.fastq = ""

fastq = file(params.fastq)

process cutadapt {
    publishDir "."

    input:
    file fastq from fastq

    output:
    file '*.fastq.gz' into outputfastq
    file '*.html' into html
    file '*.zip' into zip

    """
    filename=\$(basename "$params.fastq")
    cutadapt -o "cutadapt-test.fastq.gz" "$params.fastq"


    VERSION=`fastqc -v`
    printf "fastqc_version: \${VERSION}\\n"
    fastqc "cutadapt-test.fastq.gz" 2>&1
    """
}

outputfastq.collectFile()
html.collectFile()
zip.collectFile()