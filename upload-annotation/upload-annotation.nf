#!/usr/bin/env nextflow

params.src = "*.xlsx"
params.merge_lanes = true

src = file(params.src)

process uploadFastq {
    publishDir "."

    input:
    file input from src

    output:
    file 'output.txt' into output


    """
    ls >> output.txt 2>&1
    """
}

output.collectFile()