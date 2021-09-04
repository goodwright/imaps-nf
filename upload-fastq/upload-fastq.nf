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
    fastqc "$params.src" 2>&1

    zipname=\$(ls *.zip)
    filename="\${zipname%.*}"
    unzip -q "\$zipname"

    line=\$(head -n 1 "\$filename/summary.txt")
    declare -a linearr=(`echo \$line`)
    pass_or_fail=\$(echo \${linearr[0]})
    if [[ "\$pass_or_fail" == "PASS" ]]; then
        echo "quality_pass: true"
    else
        echo "quality_pass: false"
    fi
    """
}

html.collectFile()
zip.collectFile()