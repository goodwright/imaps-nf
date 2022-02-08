
// Authors: Alex Harston, Charlotte Capitanchik

// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CSV_TO_BARCODE {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'csv_to_barcode', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::pandas=1.1.5" : null)
    container "quay.io/biocontainers/pandas:1.1.5"

    cache false

    input:
    path annotation

    output:
    path "barcode.csv"     , emit: csv

    script:
    """
    #!/usr/bin/env python3

    from pandas import read_csv
    from sys import exit

    data = read_csv("$annotation", dtype=str, keep_default_na=False)

    five_prime = data["5' Barcode"]
    three_prime = data["3' Barcode (optional)"]
    sample_names = data["Sample Name"]

    barcode_dict = {}

    for idx in range(len(five_prime)):
        barcode_dict.setdefault(five_prime[idx], [])
        barcode_dict[five_prime[idx]].append(
            three_prime[idx] + ":" + sample_names[idx]
        )

    with open("barcode.csv", "w") as out_f:
        for five, threes in barcode_dict.items():
            if len(threes) > 1:
                if any([three.startswith(":") for three in threes]):
                    exit("5' barcode ambiguity between samples")
                out_f.write(",".join([five] + threes) + "\\n")
            else:
                if not threes[0].startswith(":"):
                    threes[0] = "," + threes[0]
                out_f.write(five + threes[0] + "\\n")
    """
}
