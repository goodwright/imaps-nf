
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

    import pandas as pd

    data = pd.read_csv("$annotation")
    data = data[["5' Barcode", 'Sample Name']]
    new = [f"{x}:{y}" for x,y in zip(data["5' Barcode"], data["Sample Name"])]
    newdf = pd.DataFrame(new)
    newdf.to_csv("barcode.csv", index=False, header=False)
    """
}