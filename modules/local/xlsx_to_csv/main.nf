
// Authors: Sam Ireland
include { initOptions; saveFiles } from './functions'

params.options = [:]
options        = initOptions(params.options)

process XLSX_TO_CSV {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'xlsx_to_csv', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::pandas=1.1.5" : null)
    container "quay.io/biocontainers/pandas:1.1.5"

    cache false

    input:
    path xlsx

    output:
    path "*.csv"     , emit: csv

    script:
    """
    pip install openpyxl

    python -c "import pandas as pd; data = pd.read_excel('$xlsx', engine='openpyxl'); data.to_csv('$xlsx' + '.csv', index=False)"

    """
}