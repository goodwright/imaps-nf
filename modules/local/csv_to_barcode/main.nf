
// Authors: Alex Harston, Charlotte Capitanchik

// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CSV_TO_BARCODE {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'ultraplex', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::pandas=1.1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pandas:1.1.5"
    } else {
        container "quay.io/biocontainers/pandas:1.1.5"
    }

    cache false

    input:
    path annotation

    output:
    path "barcode.csv"     , emit: csv

    script: // This script is bundled with the pipeline, in nf-core/cutandrun/bin/
    """
    python ../../../bin/csv_to_barcode.py --annotation $annotation
    """
}