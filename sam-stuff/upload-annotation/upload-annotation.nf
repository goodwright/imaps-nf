#!/usr/bin/env nextflow

params.xlsx = "*.xlsx"

xlsx = file(params.xlsx)

process uploadAnnotation {
    publishDir "."

    input:
    file input from xlsx


    """
    #!/usr/bin/env python3

    import pandas as pd

    if not "$params.xlsx".endswith(".xlsx"):
        raise Exception("Uploaded file is not an .xlsx file.")
    
    try:
        dfs = pd.read_excel("$params.xlsx", sheet_name=None)
    except ValueError:
        raise Exception(".xlsx file could not be parsed.")

    # Correct columns there?
    REQUIRED = [
        "Sample name", "collection name", "experimental design", "Scientist", "PI",
        "Method", "Protocol document", "Protein", "cells/tissue", "mapto", "5' barcode",
        "3' adapter", "sequencer type", "antibody used", "RT primer", "linker"
    ]
    keys = list(dfs.values())[0].keys()
    for key in REQUIRED:
        if key not in keys:
            raise Exception(".xlsx does not have a '{}' column".format(key))
    """
}