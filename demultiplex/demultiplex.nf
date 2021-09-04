#!/usr/bin/env nextflow

params.fastq = ""
params.annotation = ""

fastq = file(params.fastq)
annotation = file(params.annotation)

process createBarcodeCSV {

    input:
    file annotation

    output:
    file 'barcodes.csv' into barcodes

    """
    #!/usr/bin/env python3

    import pandas as pd
    import re

    dfs = pd.read_excel("$annotation", sheet_name=None)
    matrix = list(dfs.values())[0].values

    data = [",".join([
        f'{row[11].replace("_0", "").replace(",", "")}',
    ]) for row in matrix]

    with open("barcodes.csv", "w") as f:
        f.write("\\n".join(data))
    """
}

process demultiplex {
    publishDir "."

    input:
    file barcodes from barcodes
    file fastq from fastq

    output:
    file '*.fastq.gz' into outputs

    """
    ultraplex -i $fastq -b $barcodes
    COUNT=\$(ls -lR *.fastq.gz | wc -l | xargs)
    echo "output_count: \$COUNT"
    """
}