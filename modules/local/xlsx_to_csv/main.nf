process XLSX_TO_CSV {
    tag "$xlsx"
    label "process_low"

    container "quay.io/biocontainers/pandas:1.1.5"

    input:
    path xlsx

    output:
    path "*.csv", emit: csv

    script:
    """
    pip install openpyxl
    python -c "import pandas as pd; data = pd.read_excel('$xlsx', engine='openpyxl'); data.to_csv('$xlsx' + '.csv', index=False)"
    """
}