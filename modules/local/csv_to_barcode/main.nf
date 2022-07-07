process CSV_TO_BARCODE {
    tag "$annotation"
    label "process_low"

    container "quay.io/biocontainers/pandas:1.1.5"

    input:
    path annotation

    output:
    path "barcode.csv", emit: csv

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