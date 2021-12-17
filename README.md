# imaps-nf

## Modules

### `CSV_TO_BARCODE`

This module takes a CSV file describing a multiplexed reads file, and produces a
CSV file containing just the barcode information. That is, for each sample it
extracts the sample name and the barcode sequence that identifies it, formatted
for the needs of Ultraplex. It uses pandas to do this.
