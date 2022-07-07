# Demultiplex

The `demultiplex.nf` pipeline demultiplexes a multiplexed FASTQ reads file using barcodes in a provided spreadsheet.

## Inputs

Required files are:

- `multiplexed_fastq` - a FASTQ file containing reads from multiple samples, each with a barcode sequence at the start.
- `annotation` - a CSV or XLSX sheet of data, which needs at least a "Sample Name" column and a "5' Barcode" column.

## Processes

### `XLSX_TO_CSV`

If the `annotation` provided is an Excel spreadsheet, this process wil convert it to a CSV file.

### `CSV_TO_BARCODE`

Converts a full annotation spreadsheet to one containing just the columns that Ultraplex needs, mostly the barcode information.

### `ULTRAPLEX`

Runs the Ultraplex demultiplexer to create multiple demultiplexed reads file from the `multiplexed_fastq` file, using the generated barcodes file.
Output filenames will be the sample names from the spreadsheet with the `ultraplex_demux_` prefix added.
A log file will also be produced.

The barcode section at the start of the reads that Ultraplex looks for has 'N' bases at the start and end, and the specific barcode sequence in the middle.
The barcode sequence is used to demultiplex, the start and end sections are combined and added to the reads header as an `rbc`, and the whole barcode section is removed.

## Outputs

The main output files are the multiple `.fastq.gz` files representing demultiplexed sample files.