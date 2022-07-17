# FASTQC

`fastqc.nf` performs quality checks on a reads file.

## Inputs

Required files are:

- `fastq` - a FASTQ reads file to be assessed.

## Processes

### `FASTQC`

FASTQC performs quality control on a reads file, using the information in the base calling confidence line to determine how reliable overall the reads file is.

Some metrics calculated are average quality per read, distribution of quality scores, and level of overrepresented sequences.

## Outputs

FASTQC produces a HTML report, which shows the stats produced visually and in a human readable way.
It also produces a zip file of machine readable text files containing the same information.