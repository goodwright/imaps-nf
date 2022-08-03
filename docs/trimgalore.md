# TrimGalore

`trimgalore.nf` removes sequencing artifacts from reads files.

## Inputs

Required files are:

- `fastq` - a FASTQ reads file to be trimmed.

## Processes

### `TRIMGALORE`

Each read in a reads file represents some fragment of a DNA sequence sent for sequencing, but it will also contain bits of DNA added to it as part of the sequencing process, such as 3; sequencing adaptors.
TrimGalore removes those to leave only the biological DNA sequence of interest.

Sequences shorter than 10 bases will be removed as well.

## Outputs

A 'trimmed' FASTQC file will be produced, as well as a trimming report text file.