# Paraclu Convert

The `paraclu_convert.nf` pipeline cconverts a Paraclu TSV file to BED format

## Inputs

Required files are:

- `peaks` - a TSV output file of a Paraclu pipeline.

## Processes

### `PARACLU_CONVERT`

Paraclu produces TSV files as outputs, but typically genome peaks should be processed as BED files.
Paraclu Convert turns a Paraclu TSV file into a BED file.

## Outputs

A compressed BED file is produced.