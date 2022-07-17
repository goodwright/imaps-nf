# Bowtie Build

`bowtie_build.nf` creates a BOWTIE genome index.

## Inputs

Required files are:

- `fasta` - a FASTA file, typically a genome FASTA file.

## Processes

### `BOWTIE_BUILD`

A BOWTIE index is a binary representation of a genome that BOWTIE can use to do read alignment.
This process creates that index from a plain FASTA file.

## Outputs

A `bowtie` directory will be produced, containing the `.ebwt` files that comprise the index.