# Faidx

`faidx.nf` creates a FAIDX genome index.

## Inputs

Required files are:

- `fasta` - a FASTA file, typically a genome FASTA file split into chromosomes.

## Processes

### `SAMTOOLS_FAIDX`

A FAIDX genome index is relatively simple. It is a plain text file, with one row per sequence in the original FASTA file.
Each sequence will give the name of the sequence, its offset from the start of the file, and its length.

It is useful in efficiently accessing a section of a large FASTA file in a computationally efficient way.

## Outputs

A `.fai` file will be produced.