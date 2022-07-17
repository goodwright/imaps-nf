# STAR Build

`star_build.nf` creates a STAR genome index.

## Inputs

Required files are:

- `fasta` - a FASTA file, typically a genome FASTA file.
- `gtf` - a GTF annotation file.

## Processes

### `STAR_GENOMEGENERATE`

A STAR index is a binary representation of a genome that STAR can use to do read alignment.
This process creates that index from a plain FASTA file.

The GTF file is needed because it contains information on splice junctions, which the binary index incorporates.

## Outputs

A `star` directory will be produced, containing the various files that comprise the index.