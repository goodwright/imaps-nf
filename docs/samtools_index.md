# Samtools Index

`samtools_index.nf` creates a BAI alignment index.

## Inputs

Required files are:

- `bam` - a BAM file.

## Processes

### `SAMTOOLS_INDEX`

BAI files are indexes to a BAI alignment file, which allow for memory-efficient access to regions within it.
This process produces that index for a given BAM file.

## Outputs

A `.bai` file will be produced.