# Bowtie Align

`bowtie_align.nf` aligns FASTQ reads files against a Bowtie genome index.

## Inputs

Required files are:

- `fastq` - a FASTQ reads file to be aligned.
- `index` - a Bowtie genome index to align against.

## Processes

### `BOWTIE_ALIGN`

Bowtie is an ultrafast, memory-efficient short read aligner.
It aligns short DNA sequences (reads) to the human genome at a rate of over 25 million 35-bp reads per hour.
Bowtie indexes the genome with a Burrows-Wheeler index to keep its memory footprint small: typically about 2.2 GB for the human genome (2.9 GB for paired-end).

## Outputs

The aligned reads are output as a BAM file, while unaligned reads are output as a new reads FASTQ file.