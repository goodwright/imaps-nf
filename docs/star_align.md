# STAR Align

`star_align.nf` aligns FASTQ reads files against a STAR genome index.

## Inputs

Required files are:

- `fastq` - a FASTQ reads file to be aligned.
- `index` - a STAR genome index to align against.
- `gtf` - a GTF annotation describing the genome.

## Processes

### `STAR_ALIGN`

STAR aligns reads to a genome (in the form of a STAR genome index), and reports the aligned reads as a BAM file.
In this case the reads are sorted by coordinate.

STAR also translates the aligned reads into an alignment against the transcriptome only, and outputs this as a separate BAM file.

## Outputs

Two BAM files are produced - one basic alignment sorted by coordinate, and one transcriptome alignment.
It also produces various log files.