# Filter Transcriptome BAM

`filter_transcriptome_bam.nf` trims a BAM file down based on a list of transcript IDs.

## Inputs

Required files are:

- `bam` - a BAM file with reads aligned to a transcriptome.
- `transcripts` - a list of transcript IDs to preserve.

## Processes

### `FILTER_TRANSCRIPTOME_BAM`

A BAM file lists a series of reads aligned to a genome, and this process will remove those which aren't to to a region of the genome within certain transcripts.
These transcripts are provided as a list of transcript IDs.

## Outputs

A 'post-filtering' BAM file is produced.