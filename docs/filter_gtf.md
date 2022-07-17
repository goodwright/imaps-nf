# Filter GTF

`filter_gtf.nf` trims a GTF file down to high-quality sequences.

## Inputs

Required files are:

- `gtf` - a GTF annotation describing a genome.

## Processes

### `FILTER_GTF`

Certain transcripts will be removed from the GTF.

If the GTF file contains any transcripts marked as 'basic', all transcripts not marked as basic are removed.
The 'basic' tag is used to mark certain transcripts as representative when a gene has more than one transcript.

If the GTF uses the 'transcript support level' scoring system, those transcripts that aren't level 1 or 2 are removed.
These are considered to not have sufficient experimental evidence.

## Outputs

A 'post-filtering' GTF is produced.