# Find Longest Transcript

`find_longest_transcript.nf` finds the longest transcript for each gene.

## Inputs

Required files are:

- `gtf` - a GTF annotation describing a genome.

## Processes

### `FIND_LONGEST_TRANSCRIPT`

For each gene, one transcript is identified as the longest for that gene and output as a list of transcript IDs.

The total length of all CDS sequences is used to determine the length.
Where there are multiple transcripts with the same CDS length, the total exon length is used.

## Outputs

The transcript IDs are produced as a txt file, and a genome FAI index describing where the transcripts can be found is also produced.