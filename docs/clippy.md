# Clippy

The `clippy.nf` pipeline runs the Clippy CLIP-specific Peak Caller.

## Inputs

Required files are:

- `crosslinks` - a BED file describing the location and quantity of crosslinks.
- `gtf` - an annotation for the original genome.
- `fai` - a Faidx index for the original genome.

## Processes

### `CLIPPY`

Clippy is a peak caller - it takes a BED file identifying individual genome locations with their counts (transcription stop sites in this case, as it is a CLIP specific tool) and creates a smoothed set of peaks, using scipy's "find_peaks".

## Outputs

TSV files are produced for the peaks (nucleotide regions) and summits (individual nucleotides).

