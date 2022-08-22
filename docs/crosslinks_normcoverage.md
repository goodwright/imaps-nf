# Crosslinks Normalised Coverage

The `crosslinks_normcoverage.nf` pipeline creates a normalised Bedgraph file from a crosslinks alignment file.

## Inputs

Required files are:

- `crosslinks` - a BED file describing the location and quantity of crosslinks.

## Processes

### `CROSSLINKS_NORMCOVERAGE`

The BED file is converted to Bedgraph format, retaining the number of crosslinks per genome site normalised against each other.

## Outputs

A Bedgraph file is produced counting the number of crosslinks at each genome site.