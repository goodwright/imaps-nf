# Demultiplex

The `crosslinks_coverage.nf` pipeline 

## Inputs

Required files are:

- `crosslinks` - a BED file describing the location and quantity of crosslinks."

## Processes

### `CROSSLINKS_COVERAGE`

The BED file is converted to Bedgraph format, retaining the number of crosslinks per genome site.

## Outputs

A Bedgraph file is produced counting the number of crosslinks at each genome site.