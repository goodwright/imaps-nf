# Paraclu

The `paraclu.nf` pipeline determines the significant crosslinks from a BED file.

## Inputs

Required files are:

- `crosslinks` - a BED file describing the location and quantity of crosslinks.

## Processes

### `PARACLU`

Paraclu is a peak caller - it takes the raw list of genome locations in a BED file and determines which of these are statistically significant.

Paraclu can identify different levels of clustering - clusters within clusters.

## Outputs

A compressed TSV file containing the significant crosslinks is produced.