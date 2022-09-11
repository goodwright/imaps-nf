# iCount RnaMaps

The `icount_rnamaps.nf` computes distribution of cross-links relative to genomic landmarks.

## Inputs

Required files are:

- `crosslinks` - a BED file describing the location and quantity of crosslinks.
- `regions` - a regions GTF file created by iCount segment.

## Processes

### `ICOUNT_RNAMAPS`

iCount RnaMaps computes distribution of cross-links relative to genomic landmarks, from an initial crosslinks file.

## Outputs

A directory of useful stats and figures is produced.