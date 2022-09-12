# iCount Summary

The `icount_summary.nf` reports count of cross-link events in each region type.

## Inputs

Required files are:

- `crosslinks` - a BED file describing the location and quantity of crosslinks.
- `regions` - a regions GTF file created by iCount segment.

## Processes

### `ICOUNT_SUMMARY`

In addition to peak calling, iCount summary can also summarise the overall information content of a crosslinks file.

## Outputs

Returns metrics for genes and subtypes in the form of TSV files.