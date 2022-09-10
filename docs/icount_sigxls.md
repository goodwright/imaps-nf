# iCount Sigxls

The `icount_sigxls.nf` pipeline determines the significant crosslinks from a BED file.

## Inputs

Required files are:

- `crosslinks` - a BED file describing the location and quantity of crosslinks.
- `segmentation` - a segemntation GTF file created by iCount segment.

## Processes

### `ICOUNT_SIGXLS`

iCount Sigxls is a peak caller - it takes the raw list of genome locations in a BED file and determines which of these are statistically significant, using parmutation analysis.

## Outputs

A compressed TSV file containing the significant crosslinks is produced, as well as accompanying scores in a separate TSV file.