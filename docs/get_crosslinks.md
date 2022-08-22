# Get Crosslinks

`get_crosslinks.nf` identifies the regions of a genome that reads pile up at.

## Inputs

Required files are:

- `bam` - a BAM file.
- `bai` - a BAI index for the BAM.
- `fai` - a Faidx index file.

## Processes

### `GET_CROSSLINKS`

The input BAM file will contain a list of reads aligned to a genome, and the region they cover.

This pipeline converts that to a BED file, and takes the starting position of each read as a single region.
The presence of some regions on the negative strand and some on the positive strand will be accounted for.

The end result is a histogram of locations within the genome that reads tend to 'pile up' on.

## Outputs

A `.bed` file will be produced - a histogram of alignment starts.