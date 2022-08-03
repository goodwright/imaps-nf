# UMICollapse

`umicollapse.nf` removes reads with duplicate UMIs in a BAM alignment file.

## Inputs

Required files are:

- `bam` - an alignment file with reads containing a UMI in the header.
- `bai` - a SAMTOOLS index for the BAM file.

Other optional inputs:

- `umi_separator` - if the UMI separator is something other than `rbc:` you should specify it here.

## Processes

### `UMICOLLAPSE`

Reads files as returned by the sequencer will contain many PCR duplicates, which are identified via the UMI sequence appended to them.
If the quantity of the original pre-PCR sequences is important, these need to be collapsed back to the original molecules.
UMICollapse does this by removing all but one read for each unique UMI.

## Outputs

A filtered BAM file will be produced, and a log file of what was done.