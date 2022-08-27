# Paraclu Cut

The `paraclu_cut.nf` pipeline cuts down the output of Paraclu to the core peaks.

## Inputs

Required files are:

- `sigxls` - a TSV output file of the main Paraclu pipeline.

## Processes

### `PARACLU_CUT`

The main Paraclu program produces a list of significant crosslinks.
Paraclu Cut cuts this down even further by removing clusters that are too long, or singletons.
It also it removes clusters whose fold-increase in density is too small, and it removes clusters that are contained in larger clusters.

The remaining sites are referred to as 'peaks'.

## Outputs

A compressed TSV file containing the peaks is produced.