# imaps-nf

![](https://github.com/goodwright/imaps-nf/actions/workflows/main.yml/badge.svg)

This repository contains the Nextflow workflows, subworkflows and modules which
are used in [iMaps](https://imaps.goodwright.com). They are primarily concerned
with CLIP analysis.

## Quick-start

To run any of the pipelines, use the associated config file and run it with the
profiles `test`, `local` and `iMaps`. For example:

```bash
nextflow -C workflows/demuxandanalyse.config run workflows/demuxandanalyse.nf -profile test,local,iMaps
```

The `test` profile will auto-add all params from the `assets` folder, and the
iMaps profile will set up all the docker config. The `local` profile sets low
resource allowance limits.


## Workflows

### Prepare Genome

Generates descriptive files and indexes for a particular genome, starting from
the raw genome in FASTA format, and an accompanying annotation GTF file.

Specifically, it will generate:

- a STAR index directory using STAR.
- a FAI index file describing the sequences/chromosomes within the genome.
- segmentation and regions annotation files using `iCount segment`.
- longest transcript information.

### Demultiplex and Analyse

Uses the Demultiplex subworkflow to split a multiplexed reads file into its
component sample reads files, then performs the Primary Analysis subworkflow on
each downstream reads file.

## Subworkflows

### Demultiplex

Takes a multiplexed reads file, and a CSV file describing the different samples
it contains, and demultiplexes them using Ultraplex. The reads files produced
are then quality checked with FASTQC.

### Primary Analysis

Takes a demultiplexed reads file and performs the primary CLIP analysis workflow
on it.

## Style Guide

### Processes/Modules

Local modules should contain only the process definition, and only set directives for `tag`, `label` and `container` - other directives can be set in config.

All local modules should be accompanied by a descriptive `meta.yml`, which contains the same fields as standard nf-core yaml files.

All local modules should obtain any command line arguments from the `ext.args` directive, set in config.

All local modules should use named outputs, including a `versions.yml` file that follows current nf-core conventions.

Non-local modules should match the root level `modules.json`.

### Workflows/subworkflows

`addParams` should not be used - process-specific params should be set in config.

Comments should be used very liberally, ideally before each process call. It is often quite unintutiive what a channel definition (for example) is doing, and comments are enormously helplful.

All variables representing channels should use the `ch_` prefix.

### Conf

The directives used for particular processes should be defined in `modules.conf`. This includes `publishDir` and `ext.args`.

Resource management should be handled using the current nf-core convention:

- Processes define their own labels indicating whether they require low, medium or high resources.
- The `base.config` should set the `cpus` etc. directives for these labels.
- This should be sensitive to a `max_cpus` etc. param which can override these directives.

Each pipeline should have its own config file, which imports the base config. Often this will be all it needs, unless its processes use settings different from the defaults set globally.

### Tests

Each pipeline should have its own test suite, checking that the pipeline runs, and that the outputs it produces are sensible.

Test files for these should go in the `assets` folder.

A `test` profile should be maintained that allows the pipelines to run with default params and low resource usage on any device.

The structure and adherence to rules in this document should also be tested.

Tests will run via GitHub actions on every push.

### Docs

In addition to any documentation on iMaps docs websites, the README should contain an overview of all workflows and subworkflows available.

## Testing

The tests need Python to run. To install additional dependencies, run:

```
pip install -r tests/requirements.txt
```

To run the tests:

```
python -m unittest discover tests
```

This may take many minutes.

Or to run one test file:

```
python -m unittest tests.test_repo_structure
```
