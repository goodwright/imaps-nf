# flow-nf

![](https://github.com/goodwright/flow-nf/actions/workflows/main.yml/badge.svg)

## Repository Structure

This repository contains:

1. The custom Flow Nextflow pipelines not available elsewhere.
2. Config for all Flow pipelines, both local and third-party.

There are four types of Nextflow pipeline here:

- Top-level FASTQ consuming pipelines, which take a FASTQ file as input and perform a full analysis on it (e.g primary_clip_analysis.nf).
- Subworkflows that make up the top-level FASTQ consuming pipelines but which are specific to them (e.g. prepare_clip_genome.nf).
- General purpose pipelines used in multiple contexts (e.g. demultiplex.nf).
- Wrappers around single modules (e.g. fastqc.nf).

Pipelines in the first category are found in `workflows/`, while those in the second three categories are found in `subworkflows/`.

Every pipeline that is intended to be run stand-alone should have:

1. A config file in `conf/`.
2. A schema file in `schema/`.
3. A test suite somewhere in `tests`.
4. A markdown documentation file in `docs/` (unless it is a module wrapper).

## Pipelines

### Top-Level

### Subworkflows

### General Purpose

- `demultiplex.nf` - demultiplexes a FASTQ file into individual FASTQ files using Ultraplex.

### Modules

- `faidx.nf` - creates a Samtools `.fai` genome index from a FASTA file.

## Config

Config that applies to all pipeline runs belongs in `base.conf`.
This includes global parameters, global process directives, and resource settings.

Each profile other than `test` will have its own config file, with the settings that apply under that profile.

Each runnable pipeline should have its own config file containing process-specific options for that run, as well as the definition of the `test` profile for that pipeline.

## Schema

Schema JSON files specify the parameters that can be passed to a pipeline at the command line, organised by category.
They are broadly based on the spec of the nf-core schema JSON files, though only fields pertaining to parameters are required.
Important things to specify about parameters are:

- Whether it is a file input or a value input (such as a string or boolean).
- Whether it is a mandatory parameter.
- Any Flow-specific requirements such as it being a multiplexed file. 


## Tests

Each top-level pipeline will have a test file, testing that it can run with all meaningful combination of parameters.

Similarly there is be a test file for the runnable subworkflows, the general purpose pipelines, and the standalone modules.

There is also a test file that checks the structure of the repository itself, to ensure that the correct files are always in place.

## Documentation

Each runnable pipeline (except those which just wrap a single module) has a markdown file explaining in plain English what the pipeline is for, what its inputs are, the general flow of what it does, and what its key outputs are.

This is in addition to the more cursory explanation in this README.