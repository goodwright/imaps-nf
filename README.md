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

### Demultiplex

A subworkflow which demultiplexes a multiplexed FASTQ reads file using barcodes
in a provided spreadsheet, and performs basic quality checks on the produced
reads using FASTQC. [Full documentation](docs/Demultiplex.md).

### Stand-alone Pipelines

Some of the individual modules have simple workflow wrappers around them,
allowing them to be run directly. This includes:

- `fastqc.nf` - runs the FASTQC module.

### Prepare Genome

#### Inputs

Genome FASTA - a FASTA file containing the DNA sequence of the entire genome. This can be gzipped. Different chromosomes should be their own sequence.

Genome Annotation - a GTF file containing the official, exhaustive annotation of the provided genome FASTA.

smRNA FASTA - a FASTA file listing all the small RNA genes, such as tRNAs.

#### Processes

##### DNA_GUNZIP

If the genome FASTA file is gzipped, this process will uncompress it with gunzip and output the raw FASTA file.

##### RNA_GUNZIP

If the smRNA FASTA file is gzipped, this process will uncompress it with gunzip and output the raw FASTA file.

##### STAR_GENOMEGENERATE

Takes the raw genome and its annotation GTF and builds a STAR genome index - a compressed binary representation of the genome that STAR can align reads to.

##### BOWTIE_BUILD

Takes the raw smRNA FASTA and builds a Bowtie genome index - a compressed binary representation of the smRNA that Bowtie can align reads to.

##### SAMTOOLS_FAIDX

Creates a FAI genome index of the raw genome file - a small reference file which gives the start positions of the sequences (chromosomes in this case) within it so that they can be accessed quickly and memory-efficiently.



##### FIND_LONGEST_TRANSCRIPT

Searches through the GTF annotation and, for each gene, identifies the longest protein coding transcript. This is output as a TXT file where each line contains an ENSEMBL gene ID, the ENSEMBL transcript ID of the transcript from this gene which is the longest, and the length of that longest transcript.

Genes can have multiple transcripts (each of which may be split into multiple exons) so knowing which transcript is the 'representative' one to use is useful in transcriptomics analysis of the genome later.

An index file is also produced giving the locations of these genes in the genome FASTA.

##### FILTER_GTF

In the attributes column of a GTF file, one of the attributes is 'tag', and often annotators will tag some lines as 'basic' meaning that it is a representative transcript for a particular gene. It is often useful to only use these representative transcripts, so this process will filter a GTF to keep only gene lines, and other lines tagged as basic. It only does this if the 'basic' tag is used in the GTF file.

Transcripts also sometimes have a 'TSL' - transcript support level - which can be 1 (high confidence) or above (decreasing confidence). For genes whose transcripts have these identifiers, this process will further filter the transcripts such that only TSL1 and TSL2 are kept.

##### RAW_ICOUNT_SEGMENT/FILTERED_ICOUNT_SEGMENT

iCount segment is a tool which takes a GTF annotation file, and produces two modified versions of it that iCount needs for its downstream analysis - a 'segmentation' GTF file and a 'regions' GTF file.

A segmentation file essentially 'fills in' certain unannotated portions of the genome. It creates GTF records for intergenic regions, and ensures that every section of transcripts have an annotation (mostly by adding records for introns). Records that were previously just called 'exons' will be given a more detailed name.

A regions file is a flat representation of the genome that gives every nucleotide one and only one 'region'.


Initial GTF:
```
              |-----------gene1(G1)-----------|
              |--------transcript1(A)---------|
              |-exon--|         |----exon-----|
                       |------------------gene2(G2)--------------------|
                       |-----------------transcript2(B)----------------|
                       |-exon--|        |----exon----|          |-exon-|
```

Transcript-wise segmentation:
```
|-intergenic-|
              |-----------gene1(G1)-----------|
              |--------transcript1(A)---------|
              |-UTR5--||-intron||-----CDS-----|
                       |------------------gene2(G2)--------------------|
                       |-----------------transcript2(B)----------------|
                       |-UTR5--||intron||-----CDS----||-intron-||-UTR3-|
                                                                        |-intergenic-|
```

Ggenome-wise segmentation into regions:

```
|-intergenic-||--UTR5-||--UTR5-||-----CDS-----||-CDS-||-intron-||-UTR3-||-intergenic-|
```

This is performed twice - one on the initial GTF, once on the filtered GTF.


##### RESOLVE_UNANNOTATED/RESOLVE_UNANNOTATED_GENIC_OTHER

When performing iCount segmentation on GTF annotation, the whole genome is split into the following genomic regions: CDS, UTR3, UTR5, ncRNA, intron, intergenic.

If the segmentation is performed on filtered genomic annotation, some regions remain unnanotated. This is due to the fact that the annotation on a "gene" level covers all transcripts related to this gene, including transcripts with low transcript_support_level, which were removed by filtering. As a result, we get unannotated regions, which are covered by a gene, but lack transcripts that are required for segmentation.

This process finds unannotated regions in the iCount genomic segment (regions.gtf), annotates them according to the corresponding gene, and adds them to the original segment under the feature "genic_other" (if `--genic_other` is provided) or else with the missing transcript segment.

This process is run twice, once with `--genic_other`, once without.

### Demultiplex and Analyse

Uses the Demultiplex subworkflow to split a multiplexed reads file into its
component sample reads files, then performs the Primary CLIP Analysis
subworkflow on each downstream reads file.

## Subworkflows

### Demultiplex

Takes a multiplexed reads file, and a CSV file describing the different samples
it contains, and demultiplexes them using Ultraplex. The reads files produced
are then quality checked with FASTQC.

### Primary CLIP Analysis

Takes a demultiplexed reads file and performs the primary CLIP analysis workflow
on it.

## Style Guide

### Processes/Modules

Local modules should contain only the process definition, and only set directives for `tag`, `label` and `container` - other directives can be set in config.

All modules should be accompanied by a descriptive `meta.yml`, which contains the same fields as standard nf-core yaml files.

All modules should obtain any command line arguments from the `ext.args` directive, set in config.

All modules should use named outputs, including a `versions.yml` file that follows current nf-core conventions.

Non-local modules should match the root level `modules.json`.

### Workflows/subworkflows

`addParams` should not be used - process-specific params should be set in config.

Comments should be used very liberally, ideally before each process call. It is often quite unintutiive what a channel definition (for example) is doing, and comments are enormously helplful.

All variables representing channels should use the `ch_` prefix.

Subworkflows that just wrap a single module should be kept in the `modules` subdirectory of `subworkflows`.

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
