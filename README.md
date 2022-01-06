# imaps-nf

## Modules

### `GUNZIP`

This module takes as its input any gzipped file, and has output channels for the
decompressed file, and a text file containing the version.

### `STAR_GENOMEGENERATE`

Uses STAR to generate a STAR-formatted genome index from a genome file (in
FASTA format) and an annotation file (in GTF format). This genome index is a
compressed representation of the information in the previous two files.

The module has output channels for the index file, and a text file containing
the version.

### `SAMTOOLS_FAIDX`

Takes a genome FASTA file and creates a FAI index from it, with information on
the chromosomes it contains. It also has a second output channel for the
version.

### `ICOUNT_SEGMENT`

Runs `iCount-Mini segment` on a GTF genome annotation file and the corresponding
FAI index file, to create a segmentation file. The resulting segmentation GTF
file is output, along with a regions file (compressed) and a version file.

### `CSV_TO_BARCODE`

This module takes a CSV file describing a multiplexed reads file, and produces a
CSV file containing just the barcode information. That is, for each sample it
extracts the sample name and the barcode sequence that identifies it, formatted
for the needs of Ultraplex. It uses pandas to do this.

### `ULTRAPLEX`

A wrapper around the command-line tool Ultraplex, which demultiplexes a reads
file into separate reads files based on barcode sequences.

Its inputs are the demultiplexed file (technically a tuple of a 'meta' object
with information about the context, and the reads file path) and a barcodes CSV
which maps 5' barcodes with sample names. Ultraplex will create a reads file for
every barcode that is represented in the multiplexed file, and a reads file for
all those which don't match. It also produces a log file.

### `FASTQC`

A wrapper around the fastqc tool, which takes a reads file and validates it,
producing a HTML report of the quality of the reads.

Its input channel takes in reads as a tuple of a meta object and the reads file
object. The meta object contains flags for whether to use single- or double-end
processing, among other things.

There is an output channel for the produced HTML file, and one for the produced
zip file - each of which is again a tuple of meta information and the actual
file. There is also an output channel for the version as a text file.

### `TRIMGALORE`

A wrapper around TrimGalore, a tool for trimming reads files. Like FASTQC, its
input channel is a tuple of meta object and reads file.

There are output channels for the reads, the report log (both as meta tuples)
and the version.

## Workflows

### Prepare Genome

Generates descriptive files and indexes for a particular genome, starting from
the raw genome in FASTA format, and an accompanying annotation GTF file.

Specifically, it will generate a STAR index file, a FAI index file, and some
segmentation GTF annotation files.