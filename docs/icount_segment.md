# iCount Segment

`icount_segment.nf` converts a GTF file into ones formatted according to iCount's needs.

## Inputs

Required files are:

- `gtf` - a genome annotation.
- `fai` - a Faidx index file.

## Processes

### `ICOUNT_SEGMENT`

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

Genome-wise segmentation into regions:

```
|-intergenic-||--UTR5-||--UTR5-||-----CDS-----||-CDS-||-intron-||-UTR3-||-intergenic-|
```

## Outputs

A regions file and a segmentation file will be produced.