# Resolve Unannotated Regions

`resolve_unannotated_regions.nf` fills in empty regions in an iCount regions file.

## Inputs

Required files are:

- `filtered_regions` - an iCount regions GTF from a filtered GTF.
- `unfiltered_regions` - an iCount regions GTF from an unfiltered GTF.
- `gtf` - a genome annotation.
- `fai` - a Faidx index file.

## Processes

### `RESOLVE_UNANNOTATED_REGIONS`

When performing iCount segmentation on a GTF annotation, the whole genome is split into the following genomic regions: CDS, UTR3, UTR5, ncRNA, intron, intergenic.

If the segmentation is performed on filtered genomic annotation, some regions remain unnanotated. This is due to the fact that the annotation on a "gene" level covers all transcripts related to this gene, including transcripts with low transcript_support_level, which were removed by filtering. As a result, we get unannotated regions, which are covered by a gene, but lack transcripts that are required for segmentation.

This process finds unannotated regions in the iCount genomic segment (regions.gtf), annotates them according to the corresponding gene, and adds them to the original segment under the feature with the missing transcript segment.

## Outputs

A fully annotated GTF is produced.