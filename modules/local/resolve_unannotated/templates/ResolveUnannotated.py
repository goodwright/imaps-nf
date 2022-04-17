#!/usr/bin/env python

import pandas as pd
import pybedtools as pbt
import csv
from plumbum.cmd import sort
import tempfile

def ReadGtf(segmentation):
    df_segment = pd.read_csv(segmentation,
                             sep='\t',
                             names=['chrom', 'source', 'feature', 'start', 'end', 'name', 'strand', 'name2', 'annotations'],
                             header=None,
                             comment='#',
                             dtype={
                                    "chrom": str,
                                    "source": str,
                                    "feature": str,
                                    "start": int,
                                    "end": int,
                                    "name": str,
                                    "strand": str,
                                    "name2": str,
                                    "annotations": str,
                                    })
    return df_segment


def Fai2Bed(fai):
    df_chromosomes = pd.read_csv(fai, sep='\t', header=None, names=['chr', 'end', 'offset', 'linebases', 'linewidth'])
    df_chromosomes = df_chromosomes[['chr', 'end']].assign(start=0, name='.', score=0)
    df_chromosomes_p = df_chromosomes.copy()
    df_chromosomes_p['strand'] = '+'
    df_chromosomes_p = df_chromosomes_p[['chr', 'start', 'end', 'name', 'score', 'strand']]

    df_chromosomes_m = df_chromosomes.copy()
    df_chromosomes_m['strand'] = '-'
    df_chromosomes_m = df_chromosomes_m[['chr', 'start', 'end', 'name', 'score', 'strand']]

    df_chromosomes = pd.concat([df_chromosomes_p, df_chromosomes_m], ignore_index=True)
    bed_chr = pbt.BedTool.from_dataframe(df_chromosomes).sort()
    return(bed_chr)

def run(filtered_segment, unfiltered_segment, gtf_annotation, fai, outputdir, genic_other=False):
    # Read filtered iCount genomic segment and convert it from GTF to BED format.
    print('Reading genomic segmentation.')
    df_segment = ReadGtf(filtered_segment)
    bed_segment = df_segment.assign(start=df_segment['start']-1, score=0)[['chrom', 'start', 'end', 'feature', 'score','strand']]
    bed_segment = pbt.BedTool.from_dataframe(bed_segment).sort()
    # Read unfiltered iCount genomic segment and convert it from GTF to BED format.
    df_unfiltered = ReadGtf(unfiltered_segment)
    bed_unfiltered = df_unfiltered.assign(start=df_unfiltered['start']-1, score=0)[['chrom', 'start', 'end', 'feature', 'score','strand']]
    bed_unfiltered = pbt.BedTool.from_dataframe(bed_unfiltered).sort()

    # Convert fasta index to BED format - one entry spans one chromosome.
    bed_fai = Fai2Bed(fai)

    # Read annotation GTF, keep only gene-level entries and convert it to BED format.
    print('Getting gene-level annotation...')
    df_annotation = ReadGtf(gtf_annotation)
    df_annotation = df_annotation.loc[df_annotation['feature']=='gene']
    bed_annotation = df_annotation.assign(start=df_annotation['start']-1, score=0)[['chrom', 'start', 'end', 'annotations', 'score', 'strand']]
    bed_annotation = pbt.BedTool.from_dataframe(bed_annotation).sort()

    # Find regions that are unannotated in the iCount genome segmentation.
    print('Getting unannotated regions...')
    bed_missing = bed_fai.subtract(bed_segment, s=True,  nonamecheck=True).sort()
    print(f'Found {len(bed_missing)} unannotated genomic regions.')
    # Map gene annotation to the missing regions using bedtools and convert it to GTF format.
    print('Annotating regions with gene information...')
    bed_missing = bed_missing.map(bed_annotation, s=True, c=4, o='collapse', nonamecheck=True)
    if not genic_other:
        # Intersect missing regions with unfiltered segment to get transcript region
        print('Annotating missing regions in iCount segment with transcript regions...')
        # 1 - Use intersect to split unnanotated regions
        intersect = bed_missing.intersect(bed_unfiltered, s=True, nonamecheck=True).sort()
        # 2 - Annotate
        missingAnnotated = intersect.map(bed_unfiltered, s=True, c=4, o='collapse', nonamecheck=True)
        df_unnanotated = pd.read_csv(missingAnnotated.fn, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'annotations', 'feature'])
        df_unnanotated = df_unnanotated.assign(start=df_unnanotated['start'] + 1,
                                               source='.',
                                               name2='.')
    else:
        print('Annotationg missing regions in iCount segment as "genic_other".')
        # If genic_other flag is enabled, missing regions are annotated as "genic_other"
        df_unnanotated = pd.read_csv(bed_missing.fn, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'annotations'])
        # Feature is genic_other.
        df_unnanotated = df_unnanotated.assign(feature='genic_other',
                                           start=df_unnanotated['start'] + 1,
                                           source='.',
                                           name2='.')
    df_unnanotated = df_unnanotated[['chrom', 'source', 'feature', 'start', 'end', 'name', 'strand', 'name2', 'annotations']]
    #Add missing regions to original iCount segment.
    print('Adding annotated missisng regions to iCount segment...')
    df_segment = pd.concat([df_segment, df_unnanotated], ignore_index=True)
    # Sort GTF segment and write it to file
    if genic_other:
        identifier = 'genic_other'
    else:
        identifier = 'annotated'
    outfile = f"{outputdir}/sorted.{identifier}.{filtered_segment.split('/')[-1].replace('.gz', '')}"
    with tempfile.NamedTemporaryFile(mode='w') as tmpfile:
        df_segment.to_csv(tmpfile.name, index=False, header=False, sep='\t', quoting=csv.QUOTE_NONE)
        cmd = (sort["-t\t", "-k1,1", "-k4,4n", tmpfile.name]) > outfile
        print(cmd())
    print(f'Saved the segment as {outfile}')
    return 0

if __name__ == '__main__':
    args = "${task.ext.args}".strip().split()
    run(
        "$filtered_segmentation",
        "$unfiltered_segmentation",
        "$gtf",
        "$fai",
        "./",
        "-go" in args or "--genic_other" in args,
    )
