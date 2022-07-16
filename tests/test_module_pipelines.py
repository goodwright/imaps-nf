import os
from pathlib import Path
import nextflow
from tests.base import PipelineTest

class FaidxTests(PipelineTest):

    def setUp(self):
        PipelineTest.setUp(self)
        self.pipeline = nextflow.Pipeline(
            "subworkflows/modules/faidx.nf",
            config="conf/faidx.config"
        )


    def test_can_run_faidx(self):
        execution = self.pipeline.run(params={
            "fasta": os.path.abspath("assets/genes.fasta"),
        }, profile=["docker"], location="testlocation")
        self.check_execution_ok(execution, 1)
        self.check_process(execution, "SAMTOOLS_FAIDX", ["genes.fasta"], ["genes.fasta.fai"])



class BowtieBuildTests(PipelineTest):

    def setUp(self):
        PipelineTest.setUp(self)
        self.pipeline = nextflow.Pipeline(
            "subworkflows/modules/bowtie_build.nf",
            config="conf/bowtie_build.config"
        )


    def test_can_run_bowtie_build(self):
        execution = self.pipeline.run(params={
            "fasta": os.path.abspath("assets/genome.fasta"),
        }, profile=["docker"], location="testlocation")
        self.check_execution_ok(execution, 1)
        self.check_process(execution, "BOWTIE_BUILD", ["genome.fasta"], ["bowtie"])



class BowtieAlignTests(PipelineTest):

    def setUp(self):
        PipelineTest.setUp(self)
        self.pipeline = nextflow.Pipeline(
            "subworkflows/modules/bowtie_align.nf",
            config="conf/bowtie_align.config"
        )


    def test_can_run_bowtie_align(self):
        execution = self.pipeline.run(params={
            "fastq": os.path.abspath("assets/demultiplexed.fastq.gz"),
            "index": os.path.abspath("assets/bowtie"),
        }, profile=["docker"], location="testlocation")
        self.check_execution_ok(execution, 1)
        self.check_process(
            execution, "BOWTIE_ALIGN", ["demultiplexed.fastq.gz", "bowtie"],
            ["demultiplexed.fastq.gz.bam", "demultiplexed.fastq.gz.unmapped.fastq.gz"]
        )



class StarBuildTests(PipelineTest):

    def setUp(self):
        PipelineTest.setUp(self)
        self.pipeline = nextflow.Pipeline(
            "subworkflows/modules/star_build.nf",
            config="conf/star_build.config"
        )


    def test_can_run_star_build(self):
        execution = self.pipeline.run(params={
            "fasta": os.path.abspath("assets/genome.fasta"),
            "gtf": os.path.abspath("assets/genome.gtf"),
        }, profile=["docker"], location="testlocation")
        self.check_execution_ok(execution, 1)
        self.check_process(execution, "STAR_GENOMEGENERATE", ["genome.fasta", "genome.gtf"], ["star"])



class StarAlignTests(PipelineTest):

    def setUp(self):
        PipelineTest.setUp(self)
        self.pipeline = nextflow.Pipeline(
            "subworkflows/modules/star_align.nf",
            config="conf/star_align.config"
        )


    def test_can_run_star_align(self):
        execution = self.pipeline.run(params={
            "fastq": os.path.abspath("assets/demultiplexed.fastq.gz"),
            "index": os.path.abspath("assets/star"),
            "gtf": os.path.abspath("assets/genome.gtf"),
        }, profile=["docker"], location="testlocation")
        self.check_execution_ok(execution, 1)
        self.check_process(
            execution, "STAR_ALIGN", ["demultiplexed.fastq.gz", "star", "genome.gtf"],
            ["demultiplexed.fastq.gz.Aligned.sortedByCoord.out.bam", "demultiplexed.fastq.gz.Aligned.toTranscriptome.out.bam"]
        )



class FastqcTests(PipelineTest):

    def setUp(self):
        PipelineTest.setUp(self)
        self.pipeline = nextflow.Pipeline(
            "subworkflows/modules/fastqc.nf",
            config="conf/fastqc.config"
        )


    def test_can_run_fastqc(self):
        execution = self.pipeline.run(params={
            "fastq": os.path.abspath("assets/demultiplexed.fastq.gz"),
        }, profile=["docker"], location="testlocation")
        self.check_execution_ok(execution, 1)
        self.check_process(
            execution, "FASTQC", ["demultiplexed.fastq.gz"],
            ["demultiplexed.fastq.gz_fastqc.html", "demultiplexed.fastq.gz_fastqc.zip"]
        )



class TrimgaloreTests(PipelineTest):

    def setUp(self):
        PipelineTest.setUp(self)
        self.pipeline = nextflow.Pipeline(
            "subworkflows/modules/trimgalore.nf",
            config="conf/trimgalore.config"
        )


    def test_can_run_trimgalore(self):
        execution = self.pipeline.run(params={
            "fastq": os.path.abspath("assets/demultiplexed.fastq.gz"),
        }, profile=["docker"], location="testlocation")
        self.check_execution_ok(execution, 1)
        self.check_process(
            execution, "TRIMGALORE", ["demultiplexed.fastq.gz"],
            ["demultiplexed.fastq.gz_trimmed.fq.gz", "demultiplexed.fastq.gz.fastq.gz_trimming_report.txt"]
        )



class FilterGtfTests(PipelineTest):

    def setUp(self):
        PipelineTest.setUp(self)
        self.pipeline = nextflow.Pipeline(
            "subworkflows/modules/filter_gtf.nf",
            config="conf/filter_gtf.config"
        )
    

    def test_can_run_filter_gtf(self):
        execution = self.pipeline.run(params={
            "gtf": os.path.abspath("assets/genome.basic.gtf"),
        }, profile=["docker"], location="testlocation")
        self.check_execution_ok(execution, 1)
        self.check_process(execution, "FILTER_GTF", ["genome.basic.gtf"], ["post_filtering.genome.basic.gtf"])
        with open(Path("testlocation/results/filter_gtf/post_filtering.genome.basic.gtf")) as f:
            gtf_lines = f.read().splitlines()
        self.assertEqual(len(gtf_lines), 298)
    

    def test_can_run_filter_gtf_without_basic_tags(self):
        execution = self.pipeline.run(params={
            "gtf": os.path.abspath("assets/genome.gtf"),
        }, profile=["docker"], location="testlocation")
        self.check_execution_ok(execution, 1)
        self.check_process(execution, "FILTER_GTF", ["genome.gtf"], ["post_filtering.genome.gtf"])
        with open(Path("testlocation/results/filter_gtf/post_filtering.genome.gtf")) as f:
            gtf_lines = f.read().splitlines()
        self.assertEqual(len(gtf_lines), 600)



class LongestTranscriptTests(PipelineTest):

    def setUp(self):
        PipelineTest.setUp(self)
        self.pipeline = nextflow.Pipeline(
            "subworkflows/modules/find_longest_transcript.nf",
            config="conf/find_longest_transcript.config"
        )
    

    def test_can_run_find_longest_transcript_with_transcript_type(self):
        execution = self.pipeline.run(params={
            "gtf": os.path.abspath("assets/genome.basic.gtf"),
        }, profile=["docker"], location="testlocation")
        self.check_execution_ok(execution, 1)
        self.check_process(
            execution, "FIND_LONGEST_TRANSCRIPT", ["genome.basic.gtf"],
            ["genome.basic.gtf.longest_transcript.txt", "genome.basic.gtf.transcriptome_index.fa.fai"]
        )
        with open(Path("testlocation/results/find_longest_transcript/genome.basic.gtf.longest_transcript.txt")) as f:
            txt_lines = f.read().splitlines()
        self.assertEqual(len(txt_lines), 31)
        with open(Path("testlocation/results/find_longest_transcript/genome.basic.gtf.transcriptome_index.fa.fai")) as f:
            txt_lines = f.read().splitlines()
        self.assertEqual(len(txt_lines), 114)
    

    def test_can_run_find_longest_transcript_with_transcript_biotype(self):
        execution = self.pipeline.run(params={
            "gtf": os.path.abspath("assets/genome.gtf"),
        }, profile=["docker"], location="testlocation")
        self.check_execution_ok(execution, 1)
        self.check_process(
            execution, "FIND_LONGEST_TRANSCRIPT", ["genome.gtf"],
            ["genome.gtf.longest_transcript.txt", "genome.gtf.transcriptome_index.fa.fai"]
        )
        with open(Path("testlocation/results/find_longest_transcript/genome.gtf.longest_transcript.txt")) as f:
            txt_lines = f.read().splitlines()
        self.assertEqual(len(txt_lines), 116)
        with open(Path("testlocation/results/find_longest_transcript/genome.gtf.transcriptome_index.fa.fai")) as f:
            txt_lines = f.read().splitlines()
        self.assertEqual(len(txt_lines), 116)