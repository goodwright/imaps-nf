import os
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


    def test_can_run_bowtie_build(self):
        execution = self.pipeline.run(params={
            "fasta": os.path.abspath("assets/genome.fasta"),
            "gtf": os.path.abspath("assets/genome.gtf"),
        }, profile=["docker"], location="testlocation")
        self.check_execution_ok(execution, 1)
        self.check_process(execution, "STAR_GENOMEGENERATE", ["genome.fasta", "genome.gtf"], ["star"])



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