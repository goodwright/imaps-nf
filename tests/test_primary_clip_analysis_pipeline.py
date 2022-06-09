import os
import shutil
from unittest import TestCase
import nextflow

class PrimaryClipAnalysisRunTests(TestCase):

    def setUp(self):
        self.pipeline = nextflow.Pipeline(
            "subworkflows/primaryclipanalysis.nf",
            config="subworkflows/primaryclipanalysis.config"
        )
    

    def test_can_run_pipeline_with_genome_that_has_gzip(self):
        execution = self.pipeline.run(params={
            "fastq": "assets/ultraplex_demux_iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5.fastq.gz",
            "Hs_genome": "assets/human_genome",
        }, profile=["iMaps", "local", "test"])
        self.assertEqual(execution.status, "OK", msg=execution.stdout)
        self.assertEqual(len(execution.process_executions), 30)

        # Default UMI Separator is rbc
        for proc in execution.process_executions:
            if "UMITOOLS_DEDUP" in proc.process:
                p1, p2 = proc.hash.split("/")
                subdirs = os.listdir(os.path.join("work", p1))
                subdir = [d for d in subdirs if d.startswith(p2)][0]
                with open(os.path.join("work", p1, subdir, ".command.log")) as f:
                    self.assertIn(
                        "--umi-separator=rbc:", f.read(),
                        "Default umi-separator was not 'rbc'"
                    )
    

    def test_can_run_pipeline_with_genome_that_has_no_gzip(self):
        try:
            shutil.copytree("assets/human_genome", "assets/human_genome_no_gzip")
            os.mkdir("assets/human_genome_no_gzip/inputs")
            os.mkdir("assets/human_genome_no_gzip/inputs/fasta")
            shutil.copy(
                "assets/human_genome_no_gzip/DNA_GUNZIP/homosapien-hg37-chr21.fa",
                "assets/human_genome_no_gzip/inputs/fasta/homosapien-hg37-chr21.fa"
            )
            shutil.rmtree("assets/human_genome_no_gzip/DNA_GUNZIP")
            execution = self.pipeline.run(params={
                "fastq": "assets/ultraplex_demux_iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5.fastq.gz",
                "Hs_genome": "assets/human_genome_no_gzip",
            }, profile=["iMaps", "local", "test"])
            self.assertEqual(execution.status, "OK", msg=execution.stdout)
            self.assertEqual(len(execution.process_executions), 30)
        finally:
            if os.path.exists("assets/human_genome_no_gzip"):
                shutil.rmtree("assets/human_genome_no_gzip")
    

    def test_can_run_pipeline_with_rbc_param(self):
        execution = self.pipeline.run(params={
            "fastq": "assets/ultraplex_demux_iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5.fastq.gz",
            "Hs_genome": "assets/human_genome",
            "umi_separator": "xyz"
        }, profile=["iMaps", "local", "test"])
        
        # Custom UMI Separator is xyz
        for proc in execution.process_executions:
            if "UMITOOLS_DEDUP" in proc.process:
                p1, p2 = proc.hash.split("/")
                subdirs = os.listdir(os.path.join("work", p1))
                subdir = [d for d in subdirs if d.startswith(p2)][0]
                with open(os.path.join("work", p1, subdir, ".command.log")) as f:
                    self.assertIn(
                        "--umi-separator=xyz:", f.read(),
                        msg="Custom umi-separator was not 'xyz'"
                    )