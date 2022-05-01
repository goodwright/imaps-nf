import os
import shutil
from unittest import TestCase
import nextflow

class PrimaryAnalysisRunTests(TestCase):

    def setUp(self):
        self.pipeline = nextflow.Pipeline(
            "subworkflows/primaryanalysis.nf",
            config="subworkflows/primaryanalysis.config"
        )
    

    def test_can_run_pipeline_with_genome_that_has_gzip(self):
        execution = self.pipeline.run(params={
            "fastq": "assets/ultraplex_demux_iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5.fastq.gz",
            "Hs_genome": "assets/human_genome",
        }, profile=["iMaps", "local", "test"])
        self.assertEqual(execution.status, "OK", msg=execution.stdout)
        self.assertEqual(len(execution.process_executions), 25)
    

    def test_can_run_pipeline_with_genome_that_has_no_gzip(self):
        try:
            shutil.copytree("assets/human_genome", "assets/human_genome_no_gzip")
            shutil.rmtree("assets/human_genome_no_gzip/DNA_GUNZIP")
            execution = self.pipeline.run(params={
                "fastq": "assets/ultraplex_demux_iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5.fastq.gz",
                "Hs_genome": "assets/human_genome_no_gzip",
            }, profile=["iMaps", "local", "test"])
            self.assertEqual(execution.status, "OK", msg=execution.stdout)
            self.assertEqual(len(execution.process_executions), 25)
        finally:
            if os.path.exists("assets/human_genome_no_gzip"):
                shutil.rmtree("assets/human_genome_no_gzip")