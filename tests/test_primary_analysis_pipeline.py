import os
from unittest import TestCase
import nextflow

class PrimaryAnalysisRunTests(TestCase):

    def setUp(self):
        self.pipeline = nextflow.Pipeline(
            "subworkflows/primaryanalysis.nf",
            config="subworkflows/primaryanalysis.config"
        )
    

    def test_can_run_pipeline(self):
        execution = self.pipeline.run(params={
            "fastq": "assets/ultraplex_demux_iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5.fastq.gz",
            "Hs_genome": "assets/human_genome",
        }, profile=["iMaps", "local", "test"])
        self.assertEqual(execution.status, "OK", msg=execution.stdout)
        self.assertEqual(len(execution.process_executions), 25)