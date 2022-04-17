import os
from unittest import TestCase
import nextflow

class FastqcRunTests(TestCase):

    def setUp(self):
        self.pipeline = nextflow.Pipeline(
            "subworkflows/modules/fastqc.nf",
            config="subworkflows/modules/fastqc.config"
        )


    def test_can_run_pipeline(self):
        execution = self.pipeline.run(params={
            "fastq": "assets/ultraplex_demux_iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5.fastq.gz",
        }, profile=["iMaps", "local"])
        self.assertEqual(execution.status, "OK", msg=execution.stdout)
        self.assertEqual(len(execution.process_executions), 1)