import os
from pathlib import Path
from unittest import TestCase
import nextflow

class DemultiplexAndAnalyseRunTests(TestCase):

    def setUp(self):
        self.pipeline = nextflow.Pipeline(
            "workflows/demuxandanalyse.nf",
            config="workflows/demuxandanalyse.config"
        )


    def test_can_run_pipeline(self):
        execution = self.pipeline.run(params={
            "annotation": "assets/TEST_ANNOTATION.csv",
            "multiplexed_fastq": "assets/SmB_multiplexed.fq.gz",
            "Hs_genome": "assets/human_genome"
        }, profile=["iMaps", "local", "test"])
        self.assertEqual(execution.status, "OK", msg=execution.stdout)