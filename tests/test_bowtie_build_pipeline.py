import os
from pathlib import Path
from unittest import TestCase
import nextflow

class BowtieRunTests(TestCase):

    def setUp(self):
        self.pipeline = nextflow.Pipeline(
            "subworkflows/modules/bowtiebuild.nf",
            config="subworkflows/modules/bowtiebuild.config"
        )


    def test_can_run_pipeline(self):
        execution = self.pipeline.run(params={
            "fasta": os.path.abspath(Path("assets/Homo_sapiens.smallRNA.fa")),
        }, profile=["iMaps", "local"])
        self.assertEqual(execution.status, "OK", msg=execution.stdout)
        self.assertEqual(len(execution.process_executions), 1)