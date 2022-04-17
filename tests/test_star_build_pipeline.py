import os
from pathlib import Path
from unittest import TestCase
import nextflow

class StarBuildRunTests(TestCase):

    def setUp(self):
        self.pipeline = nextflow.Pipeline(
            "subworkflows/modules/starbuild.nf",
            config="subworkflows/modules/starbuild.config"
        )


    def test_can_run_pipeline(self):
        execution = self.pipeline.run(params={
            "fasta": os.path.abspath(Path("assets/hs.GRCh38.chr21.fa")),
            "gtf": os.path.abspath(Path("assets/hs.GRCh38.chr21.gtf")),
        }, profile=["iMaps", "local"])
        self.assertEqual(execution.status, "OK", msg=execution.stdout)
        self.assertEqual(len(execution.process_executions), 1)