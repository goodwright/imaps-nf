import os
from pathlib import Path
from unittest import TestCase
import nextflow

class FastqcRunTests(TestCase):

    def setUp(self):
        self.pipeline = nextflow.Pipeline(
            "subworkflows/modules/icountsegment.nf",
            config="subworkflows/modules/icountsegment.config"
        )


    def test_can_run_pipeline(self):
        execution = self.pipeline.run(params={
            "gtf": os.path.abspath(Path("assets/hs.GRCh38.chr21.gtf")),
            "fai": os.path.abspath(Path("assets/hs_21.fa.fai")),
        }, profile=["iMaps", "local"])
        self.assertEqual(execution.status, "OK", msg=execution.stdout)
        self.assertEqual(len(execution.process_executions), 1)