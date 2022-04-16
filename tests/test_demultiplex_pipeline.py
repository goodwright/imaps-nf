import os
from unittest import TestCase
import nextflow

class DemultiplexRunTests(TestCase):

    def setUp(self):
        self.pipeline = nextflow.Pipeline(
            "subworkflows/demultiplex.nf",
            config="subworkflows/demultiplex.config"
        )


    def test_can_run_pipeline_with_csv(self):
        execution = self.pipeline.run(params={
            "csv": "assets/TEST_ANNOTATION.csv",
            "multiplexed_fastq": "assets/SmB_multiplexed.fq.gz"
        }, profile=["iMaps", "local"])
        self.assertEqual(execution.status, "OK", msg=execution.stdout)
        self.assertEqual(len(execution.process_executions), 4)
    

    def test_can_run_pipeline_with_xlsx(self):
        execution = self.pipeline.run(params={
            "csv": os.path.abspath("assets/TEST_ANNOTATION.xlsx"),
            "multiplexed_fastq": "assets/SmB_multiplexed.fq.gz"
        }, profile=["iMaps", "local"])
        self.assertEqual(execution.status, "OK", msg=execution.stdout)
        self.assertEqual(len(execution.process_executions), 5)