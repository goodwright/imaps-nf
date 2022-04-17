import os
from unittest import TestCase
import nextflow

class PrepareGenomeRunTests(TestCase):

    def setUp(self):
        self.pipeline = nextflow.Pipeline(
            "workflows/preparegenome.nf",
            config="workflows/preparegenome.config"
        )
    

    def test_can_run_pipeline(self):
        execution = self.pipeline.run(params={
            "fasta": "https://github.com/luslab/nf-core-test-data/raw/main/data/fasta/homosapien-hg37-chr21.fa.gz",
            "gtf": "https://github.com/luslab/nf-core-test-data/raw/main/data/gtf/gencode.v35.chr21.gtf",
            "smrna_fasta": "assets/Homo_sapiens.smallRNA.fa.gz",
        }, profile=["iMaps", "local"])
        self.assertEqual(execution.status, "OK", msg=execution.stdout)
        self.assertEqual(len(execution.process_executions), 11)