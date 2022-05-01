import os
import gzip
import requests
from unittest import TestCase
import nextflow

class PrepareGenomeRunTests(TestCase):

    def setUp(self):
        self.pipeline = nextflow.Pipeline(
            "workflows/preparegenome.nf",
            config="workflows/preparegenome.config"
        )
    

    def test_can_run_pipeline_with_compressed_files(self):
        execution = self.pipeline.run(params={
            "fasta": "https://github.com/luslab/nf-core-test-data/raw/main/data/fasta/homosapien-hg37-chr21.fa.gz",
            "gtf": "https://github.com/luslab/nf-core-test-data/raw/main/data/gtf/gencode.v35.chr21.gtf",
            "smrna_fasta": "assets/Homo_sapiens.smallRNA.fa.gz",
        }, profile=["iMaps", "local"])
        self.assertEqual(execution.status, "OK", msg=execution.stdout)
        self.assertEqual(len(execution.process_executions), 11)
    

    def test_can_run_pipeline_with_uncompressed_files(self):
        try:
            resp = requests.get(
                "https://github.com/luslab/nf-core-test-data/raw/main/data/fasta/homosapien-hg37-chr21.fa.gz"
            )
            with open("homosapien-hg37-chr21.fa", "wb") as f:
                f.write(gzip.decompress(resp.content))
            with open("assets/Homo_sapiens.smallRNA.fa.gz", "rb") as f1:
                with open("Homo_sapiens.smallRNA.fa", "wb") as f2:
                    f2.write(gzip.decompress(f1.read()))
            execution = self.pipeline.run(params={
                "fasta": "homosapien-hg37-chr21.fa",
                "gtf": "https://github.com/luslab/nf-core-test-data/raw/main/data/gtf/gencode.v35.chr21.gtf",
                "smrna_fasta": "Homo_sapiens.smallRNA.fa",
            }, profile=["iMaps", "local"])
            self.assertEqual(execution.status, "OK", msg=execution.stdout)
            self.assertEqual(len(execution.process_executions), 9)
        finally:
            if os.path.exists("Homo_sapiens.smallRNA.fa"):
                os.remove("Homo_sapiens.smallRNA.fa")
            if os.path.exists("homosapien-hg37-chr21.fa"):
                os.remove("homosapien-hg37-chr21.fa")
        