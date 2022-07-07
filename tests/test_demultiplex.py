import os
import nextflow
import pandas as pd
from tests.base import PipelineTest

class DemultiplexTests(PipelineTest):

    def setUp(self):
        PipelineTest.setUp(self)
        self.pipeline = nextflow.Pipeline(
            "subworkflows/demultiplex.nf",
            config="conf/demultiplex.config"
        )


    def test_can_run_with_csv(self):
        execution = self.pipeline.run(params={
            "annotation": os.path.abspath("assets/clip-annotation.csv"),
            "multiplexed_fastq": os.path.abspath("assets/multiplexed.fastq.gz")
        }, profile=["docker"], location="testlocation")
        self.check_execution_ok(execution, 2)
        self.check_process(execution, "CSV_TO_BARCODE", ["clip-annotation.csv"], ["barcode.csv"])
        self.check_process(
            execution, "ULTRAPLEX",
            ["barcode.csv", "multiplexed.fastq.gz"],
            ["ultraplex_demux_Sample1.fastq.gz", "ultraplex_demux_Sample2.fastq.gz", "versions.yml"]
        )
    

    def test_can_run_with_xlsx(self):
        pd.read_csv("assets/clip-annotation.csv").to_excel("testlocation/annotation.xlsx")
        execution = self.pipeline.run(params={
            "annotation": os.path.abspath("testlocation/annotation.xlsx"),
            "multiplexed_fastq": os.path.abspath("assets/multiplexed.fastq.gz")
        }, profile=["docker"], location="testlocation")
        self.check_execution_ok(execution, 3)
        self.check_process(execution, "XLSX_TO_CSV", ["annotation.xlsx"], ["annotation.xlsx.csv"])
        self.check_process(execution, "CSV_TO_BARCODE", ["annotation.xlsx.csv"], ["barcode.csv"])
        self.check_process(
            execution, "ULTRAPLEX",
            ["barcode.csv", "multiplexed.fastq.gz"],
            ["ultraplex_demux_Sample1.fastq.gz", "ultraplex_demux_Sample2.fastq.gz", "versions.yml"]
        )