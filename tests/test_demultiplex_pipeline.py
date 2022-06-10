import os
from pathlib import Path
from tests.base import PipelineTest
import nextflow

class DemultiplexRunTests(PipelineTest):

    def setUp(self):
        PipelineTest.setUp(self)
        self.pipeline = nextflow.Pipeline(
            "subworkflows/demultiplex.nf",
            config="subworkflows/demultiplex.config"
        )


    def test_can_run_pipeline_with_csv(self):
        execution = self.pipeline.run(params={
            "annotation": os.path.abspath(Path("assets/TEST_ANNOTATION.csv")),
            "multiplexed_fastq": os.path.abspath(Path("assets/SmB_multiplexed.fq.gz"))
        }, profile=["iMaps", "local"], location="testlocation")
        self.check_execution_ok(execution, 4)
        self.check_results("fastqc", [
            "iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_4_fastqc.zip",
            "iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_4_fastqc.html",
            "iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5_fastqc.zip",
            "iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5_fastqc.html"
        ])
        self.check_results("ultraplex", [
            "ultraplex_demux_iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5.fastq.gz",
            "ultraplex_demux_iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_4.fastq.gz",
            r"ultraplex_\d+\.\d+\.log"
        ])
        self.check_results("csv_to_barcode", ["barcode.csv"])
    

    def test_can_run_pipeline_with_xlsx(self):
        execution = self.pipeline.run(params={
            "annotation": os.path.abspath(Path("assets/TEST_ANNOTATION.xlsx")),
            "multiplexed_fastq": os.path.abspath(Path("assets/SmB_multiplexed.fq.gz"))
        }, profile=["iMaps", "local"], location="testlocation")
        self.check_execution_ok(execution, 5)
        self.check_results("fastqc", [
            "iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_4_fastqc.zip",
            "iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_4_fastqc.html",
            "iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5_fastqc.zip",
            "iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5_fastqc.html"
        ])
        self.check_results("ultraplex", [
            "ultraplex_demux_iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5.fastq.gz",
            "ultraplex_demux_iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_4.fastq.gz",
            r"ultraplex_\d+\.\d+\.log"
        ])
        self.check_results("csv_to_barcode", ["barcode.csv"])
        self.check_results("xlsx_to_csv", ["TEST_ANNOTATION.xlsx.csv"])
    

    def test_can_run_pipeline_single_end_fastqc(self):
        execution = self.pipeline.run(params={
            "annotation": os.path.abspath(Path("assets/TEST_ANNOTATION.csv")),
            "multiplexed_fastq": os.path.abspath(Path("assets/SmB_multiplexed.fq.gz")),
            "fastqc_single_end": "true"
        }, profile=["iMaps", "local"], location="testlocation")
        self.check_execution_ok(execution, 4)
        self.check_results("fastqc", [
            "iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_4_fastqc.zip",
            "iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_4_fastqc.html",
            "iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5_fastqc.zip",
            "iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5_fastqc.html"
        ])
        self.check_results("ultraplex", [
            "ultraplex_demux_iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5.fastq.gz",
            "ultraplex_demux_iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_4.fastq.gz",
            r"ultraplex_\d+\.\d+\.log"
        ])
        self.check_results("csv_to_barcode", ["barcode.csv"])
    

    def test_can_run_pipeline_paired_end_fastqc(self):
        execution = self.pipeline.run(params={
            "annotation": os.path.abspath(Path("assets/TEST_ANNOTATION.csv")),
            "multiplexed_fastq": os.path.abspath(Path("assets/SmB_multiplexed.fq.gz")),
            "fastqc_single_end": "false"
        }, profile=["iMaps", "local"], location="testlocation")
        self.check_execution_ok(execution, 4)
        self.check_results("fastqc", [
            "iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_4_1_fastqc.zip",
            "iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_4_1_fastqc.html",
            "iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5_1_fastqc.zip",
            "iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5_1_fastqc.html"
        ])
        self.check_results("ultraplex", [
            "ultraplex_demux_iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5.fastq.gz",
            "ultraplex_demux_iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_4.fastq.gz",
            r"ultraplex_\d+\.\d+\.log"
        ])
        self.check_results("csv_to_barcode", ["barcode.csv"])