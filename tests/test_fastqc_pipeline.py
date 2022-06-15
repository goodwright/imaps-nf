import os
import nextflow
from tests.base import PipelineTest

class FastqcRunTests(PipelineTest):

    def setUp(self):
        PipelineTest.setUp(self)
        self.pipeline = nextflow.Pipeline(
            "subworkflows/modules/fastqc.nf",
            config="subworkflows/modules/fastqc.config",
        )


    def test_can_run_pipeline(self):
        execution = self.pipeline.run(params={
            "fastq": os.path.abspath(
                "assets/ultraplex_demux_iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5.fastq.gz"
            ),
        }, profile=["iMaps", "local"], location="testlocation")
        self.check_execution_ok(execution, 1)
        self.check_results("fastqc", [
            "iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5_fastqc.zip",
            "iCLIP_SmB_Cal51_NSsiRNA_20130808_LUc21_5_fastqc.html"
        ])