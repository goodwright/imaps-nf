import os
import re
import shutil
from unittest import TestCase

class PipelineTest(TestCase):

    def setUp(self):
        if os.path.exists("testlocation"): shutil.rmtree("testlocation")
        os.mkdir("testlocation")


    def tearDown(self):
        if os.path.exists("testlocation"):
            shutil.rmtree("testlocation")
    

    def check_execution_ok(self, execution, process_count):
        self.assertEqual(execution.status, "OK", msg=execution.stdout)
        self.assertNotIn("WARN: ", execution.stdout, msg=execution.stdout)
        self.assertEqual(len(execution.process_executions), process_count)
    

    def check_process(self, execution, name, inputs, outputs):
        try:
            process = [p for p in execution.process_executions if p.process == name][0]
        except IndexError:
            raise Exception(f"No process with name {name}")
        self.assertEqual(set(inputs), set(process.input_data(include_path=False)))
        all_outputs = set(process.all_output_data(include_path=False))
        for output in outputs:
            self.assertIn(output, all_outputs)