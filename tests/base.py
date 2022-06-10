import os
import re
import shutil
from unittest import TestCase

class PipelineTest(TestCase):

    def setUp(self):
        if os.path.exists("testlocation"):
            shutil.rmtree("testlocation")
        os.mkdir("testlocation")


    def tearDown(self):
        if os.path.exists("testlocation"):
            shutil.rmtree("testlocation")
    

    def check_execution_ok(self, execution, process_count):
        self.assertEqual(execution.status, "OK", msg=execution.stdout)
        self.assertNotIn("WARN: ", execution.stdout, msg=execution.stdout)
        self.assertEqual(len(execution.process_executions), process_count)


    def check_results(self, directory, files):
        results_files = os.listdir(f"testlocation/results/{directory}")

        # All the files given are present
        for f in files:
            self.assertTrue(
                any(re.match(f, result) for result in results_files),
                msg=f"{f} not produced - {results_files}"
            )
        
        # There are no extra files
        for result in results_files:
            self.assertTrue(
                any(re.match(f, result) for f in files),
                msg=f"{result} unexpectedly produced - {results_files}"
            )
        
        # Matching number
        self.assertEqual(len(files), len(results_files))