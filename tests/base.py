import os
import shutil
from unittest import TestCase

class PipelineTest(TestCase):

    def setUp(self):
        os.mkdir("testlocation")


    def tearDown(self):
        if os.path.exists("testlocation"):
            shutil.rmtree("testlocation")