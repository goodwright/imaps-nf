import os
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