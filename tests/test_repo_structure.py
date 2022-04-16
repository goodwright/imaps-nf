import os
import json
from pathlib import Path
from unittest import TestCase

class ModuleTests(TestCase):

    longMessage = False

    def load_modules(self):
        with open("modules.json") as f:
            return json.load(f)


    def test_repos_exist(self):
        modules = self.load_modules()
        for repo_name in modules["repos"]:
            self.assertTrue(
                os.path.exists(Path(f"modules/{repo_name}")),
                msg=f"{repo_name} is not in the modules directory"
            )