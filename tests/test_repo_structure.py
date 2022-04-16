import os
import json
from pathlib import Path
from unittest import TestCase

class RepoTest(TestCase):
    longMessage = False

class RepoModuleTests(RepoTest):

    def load_modules(self):
        with open("modules.json") as f:
            return json.load(f)


    def test_repos_exist(self):
        """Every repo mentioned in modules.json should be in modules/"""

        modules = self.load_modules()
        for repo_name in modules["repos"]:
            self.assertTrue(
                os.path.exists(Path(f"modules/{repo_name}")),
                msg=f"{repo_name} is not in the modules directory"
            )
    

    def test_all_repos_modules_json(self):
        """Every repo in modules/ needs an entry in modules.json"""
        
        modules = self.load_modules()
        directories = os.listdir(Path("modules"))
        for directory in directories:
            if directory != "local":
                matching_keys = [
                    k for k in modules["repos"] if k.startswith(directory)
                ]
                self.assertEqual(
                    len(matching_keys), 1,
                    msg=f"{directory} directory in modules isn't in modules.json"
                )
                
    
    def test_modules_exist(self):
        """Every module mentioned in modules.json should be in modules/"""

        modules = self.load_modules()
        for repo_name, repo in modules["repos"].items():
            for module in repo:
                self.assertTrue(
                    os.path.exists(Path(f"modules/{repo_name}/{module}")),
                    msg=f"{module} could not be found in {repo_name}"
                )
    

    def test_all_modules_in_modules_json(self):
        """Every module in modules/ needs an entry in modules.json"""

        modules = self.load_modules()
        import glob
        for repo_name, repo in modules["repos"].items():
            for directory in [f for f in glob.glob(
                f"modules/{repo_name}/**/*.nf", recursive=True
            )]:
                module_name = "/".join(directory[
                    directory.find(repo_name) + len(repo_name) + 1:
                ].split("/")[:-1])
                self.assertIn(
                    module_name, repo,
                    f"The {module_name} module could not be found in "
                    f"modules.json under {repo_name}"
                )



class LocalModuleTests(RepoTest):

    def test_main_files(self):
        """Every local module should have a valid main.nf file."""
        
        for directory in os.listdir(Path("modules/local")):
            if "." not in directory:
                self.assertIn(
                    "main.nf", os.listdir(Path(f"modules/local/{directory}")),
                    msg=f"{directory} doesn't have a main.nf file"
                )