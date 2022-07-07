import os
import re
import json
import yaml
import glob
import subprocess
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
    

    def test_meta_files(self):
        """Every local module should have a valid meta.yml file."""

        for directory in os.listdir(Path("modules/local")):
            if "." not in directory:
                self.assertIn(
                    "meta.yml", os.listdir(Path(f"modules/local/{directory}")),
                    msg=f"{directory} doesn't have a meta.yml file"
                )
                with open(f"modules/local/{directory}/meta.yml") as f:
                    meta = yaml.safe_load(f)
                    for key in ["name", "description", "input", "output"]:
                        self.assertIn(
                            key, meta,
                            msg=f"{directory}'s meta.yml is missing {key}"
                        )
                    for key in ["input", "output"]:
                        for property in meta[key]:
                            property_name = list(property.keys())[0]
                            for field in ["type", "description"]:
                                self.assertIn(
                                    field, property[property_name],
                                    msg=f"{property_name} in {directory}'s meta.yml {key} is missing {field}"
                                )



class PipelineFilesTests(RepoTest):

    def get_pipeline_names(self):
        """What are the pipelines that should be runnable by themselves?"""

        command = 'grep -r -i --include="*.nf" "workflow {" subworkflows'
        stdout = subprocess.run(command, shell=True, stdout=subprocess.PIPE).stdout
        subworkflow_names = [re.search(
            r"\/([a-z_]+?)\.nf", line
        )[1] for line in stdout.decode().splitlines()]
        workflow_names = [f[:-3] for f in os.listdir("workflows") if f.endswith(".nf")]
        return sorted(set(workflow_names + subworkflow_names))


    def test_every_pipeline_has_conf_file(self):
        pipeline_names = self.get_pipeline_names()
        for name in pipeline_names:
            self.assertIn(f"{name}.config", os.listdir("conf"), msg=f"{name}.nf has no config file")
    

    def test_every_pipeline_has_schema_file(self):
        pipeline_names = self.get_pipeline_names()
        for name in pipeline_names:
            self.assertIn(f"{name}.json", os.listdir("schema"), msg=f"{name}.nf has no schema file")
        for name in pipeline_names:
            with open(Path(f"schema/{name}.json")) as f:
                schema = json.load(f)
            for category_name, category in schema.items():
                self.assertIn("name", category, msg=f"{name}.json {category_name} category has no name")
                self.assertIn("description", category, msg=f"{name}.json {category_name} category has no description")
                self.assertIn("properties", category, msg=f"{name}.json {category_name} category has no properties")
                for property_name, property in category["properties"].items():
                    self.assertIn("type", property, msg=f"{name}.json {property_name} property has no type")
                    self.assertIn("description", property, msg=f"{name}.json {property_name} property has no description")