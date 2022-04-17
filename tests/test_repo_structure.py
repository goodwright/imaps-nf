import os
import json
import yaml
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

    TEMP_EXEMPT = [
        "filter_gtf",
        "filter_transcriptome_bam",
        "find_longest_transcript",
        "get_crosslinks",
        "icount_peaks",
        "icount_rnamaps",
        "icount_summary",
        "resolve_unannotated",
    ]

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
            if "." not in directory and directory not in self.TEMP_EXEMPT:
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



class WorkflowTests(RepoTest):

    TEMP_EXEMPT = [
        "ncrna.nf",
        "demuxandanalyse.nf",
        "primaryanalysis.nf",
        "primaryanalysis.nf"
    ]

    def test_config_files(self):
        """Every workflow should have an associated config file."""

        for directory in ["workflows", "subworkflows"]:
            files = os.listdir(directory)
            for f in files:
                if f.endswith(".nf"):
                    self.assertIn(
                        f"{f[:-3]}.config", files,
                        msg=f"{f} doesn't have an associated config file"
                    )
    

    def test_schema_files(self):
        """Every workflow should have a valid associated schema file."""

        for directory in ["workflows", "subworkflows"]:
            files = os.listdir(directory)
            for f in files:
                if f.endswith(".nf"):
                    filename = f"{f[:-3]}.json"
                    self.assertIn(
                        filename, files,
                        msg=f"{f} doesn't have an associated schema file"
                    )
                    with open(Path(f"{directory}/{filename}")) as f:
                        schema = json.load(f)
                    self.assertIn(
                        "description", schema,
                        msg=f"{filename} is missing a description"
                    )
                    self.assertIn(
                        "definitions", schema,
                        msg=f"{filename} is missing a definitions section"
                    )
                    for name, section in schema["definitions"].items():
                        for key in ["title", "description", "properties"]:
                            self.assertIn(
                                key, section,
                                msg=f"The {name} section of {filename} is missing {key}"
                            )
                        for prop_name, property in schema["definitions"][name]["properties"].items():
                            for key in ["type", "format", "required", "description"]:
                                self.assertIn(
                                    key, property,
                                    msg=f"The {prop_name} propery in the {name} section of {filename} is missing {key}"
                                )
    

    def test_no_add_params(self):
        for directory in ["workflows", "subworkflows"]:
            files = os.listdir(directory)
            for f in files:
                if f.endswith(".nf") and f not in self.TEMP_EXEMPT:
                    with open(Path(f"{directory}/{f}")) as fp:
                        self.assertNotIn(
                            "addParams", fp.read(),
                            f"{f} is using addParams - params should be set in modules.conf"
                        )
