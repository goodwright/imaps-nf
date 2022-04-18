import re
import os

with open(".nextflow.log") as f:
    log = f.read()

matches = re.findall("\[(../......)\] Submitted process > (.+?) ", log)
for match in matches:
    dirs = match[0].split("/")
    dir2 = [d for d in os.listdir("work/" + dirs[0]) if d.startswith(dirs[1])][0]
    print(match[1], os.path.abspath(f"work/{dirs[0]}/{dir2}"))