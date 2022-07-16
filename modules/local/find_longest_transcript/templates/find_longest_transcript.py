#!/usr/bin/env python

import re

gtf = "$gtf"
t_types = ["transcript_type", "transcript_biotype"]

# For each transcript, get the transcript ID, gene ID, total CDS length and
# total exon length, using only protein coding transcripts
transcripts = {}
with open(gtf) as f:
    for line in f:
        if line[0] == "#": continue
        values = line.split("\t")
        if any(f'{t} "protein_coding"' in values[8] for t in t_types):
            gene_id = re.search(r"gene_id \"(.+?)\";", values[8])[1]
            transcript_id = re.search(r"transcript_id \"(.+?)\";", values[8])[1]
            if transcript_id not in transcripts:
                transcripts[transcript_id] = {
                    "id": transcript_id, "gene_id": gene_id,
                    "exon_length": 0, "cds_length": 0
                }
            length = int(values[4]) - int(values[3]) + 1
            if values[2] == "CDS":
                transcripts[transcript_id]["cds_length"] += length
            if values[2] == "exon":
                transcripts[transcript_id]["exon_length"] += length
transcripts = list(transcripts.values())
print(f"There are {len(transcripts)} protein coding transcripts")

# Group the transcripts by gene ID
genes = {}
transcripts.sort(key=lambda t: [t["gene_id"], -t["exon_length"], t["id"]])
for transcript in transcripts:
    if transcript["gene_id"] not in genes: genes[transcript["gene_id"]] = []
    genes[transcript["gene_id"]].append(transcript)
print(f"These belong to {len(genes)} genes")

# Get the longest transcript by gene using CDS length, then exon length as a tie
# breaker, then transcript ID as a tie breaker for that
transcript_ids = []
for gene_id, gene_transcripts in genes.items():
    gene_transcripts.sort(key=lambda t: t["gene_id"])
    gene_transcripts.reverse()
    gene_transcripts.sort(key=lambda t: [-t["cds_length"], -t["exon_length"]])
    transcript_ids.append(gene_transcripts[0]["id"] + "\\n")

# Save to file
print("Saving longest transcript per gene...")
with open(f"{gtf}.longest_transcript.txt", "w") as f:
    f.writelines(transcript_ids)

# Create a genome index of longest transcripts
transcripts.sort(key=lambda t: [t["gene_id"], -t["exon_length"], t["id"]])
lines = []
total = 0
for t in transcripts:
    lines.append(f"{t['id']}\t{t['exon_length']}\t{total}\t80\t81\\n")
    total += t['exon_length'] + 1
with open(f"{gtf}.transcriptome_index.fa.fai", "w") as f:
    f.writelines(lines)