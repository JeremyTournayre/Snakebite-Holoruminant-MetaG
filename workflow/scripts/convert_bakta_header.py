#!/usr/bin/env python3
"""
Bakta MAG Protein Header Formatter + TSV Fixer
----------------------------------------------

This script:
1️⃣ Fixes Bakta MAG protein FASTA headers with bin, contig, gene number, protein ID, description, and alignment stats.
2️⃣ Updates the corresponding TSV file so that 'Sequence Id' matches the enriched FASTA headers.

Requirements:
- Python packages: Biopython, pandas
- BLAT installed and available in PATH
"""

import gzip
from Bio import SeqIO
import pandas as pd
import subprocess
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sample", required=True)
args = parser.parse_args()
sample = args.sample

print(f"\n==============================")
print(f"Processing sample: {sample}")
print(f"==============================")

# -------------------------
# Define file paths
# -------------------------
bakta_faa = f"results/mag_annotate/bakta_mags/bakta_{sample}.faa"
bakta_ffn = f"results/mag_annotate/bakta_mags/bakta_{sample}.ffn"
contigs_fasta = f"results/assemble/magscot/{sample}.fa.gz"
output_fasta = f"results/mag_annotate/bakta_mags/bakta_{sample}_fixed_headers.faa"
blat_out = f"results/mag_annotate/bakta_mags/bakta_{sample}_vs_contigs.psl"
contigs_tmp = f"{sample}_contigs_tmp.fa"
tsv_in = f"results/mag_annotate/bakta_mags/bakta_{sample}.tsv"
tsv_out = f"results/mag_annotate/bakta_mags/bakta_{sample}_faa_fixed.tsv"

# -------------------------
# 1️⃣ Map contig → bin
# -------------------------
contig2bin = {}
with gzip.open(contigs_fasta, "rt") as fin:
    for rec in SeqIO.parse(fin, "fasta"):
        header = rec.id
        if header.startswith("bin"):
            parts = header.split("_", 2)
            bin_id = parts[0].replace("bin","")
            contig_name = "_".join(parts[1:])
            contig2bin[contig_name] = bin_id
        else:
            contig2bin[header] = "NA"

# -------------------------
# 2️⃣ Decompress contigs if needed
# -------------------------
if not os.path.exists(contigs_tmp):
    with gzip.open(contigs_fasta, "rt") as fin, open(contigs_tmp, "w") as fout:
        for line in fin:
            fout.write(line)

# -------------------------
# 3️⃣ BLAT alignment
# -------------------------
if not os.path.exists(blat_out):
    print("Running BLAT...")
    subprocess.run(
        f"blat {contigs_tmp} {bakta_ffn} {blat_out} -out=blast8",
        shell=True, check=True)
else:
    print("BLAT output exists — skipping.")

# -------------------------
# 4️⃣ Parse BLAT output
# -------------------------
cols = ["query","subject","perc_id","aln_len","mismatches","gap_opens",
        "q_start","q_end","s_start","s_end","evalue","bitscore"]

blat_df = pd.read_csv(blat_out, sep="\t", comment="#", header=None)
blat_df.columns = cols

best_hits = (
    blat_df.sort_values("bitscore", ascending=False)
            .drop_duplicates("query")
            .set_index("query")
)

# -------------------------
# 5️⃣ Write final FASTA with standardized headers
# -------------------------
gene_counter = {}
with open(output_fasta, "w") as fout:
    for rec in SeqIO.parse(bakta_faa, "fasta"):
        prot_id = rec.id
        prot_len = len(rec.seq)
        description = rec.description.replace(prot_id, "").strip()

        if prot_id in best_hits.index:
            hit = best_hits.loc[prot_id]
            contig = hit["subject"]

            perc_id = round(hit["perc_id"], 1)
            aln_len = hit["aln_len"]
            prot_len_nt = prot_len * 3
            coverage = round((aln_len / prot_len_nt) * 100, 1)

            if perc_id < 100 or coverage < 100:
                contig = f"{sample}:binNA@contigNA"
        else:
            contig = f"{sample}:binNA@contigNA"
            coverage = 0.0
            perc_id = 0.0

        gene_counter[contig] = gene_counter.get(contig, 0) + 1
        gene_num = gene_counter[contig]

        header = f"{contig}_{gene_num} ID={prot_id} {description} match={perc_id}%_{coverage}%"
        fout.write(f">{header}\n{rec.seq}\n")

if os.path.exists(contigs_tmp):
    os.remove(contigs_tmp)

# 6️⃣ Fix TSV to match .faa headers (robuste)
# -------------------------

# Lire le fichier et identifier le header réel
with open(tsv_in, "r") as f:
    lines = f.readlines()

# Chercher la dernière ligne commençant par "#Sequence Id"
header_line_idx = None
for i, line in enumerate(lines):
    if line.startswith("#Sequence Id"):
        header_line_idx = i

if header_line_idx is None:
    raise ValueError("Impossible de trouver la ligne de header '#Sequence Id' dans le TSV")

# Extraire le header et les données
header = lines[header_line_idx].strip().lstrip("#").split("\t")
data_lines = lines[header_line_idx + 1:]

# Lire les données dans pandas
from io import StringIO
df = pd.read_csv(StringIO("".join(data_lines)), sep="\t", names=header)

# Build mapping Locus Tag -> enriched header
prot2header = {}
for rec in SeqIO.parse(output_fasta, "fasta"):
    parts = rec.description.split("ID=")
    header_name = parts[0].strip()
    prot_id = parts[1].split()[0]
    prot2header[prot_id] = header_name

# Remplacer Sequence Id par header enrichi
df["Sequence Id"] = df["Locus Tag"].map(prot2header)

# Sauvegarder TSV corrigé
df = df.dropna(subset=["Sequence Id"])
df.to_csv(tsv_out, sep="\t", index=False)

print(f"Finished sample {sample} ✅")