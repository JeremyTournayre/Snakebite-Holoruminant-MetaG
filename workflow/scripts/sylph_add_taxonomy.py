#!/usr/bin/env python3
"""
sylph_add_taxonomy.py
---------------------
Adds GTDB taxonomy to Sylph results by looking it up in the GTDB-Tk summary
files of the custom rumen database. Produces one file per sample.

Usage (Snakemake):
    python sylph_add_taxonomy.py \
        --sylph {input.profiling} \
        --taxo-dir {params.taxo_dir} \
        --output {output}
"""

import csv
import re
import argparse
from pathlib import Path


# Mapping dossier GTDB-Tk -> préfixe skDER
# Nécessaire quand GTDB-Tk a été lancé plusieurs fois avec des suffixes différents
GTDB_FOLDER_TO_SKDER_PREFIX = {
    "700MAGs_fa":    "700MAGs",
    "700MAGs_fasta": "700MAGs",
    "700MAGs_fna":   "700MAGs",
    # les autres : dossier == préfixe skDER, pas besoin de les lister
}


def clean_genome_name(path_or_name: str) -> str:
    """Extracts the base name without known extensions."""
    name = Path(path_or_name).name
    name = re.sub(r'_genomic\.fna\.gz$', '', name)
    name = re.sub(r'\.(fasta|fna|fa)$', '', name, flags=re.IGNORECASE)
    name = re.sub(r'\.ok$', '', name)
    return name


def load_taxonomy(taxo_dir: str) -> dict:
    taxo = {}
    taxo_dir = Path(taxo_dir)

    # ── 1. Fichiers GTDB-Tk standard (gtdbtk.*.summary.tsv) ──
    summary_files = sorted(taxo_dir.rglob("gtdbtk.*.summary.tsv"))
    for summary_file in summary_files:
        folder_name = summary_file.parent.name
        source_prefix = GTDB_FOLDER_TO_SKDER_PREFIX.get(folder_name, folder_name)
        print(f"[sylph_add_taxonomy] Loading {summary_file} -> prefix '{source_prefix}'")
        with open(summary_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                genome_id = row["user_genome"].strip()
                genome_id = re.sub(r'\.ok$', '', genome_id) 
                key = f"{source_prefix}__{genome_id}"
                taxo_str = row["classification"].strip()
                levels = [
                    lvl for lvl in taxo_str.split(";")
                    if lvl and "__" in lvl and lvl.split("__")[1] != ""
                ]
                if levels and key not in taxo:
                    taxo[key] = ";".join(levels)

    # ── 2. Fichiers GTDB r232 bruts (*_taxonomy_r232.tsv) ──
    # Format : RS_GCF_xxx.1\ttaxonomy  (pas de header)
    # skDER :  gtdb_232__GCA_xxx.1_genomic.fna.gz  -> clé : gtdb_232__GCA_xxx.1
    r232_files = sorted(taxo_dir.rglob("*_taxonomy_r232.tsv"))
    for r232_file in r232_files:
        source_prefix = "gtdb_232"
        print(f"[sylph_add_taxonomy] Loading {r232_file} -> prefix '{source_prefix}'")
        with open(r232_file) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                raw_acc, taxo_str = parts[0], parts[1]
                # Supprime le préfixe RS_ ou GB_
                acc = re.sub(r'^(RS_|GB_)', '', raw_acc)
                key = f"{source_prefix}__{acc}"
                levels = [
                    lvl for lvl in taxo_str.split(";")
                    if lvl and "__" in lvl and lvl.split("__")[1] != ""
                ]
                if levels and key not in taxo:
                    taxo[key] = ";".join(levels)

    if not summary_files and not r232_files:
        print(f"[sylph_add_taxonomy] Warning: no taxonomy files found in {taxo_dir}")

    print(f"[sylph_add_taxonomy] {len(taxo)} genomes loaded total")
    return taxo


def process_sylph_file(sylph_file: str, taxo: dict, output: str):
    """Joins a Sylph file with taxonomy and writes the output file."""
    out_rows = []
    not_found = set()

    with open(sylph_file) as f:
        reader = csv.DictReader(f, delimiter="\t")
        original_cols = list(reader.fieldnames)
        for row in reader:
            genome_name = clean_genome_name(row["Genome_file"])
            taxonomy = taxo.get(genome_name, "NA")
            if taxonomy == "NA":
                not_found.add(genome_name)

            tax_levels = {
                "domain": "", "phylum": "", "class": "",
                "order": "", "family": "", "genus": "", "species": ""
            }
            if taxonomy != "NA":
                for lvl in taxonomy.split(";"):
                    if   lvl.startswith("d__"): tax_levels["domain"]  = lvl[3:]
                    elif lvl.startswith("p__"): tax_levels["phylum"]  = lvl[3:]
                    elif lvl.startswith("c__"): tax_levels["class"]   = lvl[3:]
                    elif lvl.startswith("o__"): tax_levels["order"]   = lvl[3:]
                    elif lvl.startswith("f__"): tax_levels["family"]  = lvl[3:]
                    elif lvl.startswith("g__"): tax_levels["genus"]   = lvl[3:]
                    elif lvl.startswith("s__"): tax_levels["species"] = lvl[3:]

            out_row = dict(row)
            out_row["Genome_name"] = genome_name
            out_row["Taxonomy_full"] = taxonomy
            out_row.update(tax_levels)
            out_rows.append(out_row)

    if out_rows:
        new_cols = ["Genome_name", "Taxonomy_full", "domain", "phylum",
                    "class", "order", "family", "genus", "species"]
        fieldnames = original_cols + [c for c in new_cols if c not in original_cols]

        with open(output, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(out_rows)

    print(f"[sylph_add_taxonomy] {len(out_rows)} rows -> {output}")
    if not_found:
        print(f"[sylph_add_taxonomy] Warning: {len(not_found)} genomes without taxonomy: "
              + ", ".join(sorted(not_found)[:5])
              + ("..." if len(not_found) > 5 else ""))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sylph",    required=True, help="Sylph .profiling.tsv file")
    parser.add_argument("--taxo-dir", required=True, help="Directory containing gtdb_r232_taxonomy")
    parser.add_argument("--output",   required=True, help="Output TSV file")
    args = parser.parse_args()

    taxo = load_taxonomy(args.taxo_dir)
    process_sylph_file(args.sylph, taxo, args.output)


if __name__ == "__main__":
    main()