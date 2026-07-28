#!/usr/bin/env python3
"""
sylph_add_taxonomy.py
---------------------
Ajoute la taxonomie GTDB aux résultats Sylph en cherchant dans les 4 fichiers
de taxonomie de la BDD rumen custom. Produit un fichier par sample.

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


def clean_genome_name(path_or_name: str) -> str:
    """Extrait le nom de base sans extension connue (.fa, .fna, .fasta, .ok)."""
    name = Path(path_or_name).name
    name = re.sub(r'\.(fasta|fna|fa)$', '', name, flags=re.IGNORECASE)
    name = re.sub(r'\.ok$', '', name)
    return name


def load_taxonomy(taxo_dir: str) -> dict:
    """
    Retourne un dict {genome_name: taxonomy_string} en fusionnant les 4 sources.
    Priorité : SupaCow > RBAG > Hungate > RCG
    """
    taxo = {}
    taxo_dir = Path(taxo_dir)

    # ── 1. SupaCowTaxTotal.tsv ──
    supacow = taxo_dir / "SupaCowTaxTotal.tsv"
    if supacow.exists():
        with open(supacow) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                name = row["sample"].strip()
                levels = []
                for lvl, prefix in [
                    ("Kingdom", "d__"), ("Phylum", "p__"), ("Class", "c__"),
                    ("Order", "o__"), ("Family", "f__"), ("Genus", "g__"), ("Species", "s__")
                ]:
                    val = row.get(lvl, "").strip()
                    if val and val != "NA":
                        if prefix == "s__":
                            genus = row.get("Genus", "").strip()
                            val = f"{genus} {val}" if genus and genus != "NA" else val
                        levels.append(f"{prefix}{val}")
                if levels:
                    taxo[name] = ";".join(levels)

    # ── 2. RBAG_gtdbtk_summary.tsv ──
    rbag = taxo_dir / "RBAG_gtdbtk_summary.tsv"
    if rbag.exists():
        with open(rbag) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                name = clean_genome_name(row["user_genome"].strip())
                taxo_str = row["classification"].strip()
                levels = [
                    lvl for lvl in taxo_str.split(";")
                    if lvl and "__" in lvl and lvl.split("__")[1] != ""
                ]
                if levels and name not in taxo:
                    taxo[name] = ";".join(levels)

    # ── 3. Hungate (40168_*.tsv) ──
    hungate_files = list(taxo_dir.glob("40168_*.tsv"))
    if hungate_files:
        with open(hungate_files[0]) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                name = row["Original_Genome_Name"].strip()
                taxo_str = row["gtdb_taxanomy"].strip()
                levels = [
                    lvl for lvl in taxo_str.split(";")
                    if lvl and "__" in lvl and lvl.split("__")[1] != ""
                ]
                if levels and name not in taxo:
                    taxo[name] = ";".join(levels)

    # ── 4. RCG_450_summary.tsv ──
    rcg = taxo_dir / "RCG_450_summary.tsv"
    if rcg.exists():
        with open(rcg) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                name = row["Genome name"].strip()
                taxo_str = row["Taxonomy"].strip()
                levels = [
                    lvl for lvl in taxo_str.split(";")
                    if lvl and "__" in lvl and lvl.split("__")[1] != ""
                ]
                if levels and name not in taxo:
                    taxo[name] = ";".join(levels)

    return taxo


def process_sylph_file(sylph_file: str, taxo: dict, output: str):
    """Joint un fichier Sylph avec la taxonomie et écrit le fichier de sortie."""
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

    print(f"[sylph_add_taxonomy] {len(out_rows)} lignes → {output}")
    if not_found:
        print(f"[sylph_add_taxonomy] ⚠️  {len(not_found)} génomes sans taxonomie : "
              + ", ".join(sorted(not_found)[:5])
              + ("..." if len(not_found) > 5 else ""))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sylph",    required=True, help="Fichier .profiling.tsv Sylph")
    parser.add_argument("--taxo-dir", required=True, help="Dossier des 4 fichiers taxonomie")
    parser.add_argument("--output",   required=True, help="Fichier de sortie TSV")
    args = parser.parse_args()

    taxo = load_taxonomy(args.taxo_dir)
    process_sylph_file(args.sylph, taxo, args.output)


if __name__ == "__main__":
    main()