#!/usr/bin/env python3
"""
MAG and Contig Formatter for MetaX Compatibility
-------------------------------------------------
Processes MAG and contig annotation outputs from multiple tools
(EggNOG, DRAM, Bakta), merges FASTA sequences, standardizes headers,
maps functional annotations, and computes consensus taxonomy per contig.

Outputs are organized in the specified output directory (--outdir).

Usage:
    python metax_formatter.py --sample SAMPLE_ID --outdir /path/to/output --taxdb /path/to/taxdb

"""

import os
import sys
import hashlib
import argparse
from collections import defaultdict, Counter
from Bio import SeqIO
from ete3 import NCBITaxa


# ─────────────────────────────────────────────
# Argument parsing
# ─────────────────────────────────────────────

parser = argparse.ArgumentParser(
    description="MAG and Contig Formatter for MetaX Compatibility",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__
)
parser.add_argument(
    "--sample", required=True,
    help="Sample ID (e.g. EGP016_25_B2_004)"
)
parser.add_argument(
    "--outdir", required=True,
    help="Output directory for MetaX files"
)
parser.add_argument(
    "--taxdb", required=True,
    help="Path to directory containing the ete3 taxa.sqlite database"
)
args = parser.parse_args()

sample     = args.sample
output_dir = args.outdir
taxdb_path = args.taxdb

os.makedirs(output_dir, exist_ok=True)


# ─────────────────────────────────────────────
# Initialize NCBI taxonomy and caching
# ─────────────────────────────────────────────

ncbi = NCBITaxa(dbfile=os.path.join(taxdb_path, "taxa.sqlite"))

taxid_rank_cache: dict = {}
taxid_cache: dict      = {}

# ─────────────────────────────────────────────
# Helper functions
# ─────────────────────────────────────────────

def seq_hash(seq: str) -> str:
    """Return MD5 hash of a sequence string."""
    return hashlib.md5(seq.encode()).hexdigest()


def fill_uncharacterized(taxonomy_str: str) -> str:
    """
    Replace empty taxonomic ranks with 'uncharacterized_<previous_rank>'.
    E.g. 'p__Firmicutes;c__;o__' → 'p__Firmicutes;c__uncharacterized_Firmicutes;o__uncharacterized_uncharacterized_Firmicutes'
    """
    if not taxonomy_str or taxonomy_str == "unknown taxonomy":
        return taxonomy_str

    ranks  = taxonomy_str.split(";")
    filled = []
    for r in ranks:
        if "__" in r:
            prefix, name = r.split("__", 1)
            if not name.strip():
                prev_name = filled[-1].split("__")[1] if filled else "unknown"
                name = f"uncharacterized_{prev_name}"
            filled.append(f"{prefix}__{name}")
        else:
            filled.append(r)
    return ";".join(filled)


def get_rank_cached(ncbi_obj, taxid) -> str:
    """Get NCBI rank for a taxid using local cache; returns 'no_rank' if missing."""
    taxid = str(taxid)
    if taxid in taxid_rank_cache:
        return taxid_rank_cache[taxid]
    try:
        rank = ncbi_obj.get_rank([int(taxid)]).get(int(taxid), "no_rank")
    except Exception:
        rank = "no_rank"
    taxid_rank_cache[taxid] = rank
    return rank


def get_lineage_cached(ncbi_obj, taxid):
    """Return full lineage for a taxid as list of (taxid, name), using local cache."""
    taxid = int(taxid)
    if taxid in taxid_cache:
        return taxid_cache[taxid]
    try:
        lineage      = ncbi_obj.get_lineage(taxid)
        names        = ncbi_obj.get_taxid_translator(lineage)
        lineage_info = [(t, names[t]) for t in lineage]
    except Exception:
        lineage_info = [(0, "unknown")]
    taxid_cache[taxid] = lineage_info
    return lineage_info


def get_ncbi_taxonomy(taxid):
    """Return full NCBI taxonomy path for a taxid as list of (taxid, name)."""
    try:
        path_info = get_lineage_cached(ncbi, taxid)
        if not path_info or all(t == 0 for t, _ in path_info):
            return [(0, "unknown")]
        return path_info
    except Exception:
        return [(0, "unknown")]


# ─────────────────────────────────────────────
# Taxonomy constants
# ─────────────────────────────────────────────

RANK_PREFIX = {
    "domain":  "d__", "kingdom": "k__", "phylum":  "p__", "class":   "c__",
    "order":   "o__", "family":  "f__", "genus":   "g__", "species": "s__",
}
STANDARD_RANKS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]


def consensus_taxonomy(entries):
    """
    Compute bottom-up consensus taxonomy for a contig.

    Parameters
    ----------
    entries : list of (taxonomy_path, max_annot_lvl)
        taxonomy_path is a list of (taxid, name) tuples.

    Returns
    -------
    (consensus_str, num_proteins_matching, total_proteins_with_path)
    """
    cleaned_paths = []
    for taxonomy_path, _ in entries:
        if not taxonomy_path or taxonomy_path == ["unknown"]:
            continue
        filtered = [
            (tid, name) for tid, name in taxonomy_path
            if get_rank_cached(ncbi, tid) in STANDARD_RANKS
        ]
        if filtered:
            cleaned_paths.append(filtered)

    total_prots = len(cleaned_paths)
    if total_prots == 0:
        return "unknown taxonomy", 0, 0

    consensus = []
    level = 0
    while True:
        counter = Counter()
        for path in cleaned_paths:
            if len(path) > level:
                counter[path[level]] += 1
        if not counter:
            break
        taxon, count = counter.most_common(1)[0]
        if count / total_prots >= 0.5:
            prefix = RANK_PREFIX.get(get_rank_cached(ncbi, taxon[0]), "")
            consensus.append(f"{prefix}{taxon[1]}({taxon[0]})")
            level += 1
        else:
            break

    good = sum(
        1 for path in cleaned_paths
        if len(path) >= len(consensus)
        and [p[1] for p in path[:len(consensus)]] == [
            c.split("__")[1].split("(")[0] for c in consensus
        ]
    )
    return "|".join(consensus) if consensus else "unknown taxonomy", good, total_prots


# ─────────────────────────────────────────────
# Input paths
# ─────────────────────────────────────────────

contig_fasta             = f"results/contig_annotate/prodigal/{sample}/{sample}.prodigal.fa"
eggnog_contig_annotation = f"results/contig_annotate/eggnog/{sample}/eggnog_output.emapper.annotations"

# EggNOG on MAGs: now per-sample (all raw bins), not the old dereplicated global run
eggnog_mag_annotation    = f"results/mag_annotate/eggnog_mags/{sample}/eggnog.emapper.annotations"

# GTDB-Tk: now run on ALL raw bins (not dereplicated), both bacteria and archaea
gtdbtk_bac               = "results/mag_annotate/gtdbtk_allbins/gtdbtk.bac120.summary.tsv"
gtdbtk_ar                = "results/mag_annotate/gtdbtk_allbins/gtdbtk.ar53.summary.tsv"

fasta_paths = {
    "eggnog": f"results/mag_annotate/eggnog_mags/{sample}/eggnog.emapper.genepred.fasta",
    "dram":   f"results/mag_annotate/dram_mags/{sample}/annotate/genes.faa",
    "bakta":  f"results/mag_annotate/bakta_mags/bakta_{sample}_fixed_headers.faa",
}

dram_annot_file = f"results/mag_annotate/dram_mags/{sample}/annotate/annotations.tsv"
bakta_tsv_file  = f"results/mag_annotate/bakta_mags/bakta_{sample}_faa_fixed.tsv"


# ─────────────────────────────────────────────
# 1. Combine MAG FASTA sequences
# ─────────────────────────────────────────────

combined_mag_out = os.path.join(output_dir, f"{sample}_Fasta_MAGs.faa")
seq_info: dict   = {}   # key = (seq_hash, bin_part, contig_name)

for tool in ["eggnog", "dram", "bakta"]:
    fasta_file = fasta_paths.get(tool)
    if not fasta_file or not os.path.exists(fasta_file):
        print(f"⚠️  Warning: FASTA file missing for tool '{tool}' — {fasta_file}", file=sys.stderr)
        continue

    for rec in SeqIO.parse(fasta_file, "fasta"):
        seq    = str(rec.seq).replace("*", "").strip()
        h      = seq_hash(seq)
        rec_id = rec.id

        # Remove accidental double-prefix (e.g. SAMPLE_SAMPLE_...)
        double_prefix = f"{sample}_{sample}"
        if rec_id.startswith(double_prefix):
            rec_id = rec_id.replace(double_prefix, sample, 1)

        # Skip sequences from other samples
        if ":" in rec_id:
            rec_sample = rec_id.split(":")[0]
            if rec_sample != sample:
                continue

        # Extract bin and contig info
        if ":" in rec_id and "@" in rec_id:
            parts       = rec_id.split(":")
            bin_part    = parts[1].split("@")[0].replace("bin_", "bin")
            contig_full = parts[1].split("@")[1]
            contig_name, prot_num = (
                contig_full.rsplit("_", 1) if "_" in contig_full else (contig_full, "1")
            )
        else:
            bin_part    = "noBin"
            contig_name = rec_id
            prot_num    = "1"

        contig_name = contig_name.replace("contig_", "contig-")
        key = (h, bin_part, contig_name)

        if key not in seq_info:
            seq_info[key] = {"tools": set(), "prot_num": prot_num, "seq": seq, "orig_ids": []}
        seq_info[key]["tools"].add(tool)
        seq_info[key]["orig_ids"].append(rec.id)

with open(combined_mag_out, "w") as fout:
    for (h, bin_part, contig_name), info in seq_info.items():
        tools_str  = "-".join(sorted(info["tools"]))
        new_header = (
            f"{sample.replace('_', '-')}-{bin_part}-{contig_name}"
            f"_{tools_str}-Protein{info['prot_num']}"
        )
        fout.write(f">{new_header}\n{info['seq']}\n")

print(f"✅  MAG FASTA written → {combined_mag_out}")


# ─────────────────────────────────────────────
# 2. Combine contig FASTA sequences
# ─────────────────────────────────────────────

combined_contig_out = os.path.join(output_dir, f"{sample}_Fasta_Contigs.faa")

if os.path.exists(contig_fasta):
    with open(combined_contig_out, "w") as fout:
        for rec in SeqIO.parse(contig_fasta, "fasta"):
            seq     = str(rec.seq).replace("*", "").strip()
            prot_id = rec.id

            if "@" in prot_id:
                contig_full = prot_id.split("@")[1]
                contig_name, prot_num = (
                    contig_full.rsplit("_", 1) if "_" in contig_full else (contig_full, "1")
                )
            else:
                contig_name, prot_num = prot_id, "1"

            contig_name = contig_name.replace("contig_", "contig-")
            header = f"{sample.replace('_', '-')}-binNA-{contig_name}_Protein{prot_num}"
            fout.write(f">{header}\n{seq}\n")

    print(f"✅  Contig FASTA written → {combined_contig_out}")
else:
    print(f"⚠️  Warning: contig FASTA not found — {contig_fasta}", file=sys.stderr)


# ─────────────────────────────────────────────
# 3. Load DRAM annotations (full)
# ─────────────────────────────────────────────

dram_annotations: dict = {}
dram_seq_by_prot: dict = {}
dram_header: list      = []

if os.path.exists(dram_annot_file):
    with open(dram_annot_file) as fin:
        dram_header = fin.readline().strip().split("\t")
        for line in fin:
            fields  = line.strip().split("\t")
            prot_id = fields[0]
            dram_annotations[prot_id] = dict(zip(dram_header, fields))
else:
    print(f"⚠️  Warning: DRAM annotation file not found — {dram_annot_file}", file=sys.stderr)

if os.path.exists(fasta_paths["dram"]):
    for rec in SeqIO.parse(fasta_paths["dram"], "fasta"):
        dram_seq_by_prot[rec.id] = str(rec.seq).replace("*", "").strip()

seq_to_dram = {seq_hash(seq): pid for pid, seq in dram_seq_by_prot.items()}


# ─────────────────────────────────────────────
# 4. Load Bakta annotations
# ─────────────────────────────────────────────

bakta_annotations: dict = {}
bakta_seq_by_prot: dict = {}
bakta_header: list      = []

if os.path.exists(bakta_tsv_file):
    with open(bakta_tsv_file) as fin:
        header_line   = fin.readline().strip()
        bakta_header  = header_line.split("\t")[1:]   # skip "Sequence Id" column
        for line in fin:
            fields  = line.strip().split("\t")
            prot_id = fields[0]
            bakta_annotations[prot_id] = dict(zip(bakta_header, fields[1:]))
else:
    print(f"⚠️  Warning: Bakta annotation TSV not found — {bakta_tsv_file}", file=sys.stderr)

BAKTA_EXCLUDE = {"Locus Tag", "Type", "Start", "Stop", "Strand"}

bakta_cols_filtered = [
    c for c in bakta_header
    if not any(ex.lower() in c.lower() for ex in BAKTA_EXCLUDE)
]

if os.path.exists(fasta_paths["bakta"]):
    for rec in SeqIO.parse(fasta_paths["bakta"], "fasta"):
        bakta_seq_by_prot[rec.id] = str(rec.seq).replace("*", "").strip()

seq_to_bakta = {seq_hash(seq): pid for pid, seq in bakta_seq_by_prot.items()}


# ─────────────────────────────────────────────
# 5. MAG functional annotations (EggNOG + DRAM + Bakta)
# ─────────────────────────────────────────────

function_mag_out    = os.path.join(output_dir, f"{sample}_Function_MAGs.tsv")
eggnog_annotations: dict  = {}
eggnog_seq_by_prot: dict  = {}
eggnog_header: list       = []

# Load EggNOG MAG annotations.
# eggnog_mag_annotation is now per-sample (eggnog_mags/{sample}/...), so it already
# only contains this sample's proteins. The startswith filter is kept as a harmless
# safety net (all lines already start with "{sample}:").
if os.path.exists(eggnog_mag_annotation):
    with open(eggnog_mag_annotation) as fin:
        for line in fin:
            if line.startswith("#query"):
                eggnog_header = line.strip().split("\t")
                continue
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if not fields[0].startswith(sample + ":"):
                continue
            eggnog_annotations[fields[0]] = fields
else:
    print(f"⚠️  Warning: EggNOG MAG annotation file not found — {eggnog_mag_annotation}", file=sys.stderr)

if os.path.exists(fasta_paths["eggnog"]):
    for rec in SeqIO.parse(fasta_paths["eggnog"], "fasta"):
        if not rec.id.startswith(sample + ":"):
            continue
        eggnog_seq_by_prot[rec.id] = str(rec.seq).replace("*", "").strip()

seq_to_eggnog = {seq_hash(seq): pid for pid, seq in eggnog_seq_by_prot.items()}

# DRAM columns to exclude from the merged output
DRAM_EXCLUDE = {"fasta", "scaffold", "gene_position", "start_position",
                "end_position", "strandedness"}
dram_cols_filtered = [c for c in dram_header if not any(ex in c.lower() for ex in DRAM_EXCLUDE)]

# Fallback column count for empty annotation rows
_eggnog_ncols = len(next(iter(eggnog_annotations.values()))) if eggnog_annotations else 2

with open(function_mag_out, "w") as fout:
    # Header
    fout.write("\t".join(
        ["query"] +
        eggnog_header[1:] +
        [f"DRAM_{c}" for c in dram_cols_filtered] +
        [f"BAKTA_{c}" for c in bakta_cols_filtered]
    ) + "\n")

    for rec in SeqIO.parse(combined_mag_out, "fasta"):
        seq    = str(rec.seq).replace("*", "").strip()
        header = rec.id
        h      = seq_hash(seq)

        # EggNOG
        eggnog_prot = seq_to_eggnog.get(h)
        eggnog_data = (
            eggnog_annotations[eggnog_prot][1:]
            if eggnog_prot and eggnog_prot in eggnog_annotations
            else ["-"] * (_eggnog_ncols - 1)
        )

        # DRAM
        dram_prot = seq_to_dram.get(h)
        dram_values = (
            [dram_annotations[dram_prot].get(c, "-") for c in dram_cols_filtered]
            if dram_prot and dram_prot in dram_annotations
            else ["-"] * len(dram_cols_filtered)
        )

        # Bakta
        bakta_prot = seq_to_bakta.get(h)
        bakta_values = (
            [bakta_annotations[bakta_prot].get(c, "-") for c in bakta_cols_filtered]
            if bakta_prot and bakta_prot in bakta_annotations
            else ["-"] * len(bakta_cols_filtered)
        )

        fout.write("\t".join([header] + eggnog_data + dram_values + bakta_values) + "\n")

print(f"✅  MAG function table written → {function_mag_out}")


# ─────────────────────────────────────────────
# 6. MAG taxonomy from GTDB-Tk (all-bins direct: bac120 + ar53, no propagation)
# ─────────────────────────────────────────────

taxonomy_mag_out = os.path.join(output_dir, f"{sample}_Taxonomy_MAGs.tsv")

# Load direct per-bin taxonomy for THIS sample, from both bacterial and archaeal summaries.
# GTDB-Tk here was run on ALL raw bins (not just dereplicated representatives),
# so every bin can have its own direct classification -- no dRep propagation needed.
direct_taxonomy: dict = {}   # key = MetaX bin id 'SAMPLE-binXXXXXX'
for gtdb_file in (gtdbtk_bac, gtdbtk_ar):
    if os.path.exists(gtdb_file):
        with open(gtdb_file) as fin:
            for line in fin:
                if line.startswith("user_genome") or not line.strip():
                    continue
                if not line.startswith(sample):
                    continue
                fields = line.strip().split("\t")
                # user_genome e.g. 'EGP016_25_B2_002:bin_000061' -> 'EGP016-25-B2-002-bin000061'
                genome = fields[0].replace(":", "-").replace("_", "-").replace("bin-", "bin")
                direct_taxonomy[genome] = (
                    fields[1] if len(fields) > 1 else "unknown taxonomy"
                )
    else:
        print(f"⚠️  Warning: GTDB-Tk file not found — {gtdb_file}", file=sys.stderr)

# Write taxonomy for every bin present in the combined MAG FASTA, so that
# Taxonomy_MAGs.tsv always matches Fasta_MAGs.faa / Function_MAGs.tsv exactly.
sample_dash = sample.replace("_", "-")
sample_bins: set = set()
if os.path.exists(combined_mag_out):
    for rec in SeqIO.parse(combined_mag_out, "fasta"):
        sample_bins.add(rec.id.split("-contig")[0])   # 'SAMPLE-binXXXXXX'
else:
    print(f"⚠️  Warning: combined MAG FASTA not found — {combined_mag_out}", file=sys.stderr)

n_ok = n_unknown = 0
with open(taxonomy_mag_out, "w") as fout:
    fout.write("Genome\tTaxonomy\n")
    for metax_bin in sorted(sample_bins):
        taxonomy = direct_taxonomy.get(metax_bin)
        if taxonomy:
            n_ok += 1
        else:
            taxonomy = "unknown taxonomy"
            n_unknown += 1
        fout.write(f"{metax_bin}\t{fill_uncharacterized(taxonomy)}\n")

print(f"✅  MAG taxonomy written → {taxonomy_mag_out} "
      f"({n_ok} classified, {n_unknown} unknown)")



# ─────────────────────────────────────────────
# 7. Contig functional annotations (EggNOG)
# ─────────────────────────────────────────────

function_contig_out = os.path.join(output_dir, f"{sample}_Function_Contigs.tsv")
contig_to_prots: defaultdict = defaultdict(list)
prot_to_fields: dict         = {}
header_line: str             = ""

if os.path.exists(eggnog_contig_annotation):
    with open(eggnog_contig_annotation) as fin:
        for line in fin:
            if line.startswith("#query"):
                header_line = line
                continue
            if line.startswith("#") or not line.strip():
                continue
            fields  = line.strip().split("\t")
            prot_id = fields[0]
            contig_to_prots[prot_id].append(prot_id)
            prot_to_fields[prot_id] = fields
else:
    print(f"⚠️  Warning: EggNOG contig annotation not found — {eggnog_contig_annotation}", file=sys.stderr)

# Map original protein IDs → standardized headers
header_mapping_contigs: dict = {}
if os.path.exists(contig_fasta):
    for rec in SeqIO.parse(contig_fasta, "fasta"):
        prot_id = rec.id
        if "@" in prot_id:
            contig_full = prot_id.split("@")[1]
            contig_name, prot_num = (
                contig_full.rsplit("_", 1) if "_" in contig_full else (contig_full, "1")
            )
        else:
            contig_name, prot_num = prot_id, "1"
        contig_name = contig_name.replace("contig_", "contig-")
        new_id      = f"{sample.replace('_', '-')}-binNA-{contig_name}_Protein{prot_num}"
        header_mapping_contigs[prot_id] = new_id

with open(function_contig_out, "w") as fout:
    if header_line:
        fout.write(header_line)
        ncols = len(header_line.strip().split("\t"))
    else:
        ncols = 2   # fallback

    if os.path.exists(contig_fasta):
        for rec in SeqIO.parse(contig_fasta, "fasta"):
            prot_id = rec.id
            new_id  = header_mapping_contigs.get(prot_id, prot_id)
            if prot_id in prot_to_fields:
                fout.write("\t".join([new_id] + prot_to_fields[prot_id][1:]) + "\n")
            else:
                fout.write(
                    new_id + "\tunknown function\t" + "\t".join(["-"] * (ncols - 2)) + "\n"
                )

print(f"✅  Contig function table written → {function_contig_out}")


# ─────────────────────────────────────────────
# 8. Consensus taxonomy for contigs
# ─────────────────────────────────────────────

taxonomy_contig_out = os.path.join(output_dir, f"{sample}_Taxonomy_Contigs.tsv")
all_contigs: set                = set()
contig_tax: defaultdict         = defaultdict(list)

# Build taxonomic entries from EggNOG contig annotations
if os.path.exists(eggnog_contig_annotation):
    with open(eggnog_contig_annotation) as fin:
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue
            fields        = line.strip().split("\t")
            prot_id       = fields[0]
            max_annot_lvl = fields[5]   # column index 5 = max_annot_lvl

            # Resolve taxonomy path
            if max_annot_lvl.endswith("|-"):
                try:
                    taxid_candidate = int(fields[1].split(".")[0])
                    taxonomy_path   = get_ncbi_taxonomy(taxid_candidate)
                except Exception:
                    taxonomy_path = ["unknown"]
            elif "|" in max_annot_lvl:
                taxid_str, _ = max_annot_lvl.split("|", 1)
                try:
                    taxonomy_path = get_ncbi_taxonomy(int(taxid_str))
                except ValueError:
                    taxonomy_path = ["unknown"]
            else:
                taxonomy_path = ["unknown"]

            # Derive contig ID
            if "@" in prot_id:
                contig_full = prot_id.split("@")[1]
                contig_name = "_".join(contig_full.split("_")[:-1])
            else:
                contig_name = prot_id
            contig_name = contig_name.replace("contig_", "contig-")
            contig_id   = f"{sample.replace('_', '-')}-binNA-{contig_name}"

            contig_tax[contig_id].append((taxonomy_path, max_annot_lvl))

# Ensure every contig appears in the output, even with no annotated proteins
if os.path.exists(contig_fasta):
    for rec in SeqIO.parse(contig_fasta, "fasta"):
        prot_id = rec.id
        if "@" in prot_id:
            contig_full = prot_id.split("@")[1]
            contig_name = "_".join(contig_full.split("_")[:-1])
        else:
            contig_name = prot_id
        contig_name = contig_name.replace("contig_", "contig-")
        contig_id   = f"{sample.replace('_', '-')}-binNA-{contig_name}"
        all_contigs.add(contig_id)
        contig_tax.setdefault(contig_id, [])

with open(taxonomy_contig_out, "w") as fout:
    fout.write("Genome\tConsensus_taxonomy\tProteins_with_consensus\tTotal_proteins_with_max_annotation_level\n")
    for contig in sorted(all_contigs):
        entries              = contig_tax.get(contig, [])
        consensus, good, total = consensus_taxonomy(entries)
        fout.write(f"{contig}\t{consensus}\t{good}\t{total}\n")

print(f"✅  Contig taxonomy written → {taxonomy_contig_out}")
print(f"\n🎉  All outputs written to: {output_dir}")