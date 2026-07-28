# workflow/rules/read_annotate/kofam_rumen.smk

# ── Règle 1 : indexer le catalogue (à faire une seule fois) ──────────────────
rule read_annotate__kofam_rumen_bwa_index:
    """Index le catalogue de gènes rumen pour BWA"""
    input:
        catalogue=features["databases"]["kofam_rumen_fna"],
    output:
        features["databases"]["kofam_rumen_fna"] + ".bwt",
    log:
        PRE_KOFAM_RUMEN / "bwa_index.log",
    container:
        docker["bwa_msamtools"]
    shell:
        "bwa index {input.catalogue} 2> {log}"


# ── Règle 2 : aligner + compter par échantillon ──────────────────────────────
rule read_annotate__kofam_rumen_count:
    """Aligne les reads sur le catalogue rumen et compte par gène"""
    input:
        forwards=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_1.fq.gz",
        reverses=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_2.fq.gz",
        catalogue=features["databases"]["kofam_rumen_fna"],
        index=features["databases"]["kofam_rumen_fna"] + ".bwt",
    output:
        PRE_KOFAM_RUMEN / "{sample_id}.{library_id}.counts.tsv",
    log:
        PRE_KOFAM_RUMEN / "{sample_id}.{library_id}.counts.log",
    container:
        docker["bwa_msamtools"]
    threads: esc("cpus", "read_annotate__kofam_rumen_count")
    resources:
        runtime=esc("runtime", "read_annotate__kofam_rumen_count"),
        mem_mb=esc("mem_mb", "read_annotate__kofam_rumen_count"),
        cpus_per_task=esc("cpus", "read_annotate__kofam_rumen_count"),
        partition=esc("partition", "read_annotate__kofam_rumen_count"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'read_annotate__kofam_rumen_count')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("read_annotate__kofam_rumen_count"))
    shell:
        """
        echo "Starting kofam_rumen_count attempt {resources.attempt}" > {log}.{resources.attempt}

        bwa mem -t {threads} {input.catalogue} {input.forwards} {input.reverses} 2>> {log}.{resources.attempt} \
            | samtools sort -@ {threads} \
            | msamtools filter -b -l 80 -p 95 -z 5 - \
            | msamtools count --multi prop - > {output} 2>> {log}.{resources.attempt}

        mv {log}.{resources.attempt} {log}
        """


# ── Règle 3 : merger KO par échantillon ──────────────────────────────────────
rule read_annotate__kofam_rumen_merge_ko:
    """Croise les comptages avec l'annotation KofamScan → table KO x échantillon"""
    input:
        counts=PRE_KOFAM_RUMEN / "{sample_id}.{library_id}.counts.tsv",
        ko_annot=features["databases"]["kofam_rumen_annotation"],  # annotation_kegg.tsv
    output:
        PRE_KOFAM_RUMEN / "{sample_id}.{library_id}.ko_counts.tsv",
    log:
        PRE_KOFAM_RUMEN / "{sample_id}.{library_id}.ko_counts.log",
    shell:
        """
        # Joint sur le nom de gène (colonne 1 des deux fichiers)
        awk 'NR==FNR {{ko[$1]=$2; next}} ($1 in ko) {{print ko[$1], $0}}' \
            {input.ko_annot} {input.counts} \
            > {output} 2> {log}
        """


# ── Règle cible ──────────────────────────────────────────────────────────────
rule read_annotate__kofam_rumen:
    """Run KO annotation rumen custom"""
    input:
        [
            PRE_KOFAM_RUMEN / f"{sample_id}.{library_id}.ko_counts.tsv"
            for sample_id, library_id in SAMPLE_LIBRARY
        ]