rule read_annotate__sylph_profile:
    """Run Sylph profiler"""
    input:
        forwards=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_1.fq.gz",
        reverses=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_2.fq.gz",
        db=features["databases"]["sylph"],
    output:
        PRE_SYLPH / "{sample_id}.{library_id}.profiling.tsv",
    log:
        PRE_SYLPH / "{sample_id}.{library_id}.profiling_report.log",
    container:
        docker["sylph"]
    threads: esc("cpus", "read_annotate__sylph_profile")
    resources:
        runtime=esc("runtime", "read_annotate__sylph_profile"),
        mem_mb=esc("mem_mb", "read_annotate__sylph_profile"),
        cpus_per_task=esc("cpus", "read_annotate__sylph_profile"),
        partition=esc("partition", "read_annotate__sylph_profile"),
        #gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'read_annotate__sylph_profile')['nvme']}",
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'read_annotate__sylph_profile')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("read_annotate__sylph_profile"))
    shell:
        """
        
        echo "Starting Sylph rule attempt {resources.attempt}" > {log}.{resources.attempt}
        
        # Below this line is a dev part for testing the resubmission and escalation 
        # if [ {resources.attempt} -le 2 ]; then
        #    echo "Failing intentionally on attempt {resources.attempt}" >> {log}.{resources.attempt}
        #    ls --nonexistentflag
        #fi
        
        sylph profile {input.db} -1 {input.forwards} -2 {input.reverses} -t {threads} > {output} 2>> {log}.{resources.attempt}
        
        mv {log}.{resources.attempt} {log}
        """

SYLPH_IS_RUMEN_CUSTOM = "rumen_db" in str(features["databases"]["sylph"])

rule read_annotate__sylph:
    """Run Sylph"""
    input:
        [
            (
                PRE_SYLPH / f"{sample_id}.{library_id}.profiling_with_taxonomy.tsv"
                if SYLPH_IS_RUMEN_CUSTOM else
                PRE_SYLPH / f"{sample_id}.{library_id}.profiling.tsv"
            )
            for sample_id, library_id in SAMPLE_LIBRARY
        ]


# ── Règle : ajout taxonomie rumen aux résultats Sylph ────────────────────────
if SYLPH_IS_RUMEN_CUSTOM:

    rule read_annotate__sylph_add_taxonomy:
        """Ajoute la taxonomie rumen custom aux résultats Sylph (un fichier par sample)"""
        input:
            profiling=PRE_SYLPH / "{sample_id}.{library_id}.profiling.tsv",
        output:
            PRE_SYLPH / "{sample_id}.{library_id}.profiling_with_taxonomy.tsv",
        log:
            PRE_SYLPH / "{sample_id}.{library_id}.add_taxonomy.log",
        params:
            taxo_dir=str(Path(features["databases"]["sylph"]).parent.parent / "DB" / "taxonomy"),
            folder=config["pipeline_folder"],
        shell:
            """
            python {params.folder}/workflow/scripts/sylph_add_taxonomy.py \
                --sylph {input.profiling} \
                --taxo-dir {params.taxo_dir} \
                --output {output} \
                > {log} 2>&1
            """

    rule read_annotate__sylph_taxonomy:
        """Collect all Sylph taxonomy-annotated profiles"""
        input:
            [
                PRE_SYLPH / f"{sample_id}.{library_id}.profiling_with_taxonomy.tsv"
                for sample_id, library_id in SAMPLE_LIBRARY
            ]
