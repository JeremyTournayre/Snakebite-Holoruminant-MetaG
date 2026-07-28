rule mag_annotate__eggnog_mags__search:
    """EggNOG gene prediction + seed orthologs, per sample, over all raw bins."""
    input:
        contigs=MAGSCOT / "{assembly_id}.fa.gz",
    output:
        seed           = EGGNOGMAG / "{assembly_id}" / "eggnog.emapper.seed_orthologs",
        genepred_fasta = EGGNOGMAG / "{assembly_id}" / "eggnog.emapper.genepred.fasta",
        genepred_gff   = EGGNOGMAG / "{assembly_id}" / "eggnog.emapper.genepred.gff",
    log:
        EGGNOGMAG / "{assembly_id}" / "eggnog.log",
    container:
        docker["eggnog"]
    params:
        out_dir=lambda wc: f"{EGGNOGMAG}/{wc.assembly_id}",
        db=features["databases"]["eggnog"],
        prefix="eggnog",
    threads: esc("cpus", "mag_annotate__eggnog_mags__search")
    resources:
        runtime=esc("runtime", "mag_annotate__eggnog_mags__search"),
        mem_mb=esc("mem_mb", "mag_annotate__eggnog_mags__search"),
        cpus_per_task=esc("cpus", "mag_annotate__eggnog_mags__search"),
        partition=esc("partition", "mag_annotate__eggnog_mags__search"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'mag_annotate__eggnog_mags__search')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("mag_annotate__eggnog_mags__search"))
    shell:
        """
        if [ ! -s {input.contigs} ]; then
            echo "[INFO] {input.contigs} empty/missing. Skipping." > {log}
            mkdir -p {params.out_dir}
            touch {output.seed} {output.genepred_fasta} {output.genepred_gff}
            exit 0
        fi

        emapper.py -m diamond \
                   --override \
                   --data_dir {params.db} \
                   --itype metagenome \
                   --genepred prodigal \
                   --dbmem \
                   --no_annot \
                   --no_file_comments \
                   --cpu {threads} \
                   -i {input.contigs} \
                   --output_dir {params.out_dir} \
                   -o {params.prefix} \
        2>> {log} 1>&2
        """


checkpoint mag_annotate__eggnog_mags__split:
    """Split each sample's seed_orthologs into chunks."""
    input:
        seed = EGGNOGMAG / "{assembly_id}" / "eggnog.emapper.seed_orthologs",
    output:
        directory(EGGNOGMAG / "{assembly_id}" / "Chunks"),
    log:
        EGGNOGMAG / "{assembly_id}" / "eggnog_split.log",
    params:
        chunk_size = params["mag_annotate"]["eggnog"]["split"],
        prefix     = lambda wc: f"{EGGNOGMAG}/{wc.assembly_id}/Chunks/seed.chunk.",
    threads: 1
    shell:
        """
        mkdir -p {output}
        if [ ! -s {input.seed} ]; then
            echo "[INFO] empty seed for {wildcards.assembly_id}, no chunks." > {log}
            exit 0
        fi
        split -d -l {params.chunk_size} {input.seed} {params.prefix} \
            --additional-suffix=.seed_orthologs 2>> {log}
        """


rule mag_annotate__eggnog_mags__annotate_chunk:
    input:
        seed = EGGNOGMAG / "{assembly_id}" / "Chunks" / "seed.chunk.{i}.seed_orthologs",
    output:
        annot = EGGNOGMAG / "{assembly_id}" / "Chunks" / "seed.chunk.{i}.emapper.annotations",
        done  = EGGNOGMAG / "{assembly_id}" / "Chunks" / "seed.chunk.{i}.emapper.annotations.done",
    log:
        EGGNOGMAG / "{assembly_id}" / "Chunks" / "seed.chunk.{i}.annotations.log",
    params:
        fa     = features["databases"]["eggnog"],
        out    = lambda wc: f"seed.chunk.{wc.i}",
        outdir = lambda wc: f"{EGGNOGMAG}/{wc.assembly_id}/Chunks",
    threads: esc("cpus", "mag_annotate__eggnog_mags__annotate_chunk")
    resources:
        runtime=esc("runtime", "mag_annotate__eggnog_mags__annotate_chunk"),
        mem_mb=esc("mem_mb", "mag_annotate__eggnog_mags__annotate_chunk"),
        cpus_per_task=esc("cpus", "mag_annotate__eggnog_mags__annotate_chunk"),
        partition=esc("partition", "mag_annotate__eggnog_mags__annotate_chunk"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'mag_annotate__eggnog_mags__annotate_chunk')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("mag_annotate__eggnog_mags__annotate_chunk"))
    container:
        docker["eggnog"]
    shell:
        """
        DATA_DIR="/dev/shm/eggnog_data_mags"
        LOCK_FILE="$DATA_DIR/.lock"
        DONE_FILE="$DATA_DIR/.done"
        COUNTER_FILE="$DATA_DIR/.counter"
        mkdir -p $DATA_DIR

        ( flock -x 200
          COUNT=0; [ -f "$COUNTER_FILE" ] && COUNT=$(cat "$COUNTER_FILE")
          COUNT=$((COUNT+1)); echo $COUNT > "$COUNTER_FILE"
        ) 200>"$COUNTER_FILE.lock"

        if [ ! -f "$DONE_FILE" ]; then
            if mkdir "$LOCK_FILE" 2>/dev/null; then
                cp -r {params.fa}/* $DATA_DIR/
                touch "$DONE_FILE"; rmdir "$LOCK_FILE"
            else
                while [ ! -f "$DONE_FILE" ]; do sleep 2; done
            fi
        fi

        emapper.py --data_dir $DATA_DIR \
                   --annotate_hits_table {input.seed} \
                   --no_file_comments \
                   -o {params.out} \
                   --output_dir {params.outdir} \
                   --override \
                   --cpu {threads} >> {log} 2>&1

        ( flock -x 200
          COUNT=$(cat "$COUNTER_FILE"); COUNT=$((COUNT-1))
          if [ $COUNT -le 0 ]; then echo 0 > "$COUNTER_FILE"; rm -rf "$DATA_DIR"; rm -f "$DONE_FILE"
          else echo $COUNT > "$COUNTER_FILE"; fi
        ) 200>"$COUNTER_FILE.lock"

        touch {output.done}
        """


rule mag_annotate__eggnog_mags__merge:
    input:
        mag_annotate_aggregate_annotations_eggnog_mags,
    output:
        EGGNOGMAG / "{assembly_id}" / "eggnog.emapper.annotations",
    log:
        EGGNOGMAG / "{assembly_id}" / "merge_annotations.log",
    shell:
        """
        if [ -z "{input}" ]; then
            echo "[INFO] no annotation chunks for {wildcards.assembly_id}" > {log}
            touch {output}
        else
            cat {input} > {output} 2> {log}
        fi
        """


rule mag_annotate__eggnog_mags:
    input:
        expand(EGGNOGMAG / "{assembly_id}" / "eggnog.emapper.annotations", assembly_id=ASSEMBLIES),
        expand(EGGNOGMAG / "{assembly_id}" / "eggnog.emapper.genepred.fasta", assembly_id=ASSEMBLIES),