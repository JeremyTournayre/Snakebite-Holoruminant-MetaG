rule read_annotate__diamond__assign:
    """Run Diamond with SHM/NVME staging, safe for concurrent jobs."""
    input:
        forwards=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_1.fq.gz",
        reverses=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_2.fq.gz",
        database=lambda w: features["databases"]["diamond"][w.diamond_db],
    output:
        out_R1=DIAMOND / "{diamond_db}" / "{sample_id}.{library_id}_R1.out",
        out_R2=DIAMOND / "{diamond_db}" / "{sample_id}.{library_id}_R2.out",
    log:
        DIAMOND / "{diamond_db}_{sample_id}_{library_id}.log",
    benchmark:
        DIAMOND / "benchmark/{diamond_db}_{sample_id}_{library_id}.tsv",
    threads: esc("cpus", "read_annotate__diamond__assign"),
    resources:
        runtime=esc("runtime", "read_annotate__diamond__assign"),
        mem_mb=esc("mem_mb", "read_annotate__diamond__assign"),
        cpus_per_task=esc("cpus", "read_annotate__diamond__assign"),
        slurm_partition=esc("partition", "read_annotate__diamond__assign"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'read_annotate__diamond__assign')['nvme']}",
        attempt=lambda wc, attempt: attempt,
    retries: len(get_escalation_order("read_annotate__diamond__assign")) - 1,
    params:
        max_target_seqs=params["read_annotate"]["diamond"]["max-target-seqs"],
        retries=len(get_escalation_order("read_annotate__diamond__assign")) - 1,
        run_in_shm=config["diamond_shm"],
        copy_dbs=config["copy_dbs"],
        diamond_shm=DIAMONDSHM,
        diamond_nvme=DIAMONDNVME,
        diamond_path=lambda w: os.path.basename(features["databases"]["diamond"][w.diamond_db]),
        out_folder=lambda w: DIAMOND / w.diamond_db,
    container:
        docker["mag_annotate"],
    shell:
        r"""
        set -euo pipefail

        mkdir -p {params.out_folder}

        echo "Starting Diamond attempt {resources.attempt}" > {log}

        DB_SRC="{input.database}"
        DB_FILE="{params.diamond_path}"
        DB_BASE="${{DB_FILE%.dmnd}}"

        SHM="{params.diamond_shm}"
        NVME="{params.diamond_nvme}"

        USE_SHM="{params.run_in_shm}"
        COPY_DBS="{params.copy_dbs}"

        # ===============================
        # Decide staging location
        # ===============================
        DB_DST="$DB_SRC"
        DATA_DIR=""

        if [ "$COPY_DBS" = "True" ]; then
            DB_SIZE=$(du -sb "$DB_SRC" | cut -f1)
            SHM_AVAIL=$(df --output=avail -B1 "$SHM" 2>/dev/null | tail -1 || echo 0)
            NVME_AVAIL=$(df --output=avail -B1 "$NVME" 2>/dev/null | tail -1 || echo 0)

            echo "DB_SIZE=$DB_SIZE SHM_AVAIL=$SHM_AVAIL NVME_AVAIL=$NVME_AVAIL" >> {log}

            if [ "$USE_SHM" = "True" ] && [ "$DB_SIZE" -lt "$SHM_AVAIL" ]; then
                DATA_DIR="$SHM/{wildcards.diamond_db}"
            elif [ "$DB_SIZE" -lt "$NVME_AVAIL" ]; then
                DATA_DIR="$NVME/{wildcards.diamond_db}"
            fi
        fi

        # ===============================
        # Shared DB logic (LOCK / COUNTER)
        # ===============================
        if [ -n "$DATA_DIR" ]; then
            LOCK="$DATA_DIR/.lock"
            DONE="$DATA_DIR/.done"
            COUNTER="$DATA_DIR/.counter"

            mkdir -p "$DATA_DIR"

            # ---- Increment counter ----
            (
                flock -x 200
                COUNT=0
                [ -f "$COUNTER" ] && COUNT=$(cat "$COUNTER")
                COUNT=$((COUNT + 1))
                echo "$COUNT" > "$COUNTER"
                echo "Diamond DB users: $COUNT" >> {log}
            ) 200>"$COUNTER.lock"

            # ---- Copy DB once ----
            if [ ! -f "$DONE" ]; then
                if mkdir "$LOCK" 2>/dev/null; then
                    echo "This job copies Diamond DB to $DATA_DIR" >> {log}
                    cp "$DB_SRC" "$DATA_DIR/"
                    touch "$DONE"
                    rmdir "$LOCK"
                else
                    echo "Waiting for Diamond DB copy..." >> {log}
                    while [ ! -f "$DONE" ]; do
                        sleep 5
                    done
                fi
            fi

            DB_LOC="$DATA_DIR/$DB_BASE"
        else
            DB_LOC="$DB_SRC"
        fi

        echo "Running Diamond with DB_LOC=$DB_LOC" >> {log}

        diamond blastx \
            -d "$DB_LOC" \
            -q {input.forwards} \
            -o {output.out_R1} \
            -p {threads} \
            --max-target-seqs {params.max_target_seqs} >> {log} 2>&1

        diamond blastx \
            -d "$DB_LOC" \
            -q {input.reverses} \
            -o {output.out_R2} \
            -p {threads} \
            --max-target-seqs {params.max_target_seqs} >> {log} 2>&1

        # ===============================
        # Decrement counter + cleanup
        # ===============================
        if [ -n "$DATA_DIR" ]; then
            (
                flock -x 200
                COUNT=$(cat "$COUNTER")
                COUNT=$((COUNT - 1))
                if [ "$COUNT" -le 0 ]; then
                    echo "Last job, cleaning $DATA_DIR" >> {log}
                    rm -rf "$DATA_DIR"
                else
                    echo "$COUNT jobs still using Diamond DB" >> {log}
                    echo "$COUNT" > "$COUNTER"
                fi
            ) 200>"$COUNTER.lock"
        fi
        """


rule read_annotate__diamond__summarise:
    """Run R script to summarise Diamond results"""
    input:
        lambda w: [
            DIAMOND / w.diamond_db / f"{sample_id}.{library_id}_R1.out"
            for sample_id, library_id in SAMPLE_LIBRARY
        ] + [
            DIAMOND / w.diamond_db / f"{sample_id}.{library_id}_R2.out"
            for sample_id, library_id in SAMPLE_LIBRARY
        ]
    output:
        DIAMOND / "{diamond_db}" / "summary_a.tsv",
        DIAMOND / "{diamond_db}" / "summary_b.tsv"
    log:
        DIAMOND / "{diamond_db}" / "read_annotate__diamond__summarise.log"
    threads: esc("cpus", "read_annotate__diamond__summarise")
    resources:
        runtime=esc("runtime", "read_annotate__diamond__summarise"),
        mem_mb=esc("mem_mb", "read_annotate__diamond__summarise"),
        cpus_per_task=esc("cpus", "read_annotate__diamond__summarise"),
        slurm_partition=esc("partition", "read_annotate__diamond__summarise"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'read_annotate__diamond__summarise')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("read_annotate__diamond__summarise"))
    params:
        script=lambda w: config["pipeline_folder"] + "workflow/scripts/downstream_scripts/create_featuretable_diamond_" + w.diamond_db + ".R",
        project_folder=WD,
        database=lambda w: w.diamond_db,
        max_target_seqs=params["read_annotate"]["diamond"]["max-target-seqs"],
    container:
        docker["r_report"]
    shell:
        """
        R -e "project_folder <- '{params.project_folder}'; result_folder <- 'results/read_annotate/diamond/{params.database}'; max_target_seqs <- {params.max_target_seqs}; database <- '{params.database}';  source('{params.script}')" &> {log}
        """

rule read_annotate__diamond:
    """Run all Diamond steps"""
    input:
        [
            DIAMOND / diamond_db / f"{sample_id}.{library_id}_R1.out"
            for sample_id, library_id in SAMPLE_LIBRARY
            for diamond_db in features["databases"]["diamond"]
        ] + [
            DIAMOND / diamond_db / f"{sample_id}.{library_id}_R2.out"
            for sample_id, library_id in SAMPLE_LIBRARY
            for diamond_db in features["databases"]["diamond"]
        ] + [
            DIAMOND / diamond_db / "summary_a.tsv"
            for diamond_db in features["databases"]["diamond"]
        ] + [
            DIAMOND / diamond_db / "summary_b.tsv"
            for diamond_db in features["databases"]["diamond"]
        ]
