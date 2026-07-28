from glob import glob
from os.path import join

rule mag_annotate__eggnog_annotations_chunk:
    input:
        seed = EGGNOG / "Chunks/seed.chunk.{i}.seed_orthologs"
    output:
        EGGNOG / "Chunks/seed.chunk.{i}.emapper.annotations",
        done = EGGNOG / "Chunks/seed.chunk.{i}.emapper.annotations.done"
    log:
        EGGNOG / "Chunks/seed.chunk.{i}.annotations.log"
    params:
        fa = features["databases"]["eggnog"],
        out = lambda wc: f"seed.chunk.{wc.i}",
        outdir = EGGNOG / "Chunks"
    threads: esc("cpus", "mag_annotate__eggnog_annotations_chunk")
    resources:
        runtime=esc("runtime", "mag_annotate__eggnog_annotations_chunk"),
        mem_mb=esc("mem_mb", "mag_annotate__eggnog_annotations_chunk"),
        cpus_per_task=esc("cpus", "mag_annotate__eggnog_annotations_chunk"),
        partition=esc("partition", "mag_annotate__eggnog_annotations_chunk"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'mag_annotate__eggnog_annotations_chunk')['nvme']}",
        attempt=get_attempt,  
    retries: len(get_escalation_order("mag_annotate__eggnog_annotations_chunk"))
    container:
        docker["eggnog"]
    shell:
        """
        DATA_DIR="/dev/shm/eggnog_data_mags"
        LOCK_FILE="$DATA_DIR/.lock"
        DONE_FILE="$DATA_DIR/.done"
        COUNTER_FILE="$DATA_DIR/.counter"

        mkdir -p $DATA_DIR

        # increment counter
        (
            flock -x 200
            COUNT=0
            [ -f "$COUNTER_FILE" ] && COUNT=$(cat "$COUNTER_FILE")
            COUNT=$((COUNT+1))
            echo $COUNT > "$COUNTER_FILE"
        ) 200>"$COUNTER_FILE.lock"

        # copy DB once
        if [ ! -f "$DONE_FILE" ]; then
            if mkdir "$LOCK_FILE" 2>/dev/null; then
                cp -r {params.fa}/* $DATA_DIR/
                touch "$DONE_FILE"
                rmdir "$LOCK_FILE"
            else
                while [ ! -f "$DONE_FILE" ]; do sleep 2; done
            fi
        fi

        # run annotation
        emapper.py --data_dir $DATA_DIR \
                   --annotate_hits_table {input.seed} \
                   --no_file_comments \
                   -o {params.out} \
                   --output_dir {params.outdir} \
                   --override \
                   --cpu {threads} >> {log} 2>&1

        # cleanup counter
        (
            flock -x 200
            COUNT=$(cat "$COUNTER_FILE")
            COUNT=$((COUNT-1))
            if [ $COUNT -le 0 ]; then
                echo 0 > "$COUNTER_FILE"
                rm -rf "$DATA_DIR"
                rm -f "$DONE_FILE"
            else
                echo $COUNT > "$COUNTER_FILE"
            fi
        ) 200>"$COUNTER_FILE.lock"

        touch {output.done}

        """

rule mag_annotate__eggnog_merge_annotations:
    input:
        mag_annotate_aggregate_annotations_eggnog_search
    output:
        EGGNOG / "eggnog.emapper.annotations"
    log:
        EGGNOG / "merge_annotations.log"
    shell:
        """
        cat {input} > {output}
        """
