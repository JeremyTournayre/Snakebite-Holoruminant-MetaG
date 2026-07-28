rule mag_annotate__gtdbtk_allbins__classify:
    """Run GTDB-Tk over ALL raw bins (not dereplicated), one global run."""
    input:
        fasta_folder=DREP / "separated_bins",
        database=features["databases"]["gtdbtk"],
    output:
        summary=GTDBTK_ALLBINS / "gtdbtk.summary.tsv",
        bac120=GTDBTK_ALLBINS / "gtdbtk.bac120.summary.tsv",
        ar53=GTDBTK_ALLBINS / "gtdbtk.ar53.summary.tsv",
        align=GTDBTK_ALLBINS / "align.tar.gz",
        classify=GTDBTK_ALLBINS / "classify.tar.gz",
        identify=GTDBTK_ALLBINS / "identify.tar.gz",
    log:
        GTDBTK_ALLBINS / "gtdbtk_classify.log",
    container:
        docker["gtdbtk"]
    params:
        out_dir=GTDBTK_ALLBINS,
        ar53=GTDBTK_ALLBINS / "gtdbtk.ar53.summary.tsv",
        bac120=GTDBTK_ALLBINS / "gtdbtk.bac120.summary.tsv",
    threads: esc("cpus", "mag_annotate__gtdbtk_allbins__classify")
    resources:
        runtime=esc("runtime", "mag_annotate__gtdbtk_allbins__classify"),
        mem_mb=esc("mem_mb", "mag_annotate__gtdbtk_allbins__classify"),
        cpus_per_task=esc("cpus", "mag_annotate__gtdbtk_allbins__classify"),
        partition=esc("partition", "mag_annotate__gtdbtk_allbins__classify"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'mag_annotate__gtdbtk_allbins__classify')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("mag_annotate__gtdbtk_allbins__classify"))
    shell:
        """
        rm \
            --recursive \
            --force \
            {params.out_dir}/align \
            {params.out_dir}/classify \
            {params.out_dir}/identify \
            {params.out_dir}/gtdbtk_classify.log \
            {params.out_dir}/gtdbtk.json \
            {params.out_dir}/gtdbtk.summary.tsv \
            {params.out_dir}/gtdbtk.warnings.log \
            {params.ar53} \
            {params.bac120} \
        2>> {log}.{resources.attempt} 1>&2

        export GTDBTK_DATA_PATH="{input.database}"

        gtdbtk classify_wf \
            --genome_dir {input.fasta_folder} \
            --extension fa \
            --out_dir {params.out_dir} \
            --cpus {threads} \
            --skip_ani_screen \
        2>> {log}.{resources.attempt} 1>&2

        if [[ -s {params.ar53} ]] ; then
            ( csvstack --tabs {params.bac120} {params.ar53} \
            | csvformat --out-tabs > {output.summary} ) 2>> {log}.{resources.attempt}
        else
            cp {params.bac120} {output.summary} 2>> {log}.{resources.attempt}
            touch {output.ar53}
        fi

        for folder in align classify identify ; do
            tar --create --directory {params.out_dir} \
                --file {params.out_dir}/${{folder}}.tar.gz \
                --use-compress-program="pigz --processes {threads}" \
                --verbose ${{folder}} \
            2>> {log}.{resources.attempt} 1>&2
        done

        mv {log}.{resources.attempt} {log}
        """


rule mag_annotate__gtdbtk_allbins:
    input:
        rules.mag_annotate__gtdbtk_allbins__classify.output,