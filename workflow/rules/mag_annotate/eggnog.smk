rule mag_annotate__eggnog:
    """Run eggnog over the dereplicated mags"""
    input:
        contigs=DREP / "dereplicated_genomes.fa.gz",
    output:
        seed = EGGNOG / "eggnog.emapper.seed_orthologs",
        genepred_fasta = EGGNOG / "eggnog.emapper.genepred.fasta",
        genepred_gff   = EGGNOG / "eggnog.emapper.genepred.gff"
    log:
        protected(EGGNOG / "eggnog.log"),
    container:
        docker["eggnog"]
    params:
        out_dir=EGGNOG,
        db=features["databases"]["eggnog"],
        prefix="eggnog"
    threads: esc("cpus", "mag_annotate__eggnog")
    resources:
        runtime=esc("runtime", "mag_annotate__eggnog"),
        mem_mb=esc("mem_mb", "mag_annotate__eggnog"),
        cpus_per_task=esc("cpus", "mag_annotate__eggnog"),
        slurm_partition=esc("partition", "mag_annotate__eggnog"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'mag_annotate__eggnog')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("mag_annotate__eggnog"))
    shell:
        """
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
                   -o {params.prefix}  \
                   2>> {log} 1>&2;
        """

checkpoint mag_annotate__eggnog_split:
    """
    Split the MAG eggnog seed_orthologs into chunks (equivalent to contig cut_prodigal)
    """
    input:
        seed = EGGNOG / "eggnog.emapper.seed_orthologs"
    output:
        directory(EGGNOG / "Chunks")
    log:
        EGGNOG / "eggnog_split.log"
    params:
        chunk_size = params["mag_annotate"]["eggnog"]["split"],
        prefix = EGGNOG / "Chunks" / "seed.chunk."
    threads: 1
    shell:
        """
        mkdir -p {output}

        split -d -l {params.chunk_size} {input.seed} {params.prefix} --additional-suffix=.seed_orthologs

        """