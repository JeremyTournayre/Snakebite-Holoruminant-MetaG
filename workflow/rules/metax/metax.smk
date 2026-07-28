rule metax__ncbi_taxdump:
    output:
        db_dir = directory(METAX / "ncbi_taxdump"),
        done = METAX / "ncbi_taxdump/.done"
    params:
        folder=config["pipeline_folder"]        
    container:
        docker["metax"]
    shell:
        """
        python {params.folder}/workflow/scripts/update_ncbi.py {output.db_dir}
        touch {output.done}
        """

rule metax__bakta_fix:
    input:
        bakta="results/mag_annotate/bakta_mags/bakta_{sample}.faa"
    output:
        fixed="results/mag_annotate/bakta_mags/bakta_{sample}_fixed_headers.faa"
    log:
        METAX / "{sample}.bakta_fix.log"
    params:
        folder=config["pipeline_folder"]
    container:
        docker["metax"]
    shell:
        """
        python {params.folder}/workflow/scripts/convert_bakta_header.py --sample {wildcards.sample} &> {log}
        """

rule metax__export:
    """Produce MetaX formatted outputs"""
    input:
        contig_faa="results/contig_annotate/prodigal/{sample}/{sample}.prodigal.fa",
        contig_eggnog="results/contig_annotate/eggnog/{sample}/eggnog_output.emapper.annotations",
        dram="results/mag_annotate/dram_mags/{sample}/annotate/genes.faa",
        bakta=rules.metax__bakta_fix.output.fixed,
        mag_eggnog=EGGNOGMAG / "{sample}" / "eggnog.emapper.annotations",
        mag_eggnog_faa=EGGNOGMAG / "{sample}" / "eggnog.emapper.genepred.fasta",
        gtdbtk_bac=GTDBTK_ALLBINS / "gtdbtk.bac120.summary.tsv",
        gtdbtk_ar=GTDBTK_ALLBINS / "gtdbtk.ar53.summary.tsv",
        taxdb=rules.metax__ncbi_taxdump.output.db_dir,
        taxdbdone=rules.metax__ncbi_taxdump.output.done
    output:
        contigs_faa = METAX /"{sample}"/ "{sample}_Fasta_Contigs.faa",
        contigs_func = METAX /"{sample}"/ "{sample}_Function_Contigs.tsv",
        contigs_tax = METAX /"{sample}"/ "{sample}_Taxonomy_Contigs.tsv",
        mags_faa = METAX /"{sample}"/ "{sample}_Fasta_MAGs.faa",
        mags_func = METAX /"{sample}"/ "{sample}_Function_MAGs.tsv",
        mags_tax = METAX /"{sample}"/ "{sample}_Taxonomy_MAGs.tsv",
        flag = METAX /"{sample}"/ "{sample}.metax.done"
    params:
        folder=config["pipeline_folder"],        
    log:
        METAX /"{sample}"/ "{sample}.metax.log"
    container:
        docker["metax"]
    threads: esc("cpus", "metax__export")
    resources:
        runtime=esc("runtime", "metax__export"),
        mem_mb=esc("mem_mb", "metax__export"),
        cpus_per_task=esc("cpus", "metax__export"),
        partition=esc("partition", "metax__export"),
        attempt=get_attempt,
    retries: len(get_escalation_order("metax__export"))
    shell:
        """
        python {params.folder}/workflow/scripts/For_metaX.py --taxdb {input.taxdb} --sample {wildcards.sample} --outdir {METAX}/{wildcards.sample} &>> {log}
        touch {output.flag}
        """