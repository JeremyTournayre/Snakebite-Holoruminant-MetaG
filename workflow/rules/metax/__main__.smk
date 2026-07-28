include: "metax.smk"

rule metax:
    input:
        expand(METAX /"{sample}"/ "{sample}_Fasta_MAGs.faa", sample=SAMPLES),
        expand(METAX /"{sample}"/ "{sample}_Function_MAGs.tsv", sample=SAMPLES),
        expand(METAX /"{sample}"/ "{sample}_Taxonomy_MAGs.tsv", sample=SAMPLES),
        expand(METAX /"{sample}"/ "{sample}_Fasta_Contigs.faa", sample=SAMPLES),
        expand(METAX /"{sample}"/ "{sample}_Function_Contigs.tsv", sample=SAMPLES),
        expand(METAX /"{sample}"/ "{sample}_Taxonomy_Contigs.tsv", sample=SAMPLES)