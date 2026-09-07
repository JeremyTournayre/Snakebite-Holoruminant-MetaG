# workflow/rules/read_annotate/kofam_rumen.smk

rule read_annotate__kofam_rumen_count:
    """Aligne les reads sur le catalogue rumen et compte par gène"""
    input:
        forwards=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_1.fq.gz",
        reverses=PRE_BOWTIE2 / "decontaminated_reads" / "{sample_id}.{library_id}_2.fq.gz",
        catalogue=features["databases"]["kofam_rumen_fna"],
        index=features["databases"]["kofam_rumen_fna"] + ".bwt",
    output:
        counts=PRE_KOFAM_RUMEN / "{sample_id}.{library_id}.counts.tsv.gz",
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

        PREFIX=$(echo {output.counts} | sed 's/.tsv.gz//')
        SAM=$PREFIX.tmp.sam
        BAM_SORTED=$PREFIX.tmp.sorted.bam
        BAM_FILTERED=$PREFIX.tmp.filtered.bam

        bwa mem -t {threads} {input.catalogue} {input.forwards} {input.reverses} \
            > $SAM 2>> {log}.{resources.attempt}

        samtools sort -@ {threads} $SAM \
            > $BAM_SORTED 2>> {log}.{resources.attempt}
        rm $SAM

        msamtools filter -b -l 80 -p 95 -z 5 $BAM_SORTED \
            > $BAM_FILTERED 2>> {log}.{resources.attempt}
        rm $BAM_SORTED

        msamtools profile --multi prop --label {wildcards.sample_id}.{wildcards.library_id} \
            -o {output.counts} \
            $BAM_FILTERED \
            2>> {log}.{resources.attempt}
        rm $BAM_FILTERED

        mv {log}.{resources.attempt} {log}
        """


rule read_annotate__kofam_rumen_merge_ko:
    """Croise les comptages avec l'annotation KofamScan → table KO x échantillon"""
    input:
        counts=PRE_KOFAM_RUMEN / "{sample_id}.{library_id}.counts.tsv.gz",
        ko_annot=features["databases"]["kofam_rumen_annotation"],
    output:
        PRE_KOFAM_RUMEN / "{sample_id}.{library_id}.ko_counts.tsv",
    log:
        PRE_KOFAM_RUMEN / "{sample_id}.{library_id}.ko_counts.log",
    shell:
        """
        awk '
            NR==FNR && NR==1 {{next}}
            NR==FNR {{ko[$2]=$3; next}}
            /^#/ {{next}}
            ($1 in ko) {{print ko[$1], $0}}
        ' {input.ko_annot} <(zcat {input.counts} | grep -v "^#") \
            > {output} 2> {log}
        """

rule read_annotate__kofam_rumen_aggregate_ko:
    """Agrège les proportions par KO, ajoute le nom fonctionnel et une ligne Unknown"""
    input:
        ko_counts=PRE_KOFAM_RUMEN / "{sample_id}.{library_id}.ko_counts.tsv",
    output:
        PRE_KOFAM_RUMEN / "{sample_id}.{library_id}.ko_named.tsv",
    log:
        PRE_KOFAM_RUMEN / "{sample_id}.{library_id}.ko_named.log",
    params:
        ko_list=str(Path(features["databases"]["kofam_rumen_annotation"]).parent.parent / "kofam_db" / "ko_list"),
    shell:
        """
        awk '
            NR==FNR && NR==1 {{next}}
            NR==FNR {{
                def=""; for(i=12;i<=NF;i++) def=def (i>12?" ":"") $i
                name[$1]=def; next
            }}
            {{ko[$1]+=$3; total+=$3}}
            END {{
                for (k in ko) {{
                    if (k in name) print k"\t"name[k]"\t"ko[k]
                    else print k"\tUnknown function\t"ko[k]
                }}
                print "Unknown\tUnknown function\t"(1-total)
            }}
        ' {params.ko_list} {input.ko_counts} \
            | sort -k1,1 \
            > {output} 2> {log}
        """


rule read_annotate__kofam_rumen:
    """Run KO annotation rumen custom"""
    input:
        [
            PRE_KOFAM_RUMEN / f"{sample_id}.{library_id}.ko_named.tsv"
            for sample_id, library_id in SAMPLE_LIBRARY
        ]