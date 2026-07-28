def mag_annotate_aggregate_annotations_eggnog_search(wildcards):
    checkpoint_output = checkpoints.mag_annotate__eggnog_split.get(**wildcards).output[0]
    files = glob_wildcards(os.path.join(checkpoint_output, "seed.chunk.{i}.seed_orthologs")).i
    return expand(
        os.path.join(checkpoint_output, "seed.chunk.{i}.emapper.annotations"),
        i=files
    )

def mag_annotate_aggregate_annotations_eggnog_mags(wildcards):
    checkpoint_output = checkpoints.mag_annotate__eggnog_mags__split.get(**wildcards).output[0]
    files = glob_wildcards(
        os.path.join(checkpoint_output, "seed.chunk.{i}.seed_orthologs")
    ).i
    return expand(
        os.path.join(checkpoint_output, "seed.chunk.{i}.emapper.annotations"),
        i=files,
    )