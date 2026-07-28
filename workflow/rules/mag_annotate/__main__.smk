include: "__functions__.smk"
include: "bakta.smk"
include: "camper.smk"
include: "quast.smk"
include: "gtdbtk.smk"
include: "dram.smk"
include: "dram_mag.smk"
include: "eggnog.smk"
include: "eggnog_annotations.smk"
include: "checkm2.smk"
include: "proteinortho.smk"
include: "phylophlan.smk"
include: "gtdbtk_allbins.smk"
include: "eggnog_mags.smk"

rule mag_annotate:
    """Evaluate the dereplication steps"""
    input:
        rules.mag_annotate__bakta.output,
        rules.mag_annotate__bakta_mags.input,
        rules.mag_annotate__camper.input,
        rules.mag_annotate__eggnog.output,
        rules.mag_annotate__eggnog_merge_annotations.output,
        rules.mag_annotate__quast.output,
        rules.mag_annotate__checkm2__predict.output,
        rules.mag_annotate__gtdbtk__classify.output,
        rules.mag_annotate__dram.input,
        rules.mag_annotate__dram_mags.input,
        rules.mag_annotate__proteinortho.output,
        rules.mag_annotate__phylophlan.output,
        rules.mag_annotate__gtdbtk_allbins.input,
        rules.mag_annotate__eggnog_mags.input,