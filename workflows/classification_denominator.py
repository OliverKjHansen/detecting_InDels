from functools import partial

from gwf import Workflow

from pipeline import Pipeline

# ---------------------------------------------------------------------------
# A denominator per germline/somatic bin, for the samples listed under
# methylation_samples in config.yaml.
#
# Both indel pipelines append a classification to every call: the read's best
# germline likelihood minus its best somatic one, and the bin that difference
# falls in. The opportunity tables carry no such column, because they count read
# positions per contig, context, tissue and mapQ with read identity aggregated
# away. So the bins have numerators and no denominator.
#
# This workflow counts, over every read of the sample, how many reads fall in
# each bin and how many CpGs those reads carry. A bin can then be given a share
# of the sample's read positions P, which the opportunity tables do hold:
#
#     P_b ~= P * reads_b / reads_total      flat
#     P_b ~= P * CpGs_b  / CpGs_total       length-weighted
#
# The flat split assumes positions per read do not depend on the bin, which is
# not quite true: longer reads carry more CpGs, are classified more confidently
# and land in more extreme bins. The CpG column corrects most of that, so both
# are written and the distance between them is the size of the effect.
#
# The counts are a property of the sample, not of a pipeline or a repeat model,
# so one pass per sample writes the same file into both result trees, next to
# the calls each notebook already reads:
#
#     results/nonrepetitive_indels/{species}/{sample}/classification_reads.tsv
#     results/repetitive_indels/{species}/{sample}/classification_reads.tsv
#
# Nothing else is read or written, so this can run on its own, at any time,
# alongside or after the two indel pipelines:
#
#   gwf -f workflows/classification_denominator.py run
#
# Paths come from the `layout` templates in config.yaml, through pipeline.py.
# ---------------------------------------------------------------------------

pipeline = Pipeline()
scripts  = pipeline.scripts


# Resource requests. One streaming pass over a per-read table, so one core; the
# script keeps only the totals per bin, so memory is nominal.
gwf = Workflow(
    defaults={
        'cores': 1,
        'memory': '4g',
        'walltime': '06:00:00',
        'account': 'mutationalscanning'})


for species, samples in (pipeline.config.get("methylation_samples") or {}).items():
    for sample in samples:
        path        = partial(pipeline.path, species=species, sample=sample)
        likelihoods = path("likelihoods")
        nonrep      = path("nonrepetitive_classification_reads")
        rep         = path("repetitive_classification_reads")
        nonrep_dir  = pipeline.folder("nonrepetitive_classification_reads", species, sample)
        rep_dir     = pipeline.folder("repetitive_classification_reads", species, sample)

        gwf.target(
            f"ClassificationReads_{species}_{sample}",
            inputs=[likelihoods],
            outputs=[nonrep,
                     rep],
        ) << f"""
        set -e
        mkdir -p {nonrep_dir} {rep_dir}
        bash {scripts}/classification_read_counts.sh \\
            {likelihoods} \\
            {nonrep}
        cp {nonrep} {rep}
        """
