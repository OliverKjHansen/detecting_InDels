import os
from os.path import expanduser

import yaml
from gwf import Workflow

# ---------------------------------------------------------------------------
# InDels OUTSIDE repetitive DNA. Calls 1-3 bp insertions and deletions in single
# reads at positions that are not masked. The mask is the union of every repeat
# model from repetitive_dna.py (M0, M1, M2) and the unresolved-haplotype regions
# from precompute.py.
#
# Homopolymers are not in the mask, because repetitive_dna.py does not annotate
# k=1. They are handled inside the caller instead: an InDel whose sequence plus
# one flanking base on each side changes base at most once is dropped, and the
# opportunity counter applies the same rule, so numerator and denominator agree.
#
# The BAM is read once per window: one scan job writes the calls, the
# opportunities and a FASTQ of the reads carrying a call. The steps around it
# run once per sample, each looping over the contigs.
#
# Configuration is read from config.yaml. Run gwf from the repository root,
# after precompute.py and repetitive_dna.py:
#
#   gwf -f workflows/precompute.py run
#   gwf -f workflows/repetitive_dna.py run
#   gwf -f workflows/indel_mutation_call.py run
#
# output/, results/, plots/ and scripts/ are made absolute from there, because
# gwf runs the recipes from workflows/.
# ---------------------------------------------------------------------------

with open("config.yaml") as handle:
    config = yaml.safe_load(handle)

sample_dict       = config["samples"]
tissue_dict       = config["tissue_samples"]
methylation_dict  = config.get("methylation_samples", {})
inpath_classification = expanduser(config["paths"]["classification"])
assembly_type     = config["assembly_type"]
min_contig_length = config["min_contig_length"]
inpath_fa         = expanduser(config["paths"]["fasta"])
inpath_bam        = expanduser(config["paths"]["bam"])
alignment_overrides = config.get("alignment_overrides", {})
output            = os.path.abspath(config["paths"]["output"])
results           = os.path.abspath(config["paths"]["results"])
plots             = os.path.abspath(config["paths"]["plots"])
scripts           = os.path.abspath("scripts")

# Repeat models whose tracts are masked out, from repetitive_dna.py.
models = ["M0", "M1", "M2"]

# Step 2 splits each contig into windows of this many bp, one scan job per window.
chunk_size = 10_000_000

# Final filtering threshold: minimum Phred base quality across the InDel 4-mer.
min_kmerq = 93


# Resource requests.
gwf = Workflow(
    defaults={
        'cores': 4,
        'memory': '10g',
        'walltime': '12:00:00',
        'account': 'mutationalscanning'})


for species, samples in sample_dict.items():
    for sample in samples:
        assembly    = f"{inpath_fa}/{species}/{sample}/{assembly_type}_assembly/{sample}_{assembly_type}.fa"
        haploid     = f"{inpath_fa}/{species}/{sample}/haploid_assembly/{sample}_haploid.fa"
        # Samples listed under alignment_overrides in config.yaml use that BAM.
        alignment   = expanduser(alignment_overrides.get(species, {}).get(
            sample, f"{inpath_bam}/{species}/{sample}/kinetics/{sample}_kinetics_{assembly_type}.bam"))
        precompute  = f"{output}/precomute/{species}/{sample}"
        coverage    = f"{precompute}/avg_coverage.txt"
        repeats_dir = f"{output}/repetitive_dna/{species}/{sample}"
        workdir     = f"{output}/nonrepetitive_indels/{species}/{sample}"
        resultsdir  = f"{results}/nonrepetitive_indels/{species}/{sample}"
        plotsdir    = f"{plots}/nonrepetitive_indels/{species}/{sample}"
        contig_dict = {}

        with open(f"{assembly}.fai") as contigs:
            for line in contigs:
                contig_dict[line.split()[0]] = line.split()[1]

        long_contigs = {contig: int(length) for contig, length in contig_dict.items() if int(length) > min_contig_length}

        # Samples without read labels have no tissue table; the callers then
        # report tissue NA for every read.
        if sample in tissue_dict.get(species, []):
            tissue        = f"{precompute}/sorted_tissue_classification.tsv"
            tissue_inputs = [tissue]
        else:
            tissue        = "NA"
            tissue_inputs = []

        mask_dir          = f"{workdir}/mask"
        calls_dir         = f"{workdir}/calls"
        opportunities_dir = f"{workdir}/opportunities"
        fastq_dir         = f"{workdir}/fastq"
        haploid_dir       = f"{workdir}/haploid"
        merged_dir        = f"{workdir}/merged"

        masks         = {contig: f"{mask_dir}/{contig}.bed" for contig in long_contigs}
        calls         = {contig: f"{calls_dir}/{contig}.tsv" for contig in long_contigs}
        opportunities = {contig: f"{opportunities_dir}/{contig}.tsv" for contig in long_contigs}
        fastqs        = {contig: f"{fastq_dir}/{contig}.fastq" for contig in long_contigs}

        # Step 1 — the mask of every contig: every base any repeat model calls
        # repetitive, plus the unresolved-haplotype regions. The repeat BEDs
        # carry a header line, which is dropped before merging.
        repeat_beds     = [f"{repeats_dir}/{contig}_{model}.bed" for contig in long_contigs for model in models]
        unresolved_beds = [f"{precompute}/unresolved_haplotypes/{contig}.bed" for contig in long_contigs]
        contig_repeat_beds = " ".join(f"{repeats_dir}/${{contig}}_{model}.bed" for model in models)
        contig_names       = " ".join(long_contigs)

        gwf.target(
            f"Mask_{sample}_{assembly_type}",
            inputs=repeat_beds + unresolved_beds,
            outputs=list(masks.values()),
            cores=1,
            memory="8g",
            walltime="04:00:00",
        ) << f"""
        set -e
        mkdir -p {mask_dir}
        WORK=$(mktemp -d)
        trap 'rm -rf "$WORK"' EXIT
        for contig in {contig_names}; do
            awk 'FNR>1 {{print $1 "\\t" $2 "\\t" $3}}' {contig_repeat_beds} > $WORK/mask.bed
            cut -f1-3 {precompute}/unresolved_haplotypes/${{contig}}.bed >> $WORK/mask.bed
            sort -k1,1 -k2,2n $WORK/mask.bed | bedtools merge -i - > {mask_dir}/${{contig}}.bed.tmp
            mv {mask_dir}/${{contig}}.bed.tmp {mask_dir}/${{contig}}.bed
        done
        """

        # Step 2 runs per window. The windows [0, chunk_size),
        # [chunk_size, 2*chunk_size), ... tile each contig, and every pileup
        # position belongs to exactly one window, so each call and each
        # opportunity is counted once.
        window_calls         = {contig: [] for contig in long_contigs}
        window_opportunities = {contig: [] for contig in long_contigs}
        window_fastqs        = {contig: [] for contig in long_contigs}

        for contig, length in long_contigs.items():
            for start in range(0, length, chunk_size):
                end = min(start + chunk_size, length)
                w_calls         = f"{calls_dir}/chunks/{contig}_{start}_{end}.tsv"
                w_opportunities = f"{opportunities_dir}/chunks/{contig}_{start}_{end}.tsv"
                w_fastq         = f"{fastq_dir}/chunks/{contig}_{start}_{end}.fastq"
                window_calls[contig].append(w_calls)
                window_opportunities[contig].append(w_opportunities)
                window_fastqs[contig].append(w_fastq)

                # Step 2 — one pass over the window: call 1-3 bp InDels outside
                # the mask, dropping those in a homopolymer context; count the
                # opportunities (the denominator of the rate) with the same
                # homopolymer rule; and write the reads carrying a call, for
                # haploid remapping.
                gwf.target(
                    f"Scan_{contig}_{start}_{sample}_{assembly_type}",
                    inputs=[assembly, alignment, coverage, masks[contig]] + tissue_inputs,
                    outputs=[w_calls, w_opportunities, w_fastq],
                    cores=1,
                ) << f"""
                mkdir -p {calls_dir}/chunks {opportunities_dir}/chunks {fastq_dir}/chunks
                python {scripts}/nonrepetitive_scan.py {masks[contig]} {contig} {alignment} {assembly} {coverage} {tissue} {start} {end} {w_calls} {w_opportunities} {w_fastq}
                """

        # Step 3 — join the windows back into one file per contig, for every
        # contig in one job.
        #   calls:         concatenated, header kept from the first window
        #   opportunities: counts per (contig, context, tissue, mapQ, kmerQ),
        #                  so they are summed over the windows
        #   FASTQ:         concatenated, a read carrying calls in two windows is
        #                  kept once
        merge_commands = []
        for contig in long_contigs:
            merge_commands.append(
                f"awk 'FNR==1 && NR==1 || FNR>1' {' '.join(window_calls[contig])} > {calls[contig]}")
            merge_commands.append(
                f"awk -F'\\t' -v OFS='\\t' 'FNR==1 {{if (NR==1) print; next}} "
                f"{{n[$1 OFS $2 OFS $3 OFS $4 OFS $5] += $6}} END {{for (key in n) print key, n[key]}}' "
                f"{' '.join(window_opportunities[contig])} > {opportunities[contig]}")
            merge_commands.append(
                f"cat {' '.join(window_fastqs[contig])} | awk 'NR % 4 == 1 {{keep = !seen[$1]++}} keep' > {fastqs[contig]}")
        merge_commands = "\n        ".join(merge_commands)

        gwf.target(
            f"MergeWindows_{sample}_{assembly_type}",
            inputs=[path for contig in long_contigs
                    for path in window_calls[contig] + window_opportunities[contig] + window_fastqs[contig]],
            outputs=list(calls.values()) + list(opportunities.values()) + list(fastqs.values()),
            cores=1,
            memory="16g",
            walltime="12:00:00",
        ) << f"""
        set -e
        mkdir -p {calls_dir} {opportunities_dir} {fastq_dir}
        {merge_commands}
        """

        fastq_list        = " ".join(fastqs.values())
        combined_fastq    = f"{haploid_dir}/all.fastq"
        remapped_bam      = f"{haploid_dir}/remapped.bam"
        haploid_locations = f"{haploid_dir}/locations.tsv"

        # Step 4 — remap all call-carrying reads to the haploid assembly. The
        # FASTQs are listed explicitly and the combined FASTQ is written to
        # haploid/, so a re-run never reads its own output.
        gwf.target(
            f"HaploidMapping_{sample}_{assembly_type}",
            inputs=list(fastqs.values()) + [haploid],
            outputs=[combined_fastq, remapped_bam, haploid_locations],
            cores=16,
            memory="100g",
            walltime="24:00:00",
        ) << f"""
        mkdir -p {haploid_dir}
        cat {fastq_list} > {combined_fastq}
        pbmm2 align {haploid} {combined_fastq} {remapped_bam} --preset HIFI --sort -j 12 -J 4 -m 8G
        samtools view {remapped_bam} | awk -v OFS='\\t' '{{print $1,$2,$3,$4,$5,$6}}' | sort -k1,1 > {haploid_locations}
        """

        calls_list           = " ".join(calls.values())
        opportunities_list   = " ".join(opportunities.values())
        merged_indels        = f"{merged_dir}/indels.tsv"
        merged_opportunities = f"{merged_dir}/opportunities.tsv"
        plots_pdf            = f"{plotsdir}/plots.pdf"

        # Step 5 — haploid concordance for all contigs: is the same InDel also
        # in the read's alignment to the haploid assembly? locations.tsv is read
        # once for every contig. The contigs are merged into one table (header
        # kept once) with species/sample added.
        gwf.target(
            f"HaploidFilter_{sample}_{assembly_type}",
            inputs=list(calls.values()) + [haploid_locations],
            outputs=[merged_indels],
            cores=1,
            memory="20g",
            walltime="12:00:00",
        ) << f"""
        mkdir -p {merged_dir}
        python {scripts}/filtering_haploid_mappping.py {calls_list} {haploid_locations} > {merged_indels}.tmp
        awk 'NR==1 {{print $0 "\\t" "species" "\\t" "sample"}}; NR>1 {{print $0 "\\t" "{species}" "\\t" "{sample}"}}' {merged_indels}.tmp > {merged_indels}
        rm -f {merged_indels}.tmp
        """

        # Step 6 — merge the contigs' opportunities (header kept from the first
        # file only, species/sample added) and draw diagnostics of the
        # unfiltered calls.
        gwf.target(
            f"MergeAndPlot_{sample}_{assembly_type}",
            inputs=[merged_indels] + list(opportunities.values()),
            outputs=[merged_opportunities, plots_pdf],
            cores=8,
            memory="100g",
            walltime="24:00:00",
        ) << f"""
        mkdir -p {merged_dir} {plotsdir}
        awk 'FNR==1 && NR==1 || FNR>1' {opportunities_list} \\
            | awk 'NR==1 {{print $0 "\\t" "species" "\\t" "sample"}}; NR>1 {{print $0 "\\t" "{species}" "\\t" "{sample}"}}' \\
            > {merged_opportunities}
        Rscript {scripts}/plot_indels.R {merged_indels} {merged_opportunities} {plots_pdf}
        """

        # Step 7 — final filtering, with a report of rows removed per step.
        # These small tables are the results.
        gwf.target(
            f"FinalFiltering_{sample}_{assembly_type}",
            inputs=[merged_indels, merged_opportunities],
            outputs=[f"{resultsdir}/indels.tsv",
                     f"{resultsdir}/opportunities.tsv",
                     f"{resultsdir}/filter_report.tsv"],
            cores=4,
            memory="50g",
            walltime="01:00:00",
        ) << f"""
        mkdir -p {resultsdir}
        python {scripts}/filtering_small_indels.py {merged_indels} \\
            -o {resultsdir}/indels.tsv \\
            --haploid-check True \\
            --min-quality {min_kmerq} \\
            --report {resultsdir}/filter_report.tsv
        python {scripts}/filtering_opportunities.py {merged_opportunities} \\
            -o {resultsdir}/opportunities.tsv \\
            --min-kmerq {min_kmerq}
        """

        # Step 8 — cell-type classification, for the samples listed under
        # methylation_samples in config.yaml. Every read of those samples has a
        # log-likelihood per cell type; the script contrasts the germline
        # (Sperm.*) columns with a set of somatic ones and appends the
        # difference and its bin to every call. The likelihood table covers the
        # whole sample, so it is streamed once and filtered to the reads that
        # carry a call.
        if sample in methylation_dict.get(species, []):
            likelihoods = (f"{inpath_classification}/{species}/{sample}/{sample}"
                           f"/05_ReadsClassification/{sample}_All_Reads_Likelihood_merged.tsv")

            gwf.target(
                f"Classify_{sample}_{assembly_type}",
                inputs=[f"{resultsdir}/indels.tsv", likelihoods],
                outputs=[f"{resultsdir}/indels_classified.tsv"],
                cores=1,
                memory="8g",
                walltime="04:00:00",
            ) << f"""
            set -e
            bash {scripts}/join_read_classification.sh \\
                {resultsdir}/indels.tsv \\
                {likelihoods} \\
                {resultsdir}/indels_classified.tsv
            """
