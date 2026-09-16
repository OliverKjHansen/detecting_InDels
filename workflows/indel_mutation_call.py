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
assembly_type     = config["assembly_type"]
min_contig_length = config["min_contig_length"]
inpath_fa         = expanduser(config["paths"]["fasta"])
inpath_bam        = expanduser(config["paths"]["bam"])
output            = os.path.abspath(config["paths"]["output"])
results           = os.path.abspath(config["paths"]["results"])
plots             = os.path.abspath(config["paths"]["plots"])
scripts           = os.path.abspath("scripts")

# Repeat models whose tracts are masked out, from repetitive_dna.py.
models = ["M0", "M1", "M2"]

# Steps 2-4 split each contig into windows of this many bp, one job per window.
chunk_size = 40_000_000

# Final filtering threshold: minimum Phred base quality across the InDel 4-mer.
min_kmerq = 93


# Resource requests.
gwf = Workflow(
    defaults={
        'cores': 4,
        'memory': '10g',
        'walltime': '48:00:00',
        'account': 'mutationalscanning'})


for species, samples in sample_dict.items():
    for sample in samples:
        assembly    = f"{inpath_fa}/{species}/{sample}/{assembly_type}_assembly/{sample}_{assembly_type}.fa"
        haploid     = f"{inpath_fa}/{species}/{sample}/haploid_assembly/{sample}_haploid.fa"
        alignment   = f"{inpath_bam}/{species}/{sample}/kinetics/{sample}_kinetics_{assembly_type}.bam"
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
        calls_haploid_dir = f"{workdir}/calls_haploid"
        opportunities_dir = f"{workdir}/opportunities"
        fastq_dir         = f"{workdir}/fastq"
        haploid_dir       = f"{workdir}/haploid"
        merged_dir        = f"{workdir}/merged"

        for contig, length in long_contigs.items():
            repeat_beds   = [f"{repeats_dir}/{contig}_{model}.bed" for model in models]
            unresolved    = f"{precompute}/unresolved_haplotypes/{contig}.bed"
            mask          = f"{mask_dir}/{contig}.bed"
            calls         = f"{calls_dir}/{contig}.tsv"
            fastq         = f"{fastq_dir}/{contig}.fastq"
            opportunities = f"{opportunities_dir}/{contig}.tsv"

            # Step 1 — the mask: every base any repeat model calls repetitive,
            # plus the unresolved-haplotype regions. The repeat BEDs carry a
            # header line, which is dropped before merging.
            repeat_bed_list = " ".join(repeat_beds)

            gwf.target(
                f"Mask_{contig}_{sample}_{assembly_type}",
                inputs=repeat_beds + [unresolved],
                outputs=[mask],
                cores=1,
                memory="8g",
                walltime="01:00:00",
            ) << f"""
            set -e
            mkdir -p {mask_dir}
            WORK=$(mktemp -d)
            trap 'rm -rf "$WORK"' EXIT
            awk 'FNR>1 {{print $1 "\\t" $2 "\\t" $3}}' {repeat_bed_list} > $WORK/mask.bed
            cut -f1-3 {unresolved} >> $WORK/mask.bed
            sort -k1,1 -k2,2n $WORK/mask.bed | bedtools merge -i - > {mask}.tmp
            mv {mask}.tmp {mask}
            """

            # Steps 2-4 run per window. The windows [0, chunk_size),
            # [chunk_size, 2*chunk_size), ... tile the contig, and every pileup
            # position belongs to exactly one window, so each call and each
            # opportunity is counted once.
            window_calls, window_fastqs, window_opportunities = [], [], []

            for start in range(0, length, chunk_size):
                end = min(start + chunk_size, length)
                w_calls         = f"{calls_dir}/chunks/{contig}_{start}_{end}.tsv"
                w_fastq         = f"{fastq_dir}/chunks/{contig}_{start}_{end}.fastq"
                w_opportunities = f"{opportunities_dir}/chunks/{contig}_{start}_{end}.tsv"
                window_calls.append(w_calls)
                window_fastqs.append(w_fastq)
                window_opportunities.append(w_opportunities)

                # Step 2a — call 1-3 bp InDels outside the mask, dropping those
                # in a homopolymer context.
                gwf.target(
                    f"SmallIndels_{contig}_{start}_{sample}_{assembly_type}",
                    inputs=[assembly, alignment, coverage, mask] + tissue_inputs,
                    outputs=[w_calls],
                ) << f"""
                mkdir -p {calls_dir}/chunks
                python {scripts}/small_indels.py {mask} {contig} {alignment} {assembly} {coverage} {tissue} {start} {end} > {w_calls}
                """

                # Step 3a — extract the reads carrying a call, for haploid remapping.
                gwf.target(
                    f"BuildFastq_{contig}_{start}_{sample}_{assembly_type}",
                    inputs=[w_calls, alignment],
                    outputs=[w_fastq],
                ) << f"""
                mkdir -p {fastq_dir}/chunks
                python {scripts}/build_mutation_fastq.py {w_calls} {alignment} {contig} {w_fastq} {start} {end}
                """

                # Step 4a — count opportunities outside the mask (the
                # denominator of the rate), with the same homopolymer rule.
                gwf.target(
                    f"Opportunities_{contig}_{start}_{sample}_{assembly_type}",
                    inputs=[assembly, alignment, coverage, mask] + tissue_inputs,
                    outputs=[w_opportunities],
                ) << f"""
                mkdir -p {opportunities_dir}/chunks
                python {scripts}/opportunities.py {mask} {contig} {alignment} {assembly} {coverage} {tissue} {start} {end} > {w_opportunities}
                """

            window_calls_list         = " ".join(window_calls)
            window_fastqs_list        = " ".join(window_fastqs)
            window_opportunities_list = " ".join(window_opportunities)

            # Step 2b-4b — join the windows back into one file per contig.
            #   calls:         concatenated, header kept from the first window
            #   FASTQ:         concatenated, a read carrying calls in two
            #                  windows is kept once
            #   opportunities: counts per (contig, context, tissue, mapQ, kmerQ),
            #                  so they are summed over the windows
            gwf.target(
                f"MergeWindows_{contig}_{sample}_{assembly_type}",
                inputs=window_calls + window_fastqs + window_opportunities,
                outputs=[calls, fastq, opportunities],
                cores=1,
                memory="16g",
                walltime="02:00:00",
            ) << f"""
            mkdir -p {calls_dir} {fastq_dir} {opportunities_dir}
            awk 'FNR==1 && NR==1 || FNR>1' {window_calls_list} > {calls}
            cat {window_fastqs_list} | awk 'NR % 4 == 1 {{keep = !seen[$1]++}} keep' > {fastq}
            awk -F'\\t' -v OFS='\\t' 'FNR==1 {{if (NR==1) print; next}} {{n[$1 OFS $2 OFS $3 OFS $4 OFS $5] += $6}} END {{for (key in n) print key, n[key]}}' \\
                {window_opportunities_list} > {opportunities}
            """

        fastqs            = [f"{fastq_dir}/{contig}.fastq" for contig in long_contigs]
        fastq_list        = " ".join(fastqs)
        combined_fastq    = f"{haploid_dir}/all.fastq"
        remapped_bam      = f"{haploid_dir}/remapped.bam"
        haploid_locations = f"{haploid_dir}/locations.tsv"

        # Step 5 — remap all call-carrying reads to the haploid assembly. The
        # FASTQs are listed explicitly and the combined FASTQ is written to
        # haploid/, so a re-run never reads its own output.
        gwf.target(
            f"HaploidMapping_{sample}_{assembly_type}",
            inputs=fastqs + [haploid],
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

        # Step 6 — haploid concordance, per contig: is the same InDel also in
        # the read's alignment to the haploid assembly?
        for contig in long_contigs:
            calls         = f"{calls_dir}/{contig}.tsv"
            calls_haploid = f"{calls_haploid_dir}/{contig}.tsv"

            gwf.target(
                f"HaploidFilter_{contig}_{sample}_{assembly_type}",
                inputs=[calls, haploid_locations],
                outputs=[calls_haploid],
            ) << f"""
            mkdir -p {calls_haploid_dir}
            python {scripts}/filtering_haploid_mappping.py {calls} {haploid_locations} > {calls_haploid}.tmp
            awk 'NR==1 {{print $0 "\\t" "species" "\\t" "sample"}}; NR>1 {{print $0 "\\t" "{species}" "\\t" "{sample}"}}' {calls_haploid}.tmp > {calls_haploid}
            rm -f {calls_haploid}.tmp
            """

        calls_haploid_files  = [f"{calls_haploid_dir}/{contig}.tsv" for contig in long_contigs]
        opportunities_files  = [f"{opportunities_dir}/{contig}.tsv" for contig in long_contigs]
        calls_haploid_list   = " ".join(calls_haploid_files)
        opportunities_list   = " ".join(opportunities_files)
        merged_indels        = f"{merged_dir}/indels.tsv"
        merged_opportunities = f"{merged_dir}/opportunities.tsv"
        plots_pdf            = f"{plotsdir}/plots.pdf"

        # Step 7 — merge the contigs (header kept from the first file only) and
        # draw diagnostics of the unfiltered calls. species/sample are added to
        # the opportunities; the calls got them in step 6.
        gwf.target(
            f"MergeAndPlot_{sample}_{assembly_type}",
            inputs=calls_haploid_files + opportunities_files,
            outputs=[merged_indels, merged_opportunities, plots_pdf],
            cores=8,
            memory="100g",
            walltime="24:00:00",
        ) << f"""
        mkdir -p {merged_dir} {plotsdir}
        awk 'FNR==1 && NR==1 || FNR>1' {calls_haploid_list} > {merged_indels}
        awk 'FNR==1 && NR==1 || FNR>1' {opportunities_list} \\
            | awk 'NR==1 {{print $0 "\\t" "species" "\\t" "sample"}}; NR>1 {{print $0 "\\t" "{species}" "\\t" "{sample}"}}' \\
            > {merged_opportunities}
        Rscript {scripts}/plot_indels.R {merged_indels} {merged_opportunities} {plots_pdf}
        """

        # Step 8 — final filtering, with a report of rows removed per step.
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
