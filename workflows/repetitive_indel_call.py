import os
from os.path import expanduser

import yaml
from gwf import Workflow

# ---------------------------------------------------------------------------
# InDels inside repetitive DNA. For every repeat model written by
# repetitive_dna.py (M0, M1, M2), scans the positions inside that model's
# repeat tracts and calls InDels whose size is a multiple of the period k.
#
# The BAM is read once per window for all three models: one scan job writes the
# calls and opportunities of every model plus one FASTQ of the reads carrying a
# call in any of them. Those reads are remapped to the haploid assembly once
# per sample, and each model's calls are checked against that remapping.
#
# Configuration is read from config.yaml. Run gwf from the repository root,
# after precompute.py and repetitive_dna.py:
#
#   gwf -f workflows/precompute.py run
#   gwf -f workflows/repetitive_dna.py run
#   gwf -f workflows/repetitive_indel_call.py run
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

# Repeat models written by repetitive_dna.py.
models = ["M0", "M1", "M2"]

# Step 1 splits each contig into windows of this many bp, one scan job per window.
chunk_size = 10_000_000

# Final filtering thresholds.
min_kmerq = 93                 # minimum Phred base quality across the InDel 4-mer
min_read_end_distance = 100    # bp; InDels nearer a read end are unreliable


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
        sample_dir  = f"{output}/repetitive_indels/{species}/{sample}"
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

        # Shared by the models: the reads carrying calls, and their remapping.
        fastq_dir   = f"{sample_dir}/fastq"
        haploid_dir = f"{sample_dir}/haploid"

        # Per model.
        calls_dir         = {model: f"{sample_dir}/{model}/calls" for model in models}
        opportunities_dir = {model: f"{sample_dir}/{model}/opportunities" for model in models}
        merged_dir        = {model: f"{sample_dir}/{model}/merged" for model in models}

        for contig, length in long_contigs.items():
            repeat_beds   = {model: f"{repeats_dir}/{contig}_{model}.bed" for model in models}
            calls         = {model: f"{calls_dir[model]}/{contig}.tsv" for model in models}
            opportunities = {model: f"{opportunities_dir[model]}/{contig}.tsv" for model in models}
            fastq         = f"{fastq_dir}/{contig}.fastq"

            # Step 1 runs per window. The windows [0, chunk_size),
            # [chunk_size, 2*chunk_size), ... tile the contig, and every pileup
            # position belongs to exactly one window, so each call and each
            # opportunity is counted once.
            window_calls         = {model: [] for model in models}
            window_opportunities = {model: [] for model in models}
            window_fastqs        = []

            for start in range(0, length, chunk_size):
                end = min(start + chunk_size, length)
                w_calls         = {model: f"{calls_dir[model]}/chunks/{contig}_{start}_{end}.tsv" for model in models}
                w_opportunities = {model: f"{opportunities_dir[model]}/chunks/{contig}_{start}_{end}.tsv" for model in models}
                w_fastq         = f"{fastq_dir}/chunks/{contig}_{start}_{end}.fastq"
                for model in models:
                    window_calls[model].append(w_calls[model])
                    window_opportunities[model].append(w_opportunities[model])
                window_fastqs.append(w_fastq)

                chunk_dirs = " ".join(sorted({os.path.dirname(path) for path in
                                              [w_fastq] + list(w_calls.values()) + list(w_opportunities.values())}))
                model_args = " ".join(f"--model {model} {repeat_beds[model]} {w_calls[model]} {w_opportunities[model]}"
                                      for model in models)

                # Step 1 — one pass over the window for every model: call InDels
                # inside the repeat tracts whose size is a multiple of the
                # period, count the opportunities (the denominator of the rate),
                # and write the reads carrying a call, for haploid remapping.
                gwf.target(
                    f"TRScan_{contig}_{start}_{sample}",
                    inputs=[assembly, alignment, coverage] + list(repeat_beds.values()) + tissue_inputs,
                    outputs=list(w_calls.values()) + list(w_opportunities.values()) + [w_fastq],
                    cores=1,
                ) << f"""
                mkdir -p {chunk_dirs}
                python {scripts}/tandem_repeat_scan.py {contig} {alignment} {assembly} {coverage} {tissue} {start} {end} {w_fastq} {model_args}
                """

            # Step 2 — join the windows back into one file per contig.
            #   calls:         concatenated, header kept from the first window
            #   opportunities: counts per (contig, motif, k, tissue, mapQ,
            #                  kmerQ), so they are summed over the windows
            #   FASTQ:         concatenated, a read carrying calls in two
            #                  windows is kept once
            merge_commands = []
            for model in models:
                merge_commands.append(
                    f"awk 'FNR==1 && NR==1 || FNR>1' {' '.join(window_calls[model])} > {calls[model]}")
                merge_commands.append(
                    f"awk -F'\\t' -v OFS='\\t' 'FNR==1 {{if (NR==1) print; next}} "
                    f"{{n[$1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6] += $7}} END {{for (key in n) print key, n[key]}}' "
                    f"{' '.join(window_opportunities[model])} > {opportunities[model]}")
            merge_commands = "\n                ".join(merge_commands)
            contig_dirs = " ".join([fastq_dir] + list(calls_dir.values()) + list(opportunities_dir.values()))

            gwf.target(
                f"TRMergeWindows_{contig}_{sample}",
                inputs=[path for model in models for path in window_calls[model] + window_opportunities[model]] + window_fastqs,
                outputs=list(calls.values()) + list(opportunities.values()) + [fastq],
                cores=1,
                memory="16g",
                walltime="02:00:00",
            ) << f"""
                mkdir -p {contig_dirs}
                {merge_commands}
                cat {' '.join(window_fastqs)} | awk 'NR % 4 == 1 {{keep = !seen[$1]++}} keep' > {fastq}
                """

        fastqs            = [f"{fastq_dir}/{contig}.fastq" for contig in long_contigs]
        fastq_list        = " ".join(fastqs)
        combined_fastq    = f"{haploid_dir}/all.fastq"
        remapped_bam      = f"{haploid_dir}/remapped.bam"
        haploid_locations = f"{haploid_dir}/locations.tsv"

        # Step 3 — remap the reads carrying a call in any model to the haploid
        # assembly, once for all models. The FASTQs are listed explicitly and
        # the combined FASTQ is written to haploid/, so a re-run never reads its
        # own output.
        gwf.target(
            f"TRHaploidMapping_{sample}",
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

        for model in models:
            resultsdir = f"{results}/repetitive_indels/{species}/{sample}/{model}"
            plotsdir   = f"{plots}/repetitive_indels/{species}/{sample}/{model}"

            calls_files          = [f"{calls_dir[model]}/{contig}.tsv" for contig in long_contigs]
            opportunities_files  = [f"{opportunities_dir[model]}/{contig}.tsv" for contig in long_contigs]
            calls_list           = " ".join(calls_files)
            opportunities_list   = " ".join(opportunities_files)
            merged_indels        = f"{merged_dir[model]}/indels.tsv"
            merged_opportunities = f"{merged_dir[model]}/opportunities.tsv"
            plots_pdf            = f"{plotsdir}/plots.pdf"

            # Step 4 — haploid concordance for all contigs of the model: is the
            # same InDel also in the read's alignment to the haploid assembly?
            # locations.tsv is read once for every contig. The contigs are
            # merged into one table (header kept once) with species/sample added.
            gwf.target(
                f"TRHaploidFilter_{sample}_{model}",
                inputs=calls_files + [haploid_locations],
                outputs=[merged_indels],
                cores=1,
                memory="20g",
                walltime="12:00:00",
            ) << f"""
            mkdir -p {merged_dir[model]}
            python {scripts}/filtering_haploid_mappping.py {calls_list} {haploid_locations} > {merged_indels}.tmp
            awk 'NR==1 {{print $0 "\\t" "species" "\\t" "sample"}}; NR>1 {{print $0 "\\t" "{species}" "\\t" "{sample}"}}' {merged_indels}.tmp > {merged_indels}
            rm -f {merged_indels}.tmp
            """

            # Step 5 — merge the contigs' opportunities (header kept from the
            # first file only, species/sample added) and draw diagnostics of
            # the unfiltered calls.
            gwf.target(
                f"TRMergeAndPlot_{sample}_{model}",
                inputs=[merged_indels] + opportunities_files,
                outputs=[merged_opportunities, plots_pdf],
                cores=8,
                memory="200g",
                walltime="24:00:00",
            ) << f"""
            mkdir -p {merged_dir[model]} {plotsdir}
            awk 'FNR==1 && NR==1 || FNR>1' {opportunities_list} \\
                | awk 'NR==1 {{print $0 "\\t" "species" "\\t" "sample"}}; NR>1 {{print $0 "\\t" "{species}" "\\t" "{sample}"}}' \\
                > {merged_opportunities}
            Rscript {scripts}/plot_indels.R {merged_indels} {merged_opportunities} {plots_pdf}
            """

            # Step 6 — final filtering, with a report of rows removed per step.
            # These small tables are the results.
            gwf.target(
                f"TRFinalFiltering_{sample}_{model}",
                inputs=[merged_indels, merged_opportunities],
                outputs=[f"{resultsdir}/indels.tsv",
                         f"{resultsdir}/opportunities.tsv",
                         f"{resultsdir}/filter_report.tsv"],
                cores=4,
                memory="200g",
                walltime="01:00:00",
            ) << f"""
            mkdir -p {resultsdir}
            python {scripts}/filtering_small_indels.py {merged_indels} \\
                -o {resultsdir}/indels.tsv \\
                --haploid-check True \\
                --min-quality {min_kmerq} \\
                --min-read-end-distance {min_read_end_distance} \\
                --report {resultsdir}/filter_report.tsv
            python {scripts}/filtering_opportunities.py {merged_opportunities} \\
                -o {resultsdir}/opportunities.tsv \\
                --min-kmerq {min_kmerq}
            """

            # Step 7 — cell-type classification, for the samples listed under
            # methylation_samples in config.yaml. Every read of those samples
            # has a log-likelihood per cell type; the script contrasts the
            # germline (Sperm.*) columns with a set of somatic ones and appends
            # the difference and its bin to every call. The likelihood table
            # covers the whole sample, so it is streamed once and filtered to
            # the reads that carry a call.
            if sample in methylation_dict.get(species, []):
                likelihoods = (f"{inpath_classification}/{species}/{sample}/{sample}"
                               f"/05_ReadsClassification/{sample}_All_Reads_Likelihood_merged.tsv")

                gwf.target(
                    f"TRClassify_{sample}_{model}",
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
