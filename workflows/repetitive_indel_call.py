import os
from functools import partial

from gwf import Workflow

from pipeline import Pipeline, summed_counts

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
# Paths come from the `layout` templates and the parameters from config.yaml,
# through pipeline.py.
# ---------------------------------------------------------------------------

pipeline          = Pipeline()
config            = pipeline.config
assembly_type     = config["assembly_type"]
scripts           = pipeline.scripts

# Repeat models written by repetitive_dna.py.
models = config["repeats"]["models"]

# Step 1 splits each contig into windows of this many bp, one scan job per window.
chunk_size = config["scan"]["window_size"]
# The positions scanned.
scan_options = (f"--coverage-band {' '.join(map(str, config['scan']['coverage_band']))} "
                f"--flank-width {config['scan']['flank_width']} "
                f"--key-parameters {pipeline.key_parameters()}")

# Final filtering thresholds.
min_read_end_distance = config["final_filters"]["min_read_end_distance"]


# Resource requests.
gwf = Workflow(
    defaults={
        'cores': 4,
        'memory': '10g',
        'walltime': '12:00:00',
        'account': 'mutationalscanning'})


for species, samples in config["samples"].items():
    for sample in samples:
        path        = partial(pipeline.path, species=species, sample=sample)
        folder      = partial(pipeline.folder, species=species, sample=sample)
        assembly    = path("assembly")
        haploid     = path("haploid")
        # Samples listed under alignment_overrides in config.yaml use that BAM.
        alignment   = pipeline.alignment(species, sample)
        coverage    = path("coverage")
        long_contigs = pipeline.long_contigs(species, sample)

        # Samples without read labels have no tissue table; the callers then
        # report tissue NA for every read.
        if pipeline.is_listed("tissue_samples", species, sample):
            tissue        = path("tissue")
            tissue_inputs = [tissue]
        else:
            tissue        = "NA"
            tissue_inputs = []
        # Samples with methylation scores: each read's germline/somatic bin.
        if pipeline.is_listed("methylation_samples", species, sample):
            classification_inputs = [path("read_classification")]
            sample_scan_options   = f"{scan_options} --read-classification {classification_inputs[0]}"
        else:
            classification_inputs = []
            sample_scan_options   = scan_options

        # Shared by the models: the reads carrying calls, and their remapping.
        fastq_dir   = folder("repetitive_reads")
        haploid_dir = folder("repetitive_haploid_locations")

        # Per model.
        calls_dir         = {model: folder("repetitive_calls", model=model) for model in models}
        opportunities_dir = {model: folder("repetitive_opportunities", model=model) for model in models}
        merged_dir        = {model: folder("repetitive_merged_calls", model=model) for model in models}

        for contig, length in long_contigs.items():
            repeat_beds   = {model: path("tandem_repeats", contig=contig, model=model) for model in models}
            calls         = {model: path("repetitive_calls", contig=contig, model=model) for model in models}
            opportunities = {model: path("repetitive_opportunities", contig=contig, model=model) for model in models}
            fastq         = path("repetitive_reads", contig=contig)

            # Step 1 runs per window. The windows [0, chunk_size),
            # [chunk_size, 2*chunk_size), ... tile the contig, and every pileup
            # position belongs to exactly one window, so each call and each
            # opportunity is counted once.
            window_calls         = {model: [] for model in models}
            window_opportunities = {model: [] for model in models}
            window_fastqs        = []

            for start in range(0, length, chunk_size):
                end = min(start + chunk_size, length)
                w_calls         = {model: path("repetitive_window_calls", contig=contig, start=start, end=end, model=model)
                                   for model in models}
                w_opportunities = {model: path("repetitive_window_opportunities", contig=contig, start=start, end=end,
                                               model=model) for model in models}
                w_fastq         = path("repetitive_window_reads", contig=contig, start=start, end=end)
                for model in models:
                    window_calls[model].append(w_calls[model])
                    window_opportunities[model].append(w_opportunities[model])
                window_fastqs.append(w_fastq)

                chunk_dirs = " ".join(sorted({os.path.dirname(file) for file in
                                              [w_fastq] + list(w_calls.values()) + list(w_opportunities.values())}))
                model_args = " ".join(f"--model {model} {repeat_beds[model]} {w_calls[model]} {w_opportunities[model]}"
                                      for model in models)

                # Step 1 — one pass over the window for every model: call InDels
                # inside the repeat tracts whose size is a multiple of the
                # period, count the opportunities (the denominator of the rate),
                # and write the reads carrying a call, for haploid remapping.
                gwf.target(
                    f"TRScan_{contig}_{start}_{species}_{sample}",
                    inputs=[assembly, alignment, coverage] + list(repeat_beds.values()) + tissue_inputs + classification_inputs,
                    outputs=list(w_calls.values()) + list(w_opportunities.values()) + [w_fastq],
                    cores=1,
                ) << f"""
                mkdir -p {chunk_dirs}
                python {scripts}/tandem_repeat_scan.py {contig} {alignment} {assembly} {coverage} {tissue} {start} {end} {w_fastq} {sample_scan_options} {model_args}
                """

            # Step 2 — join the windows back into one file per contig.
            #   calls:         concatenated, header kept from the first window
            #   opportunities: counts per contig, motif, k, mapQ and the keys
            #                  of strata.py, so they are summed over the windows
            #   FASTQ:         concatenated, a read carrying calls in two
            #                  windows is kept once
            merge_commands = []
            for model in models:
                merge_commands.append(
                    f"awk 'FNR==1 && NR==1 || FNR>1' {' '.join(window_calls[model])} > {calls[model]}")
                merge_commands.append(summed_counts(window_opportunities[model], opportunities[model]))
            merge_commands = "\n                ".join(merge_commands)
            contig_dirs = " ".join([fastq_dir] + list(calls_dir.values()) + list(opportunities_dir.values()))

            gwf.target(
                f"TRMergeWindows_{contig}_{species}_{sample}",
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

        fastqs            = [path("repetitive_reads", contig=contig) for contig in long_contigs]
        fastq_list        = " ".join(fastqs)
        combined_fastq    = path("repetitive_haploid_reads")
        remapped_bam      = path("repetitive_haploid_bam")
        haploid_locations = path("repetitive_haploid_locations")

        # Step 3 — remap the reads carrying a call in any model to the haploid
        # assembly, once for all models. The FASTQs are listed explicitly and
        # the combined FASTQ is written to haploid/, so a re-run never reads its
        # own output.
        gwf.target(
            f"TRHaploidMapping_{species}_{sample}",
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
            resultsdir = folder("repetitive_results_calls", model=model)
            plotsdir   = folder("repetitive_plots", model=model)

            calls_files          = [path("repetitive_calls", contig=contig, model=model) for contig in long_contigs]
            opportunities_files  = [path("repetitive_opportunities", contig=contig, model=model) for contig in long_contigs]
            calls_list           = " ".join(calls_files)
            opportunities_list   = " ".join(opportunities_files)
            merged_indels        = path("repetitive_merged_calls", model=model)
            merged_opportunities = path("repetitive_merged_opportunities", model=model)
            plots_pdf            = path("repetitive_plots", model=model)

            # Step 4 — haploid concordance for all contigs of the model: is the
            # same InDel also in the read's alignment to the haploid assembly?
            # locations.tsv and the reads, whose sequence says how far each
            # InDel can slide, are read once for every contig. The contigs are
            # merged into one table (header kept once) with species/sample added.
            gwf.target(
                f"TRHaploidFilter_{species}_{sample}_{model}",
                inputs=calls_files + [haploid_locations, combined_fastq],
                outputs=[merged_indels],
                cores=1,
                memory="20g",
                walltime="12:00:00",
            ) << f"""
            mkdir -p {merged_dir[model]}
            python {scripts}/filtering_haploid_mappping.py {calls_list} {haploid_locations} {combined_fastq} > {merged_indels}.tmp
            awk 'NR==1 {{print $0 "\\t" "species" "\\t" "sample"}}; NR>1 {{print $0 "\\t" "{species}" "\\t" "{sample}"}}' {merged_indels}.tmp > {merged_indels}
            rm -f {merged_indels}.tmp
            """

            # Step 5 — merge the contigs' opportunities (header kept from the
            # first file only, species/sample added) and draw diagnostics of
            # the unfiltered calls.
            gwf.target(
                f"TRMergeAndPlot_{species}_{sample}_{model}",
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
            # These small tables are the results. The haploid check is a
            # column of the calls, not a filter; the report counts the calls
            # without one.
            results_calls         = path("repetitive_results_calls", model=model)
            results_opportunities = path("repetitive_results_opportunities", model=model)
            results_report        = path("repetitive_results_report", model=model)

            gwf.target(
                f"TRFinalFiltering_{species}_{sample}_{model}",
                inputs=[merged_indels, merged_opportunities],
                outputs=[results_calls,
                         results_opportunities,
                         results_report],
                cores=4,
                memory="200g",
                walltime="01:00:00",
            ) << f"""
            mkdir -p {resultsdir}
            python {scripts}/filtering_small_indels.py {merged_indels} \\
                -o {results_calls} \\
                --keep kmerQ_pass=pass \\
                --min-read-end-distance {min_read_end_distance} \\
                --report {results_report}
            python {scripts}/filtering_opportunities.py {merged_opportunities} \\
                -o {results_opportunities} \\
                --keep kmerQ_pass=pass
            """

            # Step 7 — cell-type classification, for the samples listed under
            # methylation_samples in config.yaml. Every read of those samples
            # has a log-likelihood per cell type; the script contrasts the
            # germline (Sperm.*) columns with a set of somatic ones and appends
            # the difference and its bin to every call. The likelihood table
            # covers the whole sample, so it is streamed once and filtered to
            # the reads that carry a call.
            if pipeline.is_listed("methylation_samples", species, sample):
                likelihoods = path("likelihoods")
                classified  = path("repetitive_results_classified", model=model)

                gwf.target(
                    f"TRClassify_{species}_{sample}_{model}",
                    inputs=[results_calls, likelihoods],
                    outputs=[classified],
                    cores=1,
                    memory="8g",
                    walltime="04:00:00",
                ) << f"""
                set -e
                bash {scripts}/join_read_classification.sh \\
                    {results_calls} \\
                    {likelihoods} \\
                    {classified} \\
                    {pipeline.classification_arguments()}
                """
