from functools import partial

from gwf import Workflow

from pipeline import (Pipeline, final_filtering_options, grouped_by_length, joined_calls, joined_reads, remapping,
                      retried, shell)

# ---------------------------------------------------------------------------
# InDels OUTSIDE repetitive DNA. Calls 1-3 bp insertions and deletions in single
# reads at positions that are not masked. The mask is the union of every repeat
# model from repetitive_dna.py (M0, M1, M2), its excluded group (homopolymer
# runs and tracts with longer periods) and the unresolved-haplotype regions
# from precompute.py.
#
# Homopolymer runs shorter than the excluded ones are not in the mask. They are
# handled inside the caller instead: an InDel whose sequence plus one flanking
# base on each side changes base at most once is dropped, and the opportunity
# counter applies the same rule, so numerator and denominator agree.
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
# Paths come from the `layout` templates and the parameters from config.yaml,
# through pipeline.py.
# ---------------------------------------------------------------------------

pipeline          = Pipeline()
config            = pipeline.config
assembly_type     = config["assembly_type"]
scripts           = pipeline.scripts

# Repeat models whose tracts are masked out, from repetitive_dna.py.
models = config["repeats"]["models"]

# Step 2 splits each contig into windows of this many bp, and runs this many
# windows in one job, at the same time, one core each: fewer jobs, each as long
# as its slowest window.
chunk_size = config["scan"]["window_size"]
windows_per_job = config["scan"]["windows_per_job"]
# Step 3 merges the windows of contigs adding up to this many bp in one job.
merge_bp_per_job = config["scan"]["merge_bp_per_job"]
# The positions scanned, and the InDel sizes called and counted.
scan_options = (f"--coverage-band {' '.join(map(str, config['scan']['coverage_band']))} "
                f"--indel-sizes {' '.join(map(str, config['scan']['indel_sizes']))} "
                f"--flank-width {config['scan']['flank_width']} "
                f"--key-parameters {pipeline.key_parameters(nonrepetitive=True)}")



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
        resultsdir  = folder("nonrepetitive_results_calls")
        plotsdir    = folder("nonrepetitive_plots")
        long_contigs = pipeline.long_contigs(species, sample)

        # Samples without read labels have no tissue table; the callers then
        # report tissue NA for every read.
        if pipeline.is_listed("tissue_samples", species, sample):
            tissue        = path("tissue")
            tissue_inputs = [tissue]
        else:
            tissue        = "NA"
            tissue_inputs = []
        # Samples with methylation scores: each read's germline/somatic bin,
        # from the compact table (ticket 44).
        if pipeline.is_listed("methylation_samples", species, sample):
            read_bins_inputs      = [path("read_bins")]
            sample_scan_options   = f"{scan_options} --read-bins {read_bins_inputs[0]}"
        else:
            read_bins_inputs      = []
            sample_scan_options   = scan_options
        # The scan keeps a count per combination of keys until its window is
        # done. The bins split every count, and on Revio reads the subread
        # support's bins do too: ob007's scans needed up to 20 GB and more,
        # every other sample's 6 GB at most (ticket 45). A tissue sample's
        # tissues (ten for R22013) split them as well, and each scan loads the
        # sample's read-to-tissue table (12 M reads, about 2.4 GB): R22013's
        # scans peaked at 35-41 GB of 40, and 1,376 of the tissue samples'
        # failed there (ticket 51). This is per window; a job asks for it once
        # for each of its windows.
        many_keys   = read_bins_inputs and not pipeline.is_listed("sequel_ii_samples", species, sample)
        window_memory_gb = 40 if many_keys else 20 if tissue_inputs else 10

        mask_dir          = folder("nonrepetitive_mask")
        calls_dir         = folder("nonrepetitive_calls")
        opportunities_dir = folder("nonrepetitive_opportunities")
        fastq_dir         = folder("nonrepetitive_reads")
        haploid_dir       = folder("nonrepetitive_haploid_locations")
        merged_dir        = folder("nonrepetitive_merged_calls")

        masks         = {contig: path("nonrepetitive_mask", contig=contig) for contig in long_contigs}
        calls         = {contig: path("nonrepetitive_calls", contig=contig) for contig in long_contigs}
        opportunities = {contig: path("nonrepetitive_opportunities", contig=contig) for contig in long_contigs}
        fastqs        = {contig: path("nonrepetitive_reads", contig=contig) for contig in long_contigs}

        # Step 1 — the mask of every contig: every base any repeat model calls
        # repetitive, the unresolved-haplotype regions, and the excluded group
        # (homopolymer runs and tracts with longer periods), merged. The
        # repeat bases are those of every tract before overlap resolution, so
        # a tract resolution drops is masked whole.
        repeat_beds     = [path("repeat_bases", contig=contig, model=model) for contig in long_contigs for model in models]
        unresolved_beds = [path("unresolved", contig=contig) for contig in long_contigs]
        excluded_beds   = [path("excluded", contig=contig) for contig in long_contigs]
        # The same files for the contig in the shell loop's $contig.
        contig_repeat_beds = " ".join(path("repeat_bases", contig="${contig}", model=model) for model in models)
        contig_unresolved  = path("unresolved", contig="${contig}")
        contig_excluded    = path("excluded", contig="${contig}")
        contig_mask        = path("nonrepetitive_mask", contig="${contig}")
        contig_names       = " ".join(long_contigs)

        gwf.target(
            f"non_rep_Mask_{species}_{sample}_{assembly_type}",
            inputs=repeat_beds + unresolved_beds + excluded_beds,
            outputs=list(masks.values()),
            cores=1,
            memory="8g",
            walltime="04:00:00",
        ) << f"""
        set -e
        mkdir -p {mask_dir}
        for contig in {contig_names}; do
            python {scripts}/merge_regions.py {contig_repeat_beds} {contig_unresolved} {contig_excluded} > {contig_mask}.tmp
            mv {contig_mask}.tmp {contig_mask}
        done
        """

        # Step 2 runs per window. The windows [0, chunk_size),
        # [chunk_size, 2*chunk_size), ... tile each contig, and every pileup
        # position belongs to exactly one window, so each call and each
        # opportunity is counted once.
        window_calls         = {contig: [] for contig in long_contigs}
        window_opportunities = {contig: [] for contig in long_contigs}
        window_fastqs        = {contig: [] for contig in long_contigs}
        windows              = []

        for contig, length in long_contigs.items():
            for start in range(0, length, chunk_size):
                end = min(start + chunk_size, length)
                w_calls         = path("nonrepetitive_window_calls", contig=contig, start=start, end=end)
                w_opportunities = path("nonrepetitive_window_opportunities", contig=contig, start=start, end=end)
                w_fastq         = path("nonrepetitive_window_reads", contig=contig, start=start, end=end)
                window_calls[contig].append(w_calls)
                window_opportunities[contig].append(w_opportunities)
                window_fastqs[contig].append(w_fastq)
                windows.append((contig, start, end, w_calls, w_opportunities, w_fastq))

        # Step 2 — one pass over each window: call 1-3 bp InDels outside the
        # mask, dropping those in a homopolymer context; count the
        # opportunities (the denominator of the rate) with the same homopolymer
        # rule; and write the reads carrying a call, for haploid remapping.
        # windows_per_job windows per job, in the order above, scanned at the
        # same time; the job fails if any of them fails.
        for first in range(0, len(windows), windows_per_job):
            group = windows[first:first + windows_per_job]
            scans = "\n".join(
                f"python {scripts}/nonrepetitive_scan.py {masks[contig]} {contig} {alignment} {assembly} {coverage} "
                f"{tissue} {start} {end} {w_calls} {w_opportunities} {w_fastq} {sample_scan_options} & pids+=($!)"
                for contig, start, end, w_calls, w_opportunities, w_fastq in group)
            group_masks = list(dict.fromkeys(masks[contig] for contig, *_ in group))
            contig, start = group[0][:2]

            gwf.target(
                f"non_rep_Scan_{contig}_{start}_{species}_{sample}_{assembly_type}",
                inputs=[assembly, alignment, coverage] + group_masks + tissue_inputs + read_bins_inputs,
                outputs=[path for window in group for path in window[3:]],
                cores=len(group),
                memory=f"{window_memory_gb * len(group)}g",
                walltime="04:00:00",
            ) << f"""
set -euo pipefail
mkdir -p {folder("nonrepetitive_window_calls")} {folder("nonrepetitive_window_opportunities")} {folder("nonrepetitive_window_reads")}
pids=()
{scans}
for pid in "${{pids[@]}}"; do wait "$pid"; done
"""

        # Step 3 — join the windows back into one file per contig (ticket 37),
        # for the contigs of a group (merge_bp_per_job of them) in one job,
        # one after the other:
        #   calls:         concatenated, header kept from the first window
        #   FASTQ:         concatenated, a read carrying calls in two windows is
        #                  kept once
        # 3b — the opportunities, counts per contig, context, mapQ and the keys
        # of strata.py, summed over the windows into the cube's Parquet (ticket
        # 40), in jobs of their own: the calls and reads, and the remapping
        # after them, don't depend on how the opportunities are kept. A sum
        # that fails its check at random (ticket 41) is retried in the job, so
        # the group's other contigs aren't summed again.
        for group in grouped_by_length(long_contigs, merge_bp_per_job):
            first = group[0]
            joins = "\n".join(f"{joined_calls(window_calls[contig], calls[contig])}\n"
                              f"{joined_reads(window_fastqs[contig], fastqs[contig])}" for contig in group)
            gwf.target(
                f"non_rep_MergeWindows_{first}_{species}_{sample}_{assembly_type}",
                inputs=[path for contig in group for path in window_calls[contig] + window_fastqs[contig]],
                outputs=[path for contig in group for path in (calls[contig], fastqs[contig])],
                cores=1,
                memory="4g",
                walltime="04:00:00",
            ) << f"""
set -e
mkdir -p {calls_dir} {fastq_dir}
{joins}
"""

            sums = "\n".join(retried(f"python {scripts}/opportunity_tables.py sum {opportunities[contig]} "
                                     f"{' '.join(window_opportunities[contig])}") for contig in group)
            gwf.target(
                f"non_rep_MergeOpportunities_{first}_{species}_{sample}_{assembly_type}",
                inputs=[path for contig in group for path in window_opportunities[contig]],
                outputs=[opportunities[contig] for contig in group],
                cores=4,   # sort and DuckDB use them
                memory="32g",   # sort, then DuckDB, gets half (ticket 50)
                walltime="04:00:00",
            ) << f"""
set -e
mkdir -p {opportunities_dir}
{sums}
"""

        fastq_list        = " ".join(fastqs.values())
        combined_fastq    = path("nonrepetitive_haploid_reads")
        remapped_bam      = path("nonrepetitive_haploid_bam")
        haploid_locations = path("nonrepetitive_haploid_locations")

        # Step 4 — remap all call-carrying reads to the haploid assembly. The
        # FASTQs are listed explicitly and the combined FASTQ is written to
        # haploid/, so a re-run never reads its own output.
        gwf.target(
            f"non_rep_HaploidMapping_{species}_{sample}_{assembly_type}",
            inputs=list(fastqs.values()) + [haploid],
            outputs=[combined_fastq, remapped_bam, haploid_locations],
            cores=16,
            memory="100g",
            walltime="08:00:00",
        ) << f"""
        mkdir -p {haploid_dir}
        cat {fastq_list} > {combined_fastq}
        {remapping(haploid, combined_fastq, remapped_bam, haploid_locations)}
        """

        # Step 4b — the same reads aligned to each haplotype of the diploid
        # assembly alone, for the hap1/hap2 test.
        haplotype_locations = {}
        for haplotype in pipeline.haplotypes:
            reference = path("haplotype_assembly", haplotype=haplotype)
            haplotype_bam = path("nonrepetitive_haplotype_bam", haplotype=haplotype)
            haplotype_locations[haplotype] = path("nonrepetitive_haplotype_locations", haplotype=haplotype)
            gwf.target(
                f"non_rep_HaplotypeMapping_{haplotype}_{species}_{sample}_{assembly_type}",
                inputs=[combined_fastq, reference],
                outputs=[haplotype_bam, haplotype_locations[haplotype]],
                cores=16,
                memory="100g",
                walltime="08:00:00",
            ) << f"""
            {remapping(reference, combined_fastq, haplotype_bam, haplotype_locations[haplotype])}
            """
        haplotype_options = " ".join(f"--haplotype {haplotype}={locations}"
                                     for haplotype, locations in haplotype_locations.items())

        calls_list           = " ".join(calls.values())
        merged_indels        = path("nonrepetitive_merged_calls")
        merged_opportunities = path("nonrepetitive_merged_opportunities")
        plots_pdf            = path("nonrepetitive_plots")

        # Step 5 — haploid concordance for all contigs: is the same InDel also
        # in the read's alignment to the haploid assembly? locations.tsv and the
        # reads, whose sequence says how far each InDel can slide, are read
        # once for every contig. The contigs are merged into one table (header
        # kept once) with species/sample added.
        gwf.target(
            f"non_rep_HaploidFilter_{species}_{sample}_{assembly_type}",
            inputs=list(calls.values()) + [haploid_locations, combined_fastq] + list(haplotype_locations.values()),
            outputs=[merged_indels],
            cores=1,
            memory="20g",
            walltime="04:00:00",
        ) << f"""
        mkdir -p {merged_dir}
        python {scripts}/filtering_haploid_mapping.py {calls_list} {haploid_locations} {combined_fastq} {haplotype_options} > {merged_indels}.tmp
        awk 'NR==1 {{print $0 "\\t" "species" "\\t" "sample"}}; NR>1 {{print $0 "\\t" "{species}" "\\t" "{sample}"}}' {merged_indels}.tmp > {merged_indels}
        rm -f {merged_indels}.tmp
        """

        # Step 6 — join the contigs' opportunity cubes into the sample's, with
        # species/sample added, and draw diagnostics of the unfiltered calls.
        gwf.target(
            f"non_rep_MergeAndPlot_{species}_{sample}_{assembly_type}",
            inputs=[merged_indels] + list(opportunities.values()),
            outputs=[merged_opportunities, plots_pdf],
            cores=8,
            memory="64g",   # most of it R reading the calls
            walltime="04:00:00",
        ) << f"""
        mkdir -p {merged_dir} {plotsdir}
        python {scripts}/opportunity_tables.py join {merged_opportunities} --species {species} --sample {sample} \\
            {' '.join(opportunities.values())}
        Rscript {scripts}/plot_indels.R {merged_indels} {plots_pdf}
        """

        # Step 7 — final filtering, with a report of the rows each step removes.
        # These small tables are the results. kmerQ by position and the base
        # after a deletion apply to calls and opportunities alike; the
        # inserted bases and the haplotype checks to the calls only (ticket 49,
        # final_filters in config.yaml). The report counts the calls without a
        # haploid check.
        results_calls         = path("nonrepetitive_results_calls")
        results_opportunities = path("nonrepetitive_results_opportunities")
        results_report        = path("nonrepetitive_results_report")

        gwf.target(
            f"non_rep_FinalFiltering_{species}_{sample}_{assembly_type}",
            inputs=[merged_indels, merged_opportunities, path("smrt_cells")],
            outputs=[results_calls,
                     results_opportunities,
                     results_report],
            cores=4,
            memory="16g",
            walltime="04:00:00",
        ) << f"""
        mkdir -p {resultsdir}
        python {scripts}/final_filtering.py {merged_indels} {merged_opportunities} \\
            {results_calls} {results_opportunities} {results_report} \\
            --smrt-cells {path("smrt_cells")} {shell(final_filtering_options(config, species, sample))}
        """

        # Step 8 — cell-type classification, for the samples listed under
        # methylation_samples in config.yaml. Every read of those samples has a
        # log-likelihood per cell type; the script contrasts the germline
        # (Sperm.*) columns with a set of somatic ones and appends the
        # difference and its bin to every call. The likelihood table covers the
        # whole sample, so it is streamed once and filtered to the reads that
        # carry a call.
        if pipeline.is_listed("methylation_samples", species, sample):
            likelihoods = path("likelihoods")
            classified  = path("nonrepetitive_results_classified")

            gwf.target(
                f"non_rep_Classify_{species}_{sample}_{assembly_type}",
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
