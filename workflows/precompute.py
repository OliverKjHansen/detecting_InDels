from gwf import Workflow

from pipeline import Pipeline

# ---------------------------------------------------------------------------
# Configuration, read from config.yaml. Run gwf from the repository root:
#
#   gwf -f workflows/precompute.py run
#
# Paths come from the `layout` templates in config.yaml, through pipeline.py.
# ---------------------------------------------------------------------------

pipeline          = Pipeline()
assembly_type     = pipeline.config["assembly_type"]
min_contig_length = pipeline.config["min_contig_length"]
scripts           = pipeline.scripts

# Step 3 finds unresolved haplotypes in windows of this many bp, one job each.
chunk_size = pipeline.config["precompute"]["unresolved_window_size"]
# What counts as evidence of an unresolved region.
unresolved_evidence = (f"--min-base-quality {pipeline.config['precompute']['unresolved_min_base_quality']} "
                       f"--min-reads {pipeline.config['precompute']['unresolved_min_reads']}")


# Resource requests.
gwf = Workflow(
    defaults={
        'cores': 4,
        'memory': '40g',
        'walltime': '12:00:00',
        'account': 'mutationalscanning'})


for species, samples in pipeline.config["samples"].items():
    for sample in samples:
        assembly  = pipeline.path("assembly", species, sample)
        # Samples listed under alignment_overrides in config.yaml use that BAM.
        alignment = pipeline.alignment(species, sample)
        outdir    = pipeline.folder("coverage", species, sample)
        coverage  = pipeline.path("coverage", species, sample)
        long_contigs = pipeline.long_contigs(species, sample)

        # Step 1 — split the assembly into per-contig FASTAs.
        contigs_dir = pipeline.folder("contig_fasta", species, sample)

        gwf.target(
            f"faSplit_{species}_{sample}_{assembly_type}",
            inputs=[assembly],
            outputs=[pipeline.path("contig_fasta", species, sample, contig=contig) for contig in long_contigs],
        ) << f"""
        mkdir -p {contigs_dir}
        faSplit byname {assembly} {contigs_dir}/
        """

        # Step 2 — average depth across the positions samtools reports.
        long_contigs_bed = pipeline.path("long_contigs", species, sample)

        gwf.target(
            f"average_coverage_{species}_{sample}",
            inputs=[alignment, pipeline.path("assembly_fai", species, sample)],
            outputs=[long_contigs_bed, coverage],
            cores=8,
            memory="100g",
            walltime="08:00:00",
        ) << f"""
        mkdir -p {outdir}
        awk -v min={min_contig_length} 'BEGIN {{OFS="\\t"}} $2 > min {{print $1, 0, $2}}' {assembly}.fai > {long_contigs_bed}
        samtools depth -b {long_contigs_bed} {alignment} | awk '{{sum += $3; count++}} END {{print sum / count}}' > {coverage}
        """

        # Step 3 — unresolved haplotypes, per contig, split into windows.
        #
        # The windows [0, chunk_size), [chunk_size, 2*chunk_size), ... tile the
        # contig, so every position is scanned by exactly one job. A region is
        # built from the full span of its reads, so it can reach past its
        # window; the pieces from neighbouring windows then overlap, and
        # bedtools merge joins them back into one region.
        chunk_dir = pipeline.folder("unresolved_window", species, sample)

        for contig, length in long_contigs.items():
            unresolved = pipeline.path("unresolved", species, sample, contig=contig)
            chunk_beds = []

            # Step 3a — one job per window.
            for start in range(0, length, chunk_size):
                end = min(start + chunk_size, length)
                chunk_bed = pipeline.path("unresolved_window", species, sample, contig=contig, start=start, end=end)
                chunk_beds.append(chunk_bed)

                gwf.target(
                    f"unresolved_haplotype_{species}_{sample}_{contig}_{start}_{assembly_type}",
                    inputs=[assembly, alignment, coverage],
                    outputs=[chunk_bed],
                ) << f"""
                mkdir -p {chunk_dir}
                python {scripts}/resolve_haplotypes.py {contig} {start} {end} {alignment} {assembly} {coverage} {unresolved_evidence} > {chunk_bed}
                """

            # Step 3b — join the windows back into one BED per contig.
            chunk_list = " ".join(chunk_beds)

            gwf.target(
                f"merge_unresolved_{species}_{sample}_{contig}_{assembly_type}",
                inputs=chunk_beds,
                outputs=[unresolved],
                cores=1,
                memory="4g",
                walltime="00:30:00",
            ) << f"""
            cat {chunk_list} | sort -k1,1 -k2,2n | bedtools merge -i - > {unresolved}
            """

        # Step 4 — read -> tissue table, sorted by read ID. Only samples listed
        # under tissue_samples in config.yaml have read labels.
        if pipeline.is_listed("tissue_samples", species, sample):
            read_labels = pipeline.path("read_labels", species, sample)
            tissue      = pipeline.path("tissue", species, sample)

            gwf.target(
                f"sorting_tissue_classification_{species}_{sample}_{assembly_type}",
                inputs=[read_labels],
                outputs=[tissue],
            ) << f"""
            mkdir -p {pipeline.folder("tissue", species, sample)}
            printf "readID\\ttissue\\n" > {tissue}
            awk -v OFS="\\t" '{{gsub(/[ \\t]+/, "\\t"); print $1, $2}}' {read_labels} | sort -k1,1 >> {tissue}
            """

        # Step 5 — every read's germline/somatic bin, for the samples listed
        # under methylation_samples in config.yaml. The scans read each read's
        # bin from it, for calls and opportunities alike.
        if pipeline.is_listed("methylation_samples", species, sample):
            likelihoods    = pipeline.path("likelihoods", species, sample)
            classification = pipeline.path("read_classification", species, sample)

            gwf.target(
                f"read_classification_{species}_{sample}_{assembly_type}",
                inputs=[likelihoods],
                outputs=[classification],
                cores=1,
                memory="4g",
                walltime="06:00:00",
            ) << f"""
            mkdir -p {pipeline.folder("read_classification", species, sample)}
            bash {scripts}/read_classification.sh {likelihoods} {classification} {pipeline.classification_arguments()}
            """
