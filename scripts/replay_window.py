"""
Replay workflows 3 and 4 on a window mirrored from the cluster, with this
checkout's scripts, and compare each table with the cluster's. For software
testing only: no data conclusions are drawn from the mirror.

    python scripts/replay_window.py mirror/human/chk/h1tg000001l_0_10000000

The window's config.yaml says which sample it holds and where its files are.
The commands are built here, with the parameters in this checkout's
config.yaml, so a window exported before a script changed its arguments still
replays, and the replay tests the pipeline as it is now. Outputs go to local/ in
the window and the scripts' messages to local/replay.log. Each table is
compared with scripts/compare_tables.py; the exit status is 1 if any differs.
"""
import json
import os
import subprocess
import sys

import yaml

from compare_tables import compare, report

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(SCRIPTS)
sys.path.insert(0, os.path.join(REPO, "workflows"))
from pipeline import Pipeline, final_filtering_options  # noqa: E402


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    window = os.path.abspath(sys.argv[1])
    contig, start, end = os.path.basename(window).rsplit("_", 2)
    with open(os.path.join(REPO, "config.yaml")) as handle:
        current = yaml.safe_load(handle)
    models = current["repeats"]["models"]
    scan_options = ["--coverage-band", *current["scan"]["coverage_band"],
                    "--flank-width", current["scan"]["flank_width"],
                    "--key-parameters", json.dumps(current["strata"])]
    data = Pipeline(f"{window}/config.yaml", base=window)
    # A window exported before a layout template existed lacks it; this
    # checkout's template places the file the same way.
    data.config["layout"] = {**current["layout"], **data.config["layout"]}
    (species, samples), = data.config["samples"].items()
    sample, = samples
    tissue = data.path("tissue", species, sample) if data.is_listed("tissue_samples", species, sample) else "NA"
    inputs = [contig, data.alignment(species, sample), data.path("assembly", species, sample),
              data.path("coverage", species, sample), tissue, start, end]
    # The reads' germline/somatic bins, if the window has them.
    classification = data.path("read_classification", species, sample)
    if os.path.exists(classification):
        scan_options += ["--read-classification", classification]

    local, cluster = f"{window}/local", f"{window}/cluster"
    os.makedirs(local, exist_ok=True)
    log = open(f"{local}/replay.log", "w")

    def run(script, *args, stdout=None):
        subprocess.run([sys.executable, os.path.join(SCRIPTS, script), *map(str, args)],
                       stdout=stdout or log, stderr=log, check=True)

    def replay(workflow, calls, reads, results):
        """
        The haploid check and the final filters, from a scan's calls,
        opportunities and reads (the FASTQ the workflows remap), with the
        hap1/hap2 test if the window has the reads' alignments to each
        haplotype.
        """
        haplotypes = []
        for haplotype in current["precompute"]["haplotype_prefixes"]:
            records = data.path(f"{workflow}_haplotype_locations", species, sample, haplotype=haplotype)
            if os.path.exists(records):
                haplotypes += ["--haplotype", f"{haplotype}={records}"]
        with open(f"{results}/merged_indels.tsv", "w") as merged:
            run("filtering_haploid_mapping.py", calls, data.path(f"{workflow}_haploid_locations", species, sample),
                reads, *haplotypes, stdout=merged)
        run("final_filtering.py", f"{results}/merged_indels.tsv", f"{results}/opportunities.tsv",
            f"{results}/indels.tsv", f"{results}/final_opportunities.tsv", f"{results}/filter_report.tsv",
            "--smrt-cells", smrt_cells,
            *final_filtering_options(current, species, sample, repeats=workflow == "repetitive"))

    # Each SMRT Cell's highest base quality, from the window's reads (workflow 1, step 2b).
    smrt_cells = f"{local}/smrt_cells.tsv"
    run("smrt_cell_qualities.py", data.alignment(species, sample), smrt_cells,
        "--reads", current["precompute"]["smrt_cell_reads"])

    # Workflow 3: the scan (step 2), the haploid check (5), final filtering (7).
    w3 = f"{local}/nonrepetitive_indels"
    os.makedirs(w3, exist_ok=True)
    run("nonrepetitive_scan.py", data.path("nonrepetitive_mask", species, sample, contig=contig), *inputs,
        f"{w3}/calls.tsv", f"{w3}/opportunities.tsv", f"{w3}/reads.fastq",
        *scan_options, "--indel-sizes", *current["scan"]["indel_sizes"])
    replay("nonrepetitive", f"{w3}/calls.tsv", f"{w3}/reads.fastq", w3)
    tables = [("workflow 3", w3, f"{cluster}/nonrepetitive_indels")]

    def tracts(model):
        """
        The tracts workflow 4 scans: those it measures, or, in a window exported
        before it dropped the tracts that overlap an excluded or unresolved
        region, all the model's tracts, which it scanned then.
        """
        measured = data.path("repetitive_tracts", species, sample, contig=contig, model=model)
        return measured if os.path.exists(measured) else \
            data.path("tandem_repeats", species, sample, contig=contig, model=model)

    # Workflow 4: the scan of every model (step 1), the haploid check (4),
    # final filtering (6).
    if os.path.isdir(f"{cluster}/repetitive_indels"):
        w4 = f"{local}/repetitive_indels"
        model_args = []
        for model in models:
            os.makedirs(f"{w4}/{model}", exist_ok=True)
            model_args += ["--model", model, tracts(model),
                           f"{w4}/{model}/calls.tsv", f"{w4}/{model}/opportunities.tsv"]
        run("tandem_repeat_scan.py", *inputs, f"{w4}/reads.fastq", *scan_options,
            "--reach", max(current["repeats"]["periods"]), "--min-tract-flank", current["scan"]["min_tract_flank"],
            *model_args)
        for model in models:
            replay("repetitive", f"{w4}/{model}/calls.tsv", f"{w4}/reads.fastq", f"{w4}/{model}")
            tables.append((f"workflow 4 {model}", f"{w4}/{model}", f"{cluster}/repetitive_indels/{model}"))

    # Opportunities only for a whole scan job: the cluster counts them per job.
    same = True
    for label, mine, theirs in tables:
        for table, name in (("calls", "calls"), ("opportunities", "opportunities"),
                            ("merged_indels", "calls after the haploid check"), ("indels", "final calls")):
            if table == "opportunities" and not os.path.exists(f"{theirs}/{table}.tsv"):
                continue
            comparison = compare(f"{mine}/{table}.tsv", f"{theirs}/{table}.tsv")
            print(report(f"{label} {name}", comparison), flush=True)
            same = same and comparison.same
    sys.exit(0 if same else 1)


if __name__ == "__main__":
    main()
