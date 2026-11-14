"""
What every workflow needs from config.yaml: the samples, their paths built
from the `layout` templates, and the contigs long enough to analyse.

    from pipeline import Pipeline

    pipeline = Pipeline()
    assembly = pipeline.path("assembly", species, sample)
    for contig, length in pipeline.long_contigs(species, sample).items():
        fasta = pipeline.path("contig_fasta", species, sample, contig=contig)

gwf puts the workflow's folder on the import path, so the workflows import
this file directly. Run gwf from the repository root: config.yaml is read from
the current directory, and relative paths in it, such as `output`, are made
absolute from there, because gwf runs the recipes from workflows/.
"""

import json
import os
import shlex
import string

import yaml

# A sub-window the long-period TRF doesn't finish (ticket 36): its window's BED
# gives its core with this in the motif column.
TIMED_OUT = "timed_out"

# The entries of the `strata` block that only workflow 3's scan reads (ticket 49).
NONREPETITIVE_STRATA = ("flank_quality_bins", "min_inserted_quality", "nonrepetitive_without", "nonrepetitive_bins")


class Pipeline:
    def __init__(self, config_path="config.yaml", base=None):
        """
        Relative paths in the config are taken from `base`, by default the
        current directory. `scripts` is always the scripts/ folder of the
        current directory, the repository root gwf runs in.
        """
        with open(config_path) as handle:
            self.config = yaml.safe_load(handle)
        base = base or os.getcwd()
        self.roots = {name: os.path.abspath(os.path.join(base, os.path.expanduser(value)))
                      for name, value in self.config["paths"].items()}
        self.scripts = os.path.abspath("scripts")

    def is_listed(self, table, species, sample):
        """Whether the sample is listed under its species in a table such as `tissue_samples`."""
        return sample in ((self.config.get(table) or {}).get(species) or [])

    def path(self, name, species, sample, **fields):
        """
        The `layout` template `name` for one sample, filled in from `paths` and
        `assembly_type` plus `fields` such as contig or model.
        """
        return self.config["layout"][name].format(
            **self.roots, assembly_type=self.config["assembly_type"],
            species=species, sample=sample, **fields)

    def folder(self, name, species, sample, **fields):
        """
        The folder the files of `layout` template `name` go in, for one sample.
        `fields` fills the placeholders in the folder part, such as model; the
        per-file ones, such as contig, may only appear in the file name.
        """
        filled = {"species", "sample", "assembly_type", *self.roots, *fields}
        blanks = {field: "\0" for _, field, _, _ in string.Formatter().parse(self.config["layout"][name])
                  if field and field not in filled}
        folder = os.path.dirname(self.path(name, species, sample, **fields, **blanks))
        if "\0" in folder:
            raise ValueError(f"layout.{name} has a per-file field in its folder")
        return folder

    def alignment(self, species, sample):
        """The sample's BAM: its entry under `alignment_overrides`, if any, else the layout's."""
        overrides = (self.config.get("alignment_overrides") or {}).get(species) or {}
        if sample in overrides:
            return os.path.expanduser(overrides[sample])
        return self.path("alignment", species, sample)

    def window_size(self, species, sample):
        """The sample's workflow 3 window, bp: its entry under scan.window_size_by_sample, if any, else scan.window_size."""
        scan = self.config["scan"]
        return ((scan.get("window_size_by_sample") or {}).get(species) or {}).get(sample, scan["window_size"])

    def classification_arguments(self):
        """
        The last four arguments of join_read_classification.sh and
        classification_read_counts.sh, from the `classification` block: the
        germline prefix, the somatic cell types, the bin edges, the bin labels.
        """
        block = self.config["classification"]
        values = [block["germline_prefix"], ",".join(block["somatic_types"]),
                  ",".join(map(str, block["bin_edges"])), ",".join(map(str, block["bin_labels"]))]
        return shell(values)

    def key_parameters(self, nonrepetitive=False):
        """
        The scans' --key-parameters: the `strata` block, as JSON quoted for the
        shell. Workflow 3's own entries only for workflow 3 (`nonrepetitive`),
        so workflow 4's commands don't change with them.
        """
        strata = {key: value for key, value in self.config["strata"].items()
                  if nonrepetitive or key not in NONREPETITIVE_STRATA}
        return shell([json.dumps(strata, separators=(",", ":"))])

    @property
    def haplotypes(self):
        """{name: prefix}: the diploid assembly's haplotypes, hap1 and hap2, and how their contigs' names start."""
        return self.config["precompute"]["haplotype_prefixes"]

    def long_period_window_beds(self, species, sample, contig, length):
        """
        The long-period BEDs repetitive_dna.py writes for a contig of `length`,
        one per M1 window, which also give the sub-windows marked TIMED_OUT.
        """
        return [self.path("long_periods_window", species, sample, contig=contig, start=start)
                for start in range(0, length, self.config["repeats"]["trf_window_size"])]

    def long_contigs(self, species, sample):
        """
        {contig: length} for the contigs of the sample's assembly longer than
        `min_contig_length`, in the order of its index.
        """
        contigs = {}
        with open(self.path("assembly_fai", species, sample)) as handle:
            for line in handle:
                contig, length = line.split()[:2]
                if int(length) > self.config["min_contig_length"]:
                    contigs[contig] = int(length)
        return contigs


def final_filtering_options(config, species, sample, repeats=False):
    """
    final_filtering.py's options for one sample, as arguments: its platform
    and the `final_filters` block of the config: the filters of both
    workflows, then workflow 4's own when `repeats`, or else workflow 3's own
    and its calls-only steps.
    """
    filters = config["final_filters"]
    sequel_ii = sample in ((config.get("sequel_ii_samples") or {}).get(species) or [])
    options = ["--platform", "Sequel II" if sequel_ii else "Revio",
               "--capped-quality", str(filters["capped_quality"])]
    own = filters["repeat_keep"] if repeats else filters["nonrepetitive_keep"]
    for column, values in {**(filters["keep"] or {}), **own}.items():
        options += ["--keep", f"{column}={','.join(map(str, values))}"]
    if not repeats:
        for column, values in filters["nonrepetitive_calls_only"].items():
            options += ["--calls-only", f"{column}={','.join(map(str, values))}"]
    return options


def shell(arguments):
    """Arguments as one string for a shell command."""
    return " ".join(shlex.quote(argument) for argument in arguments)


def remapping(reference, reads, bam, locations):
    """
    A shell command that aligns reads (FASTQ) to a reference with pbmm2 and
    writes each alignment's read, flag, contig, position, mapQ and CIGAR,
    sorted by read: the records the haploid check reads.
    """
    return (f"pbmm2 align {reference} {reads} {bam} --preset HIFI --sort -j 12 -J 4 -m 8G && "
            f"samtools view {bam} | awk -v OFS='\\t' '{{print $1,$2,$3,$4,$5,$6}}' | sort -k1,1 > {locations}")


def into_file(output):
    """
    Shell commands that write a pipe into `output`, gzipped if its name ends
    in .gz: to a temporary file first, renamed once the pipe has ended.
    """
    compress = "gzip -6" if str(output).endswith(".gz") else "cat"
    return f"{compress} > {output}.tmp; mv {output}.tmp {output}"


# The joins below run in a target's bash with `set -e`, where a pipeline's
# status is its last command's. They set pipefail, so a window that can't be
# read stops the job instead of leaving a partial table, and they write
# through a temporary file, so a job stopped part-way leaves no file gwf would
# take as done (ticket 37). The opportunity tables are joined by
# scripts/opportunity_tables.py (ticket 40).


def grouped_by_length(lengths, budget):
    """
    The contigs in their order, in groups whose lengths add up to at most
    `budget` bp; a contig longer than the budget is a group of its own.
    `lengths`: {contig: bp}.
    """
    groups, total = [], 0
    for contig, length in lengths.items():
        if groups and total + length <= budget:
            groups[-1].append(contig)
            total += length
        else:
            groups.append([contig])
            total = length
    return groups


def retried(command, attempts=3):
    """
    A shell command that runs `command` up to `attempts` times, until it
    succeeds, and fails if every attempt does: for a step that fails at random
    and leaves no output when it does (ticket 41's opportunity sums).
    """
    return (f"for attempt in $(seq {attempts}); do {command} && break; "
            f"[ \"$attempt\" -lt {attempts} ] || exit 1; done")


def joined_calls(inputs, output):
    """A shell command that joins call tables, with the first one's header once."""
    return (f"set -o pipefail; awk 'FNR==1 && NR==1 || FNR>1' {' '.join(str(path) for path in inputs)} "
            f"| {into_file(output)}")


def joined_reads(inputs, output):
    """A shell command that joins FASTQs, plain or gzipped, keeping each read once, by its name."""
    return (f"set -o pipefail; gzip -dcf {' '.join(str(path) for path in inputs)} "
            f"| awk 'NR % 4 == 1 {{keep = !seen[$1]++}} keep' | {into_file(output)}")
