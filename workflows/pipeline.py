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

    def key_parameters(self):
        """The scans' --key-parameters: the `strata` block, as JSON quoted for the shell."""
        return shell([json.dumps(self.config["strata"], separators=(",", ":"))])

    @property
    def haplotypes(self):
        """{name: prefix}: the diploid assembly's haplotypes, hap1 and hap2, and how their contigs' names start."""
        return self.config["precompute"]["haplotype_prefixes"]

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
    and the `final_filters` block of the config, with workflow 4's own filters
    when `repeats`.
    """
    filters = config["final_filters"]
    sequel_ii = sample in ((config.get("sequel_ii_samples") or {}).get(species) or [])
    options = ["--platform", "Sequel II" if sequel_ii else "Revio",
               "--capped-quality", str(filters["capped_quality"])]
    for column, values in {**filters["keep"], **(filters["repeat_keep"] if repeats else {})}.items():
        options += ["--keep", f"{column}={','.join(map(str, values))}"]
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


def compressor(output):
    """The command that writes a pipe to `output`: gzip if its name ends in .gz."""
    return "gzip -6" if str(output).endswith(".gz") else "cat"


def table_bodies(inputs):
    """A shell command that writes the rows of tables, plain or gzipped, without their headers."""
    return f"for table in {' '.join(str(path) for path in inputs)}; do gzip -dcf $table | tail -n +2; done"


def joined_tables(inputs):
    """A shell command that writes tables, plain or gzipped, as one: the first one's header, then every row."""
    return f"{{ gzip -dcf {inputs[0]} | head -n 1; {table_bodies(inputs)}; }}"


def joined_reads(inputs, output):
    """A shell command that joins FASTQs, plain or gzipped, keeping each read once, by its name."""
    return (f"gzip -dcf {' '.join(str(path) for path in inputs)} "
            f"| awk 'NR % 4 == 1 {{keep = !seen[$1]++}} keep' | {compressor(output)} > {output}")


def summed_counts(inputs, output):
    """
    A shell command that joins count tables from the windows of a contig: the
    header once, then the last column (count) summed over every distinct value
    of all the others, whatever keys the scan wrote. The tables may be
    gzipped, and the output is if its name ends in .gz. The rows are sorted
    first, so equal keys are adjacent, and summed run by run: memory doesn't
    grow with the number of distinct rows, and sort spills to temporary files
    next to the output.
    """
    folder = os.path.dirname(str(output)) or "."
    return (f"{{ gzip -dcf {inputs[0]} | head -n 1; {table_bodies(inputs)} "
            f"| LC_ALL=C sort -S 2G -T {folder} "
            "| awk -F'\\t' -v OFS='\\t' '{key = $1; for (i = 2; i < NF; i++) key = key OFS $i; "
            "if (NR > 1 && key != last) print last, n; if (NR == 1 || key != last) n = 0; n += $NF; last = key} "
            "END {if (NR) print last, n}'; } "
            f"| {compressor(output)} > {output}")
