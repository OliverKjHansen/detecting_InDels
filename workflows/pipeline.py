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
        return " ".join(shlex.quote(value) for value in values)

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


def summed_counts(inputs, output):
    """
    A shell command that joins count tables from the windows of a contig: the
    header once, then the last column (count) summed over every distinct value
    of all the others, whatever keys the scan wrote.
    """
    return ("awk -F'\\t' -v OFS='\\t' 'FNR==1 {if (NR==1) print; next} "
            "{key = $1; for (i = 2; i < NF; i++) key = key OFS $i; n[key] += $NF} "
            "END {for (key in n) print key, n[key]}' "
            f"{' '.join(str(path) for path in inputs)} > {output}")
