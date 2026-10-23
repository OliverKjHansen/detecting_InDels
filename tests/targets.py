"""
Run the workflows' own targets on a small sample, with bash, as gwf would: a
checkout under a test's directory, with config.yaml's layout and parameters,
the sample's files written where the workflow before leaves them, and
stand-ins for external tools. Shared by the tests that run targets.
"""

import os
import runpy
import subprocess
import sys

import yaml

from workflow_config import REPO_ROOT

SPECIES, SAMPLE, CONTIG = "sp", "s", "c1"


class Sample:
    """
    One sample with one contig in a checkout under `root`: the workflows'
    files at the layout's paths, its targets runnable by name.
    """

    def __init__(self, root, monkeypatch, sequence):
        self.root, self.monkeypatch = root, monkeypatch
        os.symlink(os.path.join(REPO_ROOT, "scripts"), f"{root}/scripts")
        os.makedirs(f"{root}/bin")
        self.stand_in("samtools", "#!/bin/sh\necho '>window'\n")

        with open(os.path.join(REPO_ROOT, "config.yaml")) as handle:
            config = yaml.safe_load(handle)
        config.update({
            "samples": {SPECIES: [SAMPLE]},
            "tissue_samples": {}, "methylation_samples": {}, "sperm_samples": {}, "sample_age": {},
            "alignment_overrides": {},
            "paths": {"fasta": f"{root}/fasta", "bam": f"{root}/bam", "read_labels": f"{root}/labels",
                      "classification": f"{root}/classification",
                      "output": "output", "results": "results", "plots": "plots"},
            "tools": {"mreps": f"{root}/bin/mreps", "trf": f"{root}/bin/trf"},
            "min_contig_length": 100,
        })
        with open(f"{root}/config.yaml", "w") as handle:
            yaml.safe_dump(config, handle, sort_keys=False)

        monkeypatch.chdir(root)
        monkeypatch.syspath_prepend(os.path.join(REPO_ROOT, "workflows"))
        from pipeline import Pipeline
        self.pipeline = Pipeline()
        for name in ("assembly", "contig_fasta"):
            self.write(name, [f">{CONTIG}", sequence], contig=CONTIG)
        self.write("assembly_fai", [f"{CONTIG}\t{len(sequence)}\t4\t{len(sequence)}\t{len(sequence) + 1}"])

    def path(self, name, **fields):
        return self.pipeline.path(name, SPECIES, SAMPLE, **fields)

    def write(self, name, lines, **fields):
        """Write the lines to the layout's file `name`."""
        path = self.path(name, **fields)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as handle:
            handle.writelines(f"{line}\n" for line in lines)
        return path

    def write_tracts(self, name, tracts, model="M1", **fields):
        """A repeat BED as workflow 2 writes it: header, then (start, end, motif, k) per tract."""
        self.write(name, ["contig\tstart\tend\tcanonical_motif\tk\tcopies\ttrf_purity\tmodel_purity\tmodel"
                          "\tspecies\tsample"]
                   + [f"{CONTIG}\t{start}\t{end}\t{motif}\t{k}\t{(end - start) / k:g}\t100.00\t100.00\t{model}"
                      f"\t{SPECIES}\t{SAMPLE}" for start, end, motif, k in tracts], model=model, **fields)

    def configure(self, block, **values):
        """Change values of one block of the sample's config.yaml, as a test needs them; the targets read it when run."""
        with open(f"{self.root}/config.yaml") as handle:
            config = yaml.safe_load(handle)
        config[block].update(values)
        with open(f"{self.root}/config.yaml", "w") as handle:
            yaml.safe_dump(config, handle, sort_keys=False)

    def stand_in(self, tool, script):
        path = f"{self.root}/bin/{tool}"
        with open(path, "w") as handle:
            handle.write(script)
        os.chmod(path, 0o755)

    def run(self, workflow, target):
        """Run one target of workflows/{workflow}.py with bash, as gwf would."""
        targets = runpy.run_path(os.path.join(REPO_ROOT, "workflows", f"{workflow}.py"))["gwf"].targets
        env = dict(os.environ, PATH=os.pathsep.join([f"{self.root}/bin", os.path.dirname(sys.executable),
                                                     os.environ["PATH"]]))
        done = subprocess.run(["bash", "-c", targets[target].spec], cwd=self.root, env=env,
                              capture_output=True, text=True)
        assert done.returncode == 0, f"{target} failed:\n{done.stderr}"

    def regions(self, name, **fields):
        """The (start, end) intervals of a BED written by a target, header skipped."""
        with open(self.path(name, **fields)) as handle:
            rows = [line.rstrip("\n").split("\t") for line in handle]
        return [(int(row[1]), int(row[2])) for row in rows if row[1].isdigit()]

    def table(self, name, **fields):
        """A summary table written by a target, one dict per row."""
        with open(self.path(name, **fields)) as handle:
            header, *rows = [line.rstrip("\n").split("\t") for line in handle]
        return [dict(zip(header, row)) for row in rows]


def sample(tmp_path, monkeypatch, sequence):
    return Sample(os.path.realpath(tmp_path), monkeypatch, sequence)
