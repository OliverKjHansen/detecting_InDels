"""
A config.yaml for loading the workflows without the cluster.

Every key comes from the repository's config.yaml except the samples, the
paths and the tools, which point into a fixture root holding one assembly
index (.fai) per sample. That is all a workflow reads when it is loaded, so
every target it would run can be listed. The contigs and samples are chosen to
reach every branch:

  ape/s1     tissue labels, methylation scores; c1 spans two 10 Mb scan
             windows and three 5 Mb TRF windows; c3 is exactly
             min_contig_length, which is too short (the test is strict)
  ape/s2     an alignment override, no labels
  monkey/s3  a sperm sample of another species; c3 is far too short
"""

import os

import yaml

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SAMPLES = {"ape": ["s1", "s2"], "monkey": ["s3"]}

CONTIGS = {
    ("ape", "s1"): {"c1": 12_000_000, "c2": 1_500_000, "c3": 1_000_000},
    ("ape", "s2"): {"c2": 1_500_000},
    ("monkey", "s3"): {"c2": 1_500_000, "c3": 500_000},
}


def write_config(root):
    """Write config.yaml and the assembly indexes under `root`; return the config."""
    with open(os.path.join(REPO_ROOT, "config.yaml")) as handle:
        config = yaml.safe_load(handle)

    config.update({
        "samples": SAMPLES,
        "tissue_samples": {"ape": ["s1"]},
        "methylation_samples": {"ape": ["s1"]},
        "sperm_samples": {"monkey": ["s3"]},
        "sample_age": {"ape": {"s1": 30, "s2": 12}, "monkey": {"s3": 7}},
        "alignment_overrides": {"ape": {"s2": f"{root}/override/s2_diploid.bam"}},
        "paths": {
            "fasta": f"{root}/fasta",
            "bam": f"{root}/bam",
            "read_labels": f"{root}/read_labels",
            "classification": f"{root}/classification",
            "microsat": f"{root}/microsat",
            # Relative, as in the real config: the workflows make them
            # absolute from the directory gwf runs in.
            "output": "output",
            "results": "results",
            "plots": "plots",
        },
        "tools": {"mreps": f"{root}/tools/mreps", "trf": f"{root}/tools/trf"},
    })

    kind = config["assembly_type"]
    for (species, sample), contigs in CONTIGS.items():
        fai = f"{root}/fasta/{species}/{sample}/{kind}_assembly/{sample}_{kind}.fa.fai"
        os.makedirs(os.path.dirname(fai), exist_ok=True)
        with open(fai, "w") as handle:
            for contig, length in contigs.items():
                handle.write(f"{contig}\t{length}\t0\t60\t61\n")

    with open(os.path.join(root, "config.yaml"), "w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)
    return config
