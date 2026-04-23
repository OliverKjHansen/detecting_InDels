from os.path import expanduser
from gwf import *

sample_dict = {
            'marmoset':['lambiek',#'marmoset'
            ]
                }

periods = [1, 2, 3, 4, 5, 6] #maybe move to config file later

assembly_type = "diploid"
inpath_fa = expanduser("~/mutationalscanning/DerivedData/fasta/HiFi")
scripts = expanduser("~/mutationalscanning/Workspaces/oliver/endjoining/scripts")
mreps = expanduser("~/mutationalscanning/Workspaces/oliver/microsatellites/development/mreps/mreps") #change placement of this later


gwf = Workflow(
    defaults={
        'cores': 1,
        'memory': '10g',
        'walltime': '4:00:00',
        'account': 'mutationalscanning'})


for species, samples in sample_dict.items():
    for sample in samples:
        index_assembly = f"{inpath_fa}/{species}/{sample}/diploid_assembly/{sample}_diploid.fa.fai"
        assembly = f"{inpath_fa}/{species}/{sample}/diploid_assembly/{sample}_diploid.fa"
        contig_dict = {}

        with open(index_assembly) as contigs:
            for line in contigs:
                contig_dict[line.split()[0]] = line.split()[1]

        # Step 1a: split assembly into per-contig FASTAs (shared input for all k targets)
        gwf.target(
            f"faSplit_{sample}_diploid",
            inputs=[assembly],
            outputs=[f"fasplit/{species}/{sample}/{contig}.fa"
                    for contig, length in contig_dict.items() if int(length) > 10**6]
        ) << f"""
        mkdir -p fasplit/{species}/{sample}
        faSplit byname {assembly} fasplit/{species}/{sample}/
        """

        for contig, length in contig_dict.items():
            if int(length) > 10**6: # we only want to run on contigs >1Mb.)
                contig_fa = f"fasplit/{species}/{sample}/{contig}.fa"

                for k in periods:
                    # Step 1a: run mreps in period-constrained exact mode
                    # -exp 2      → minimum exponent ≥ 2 (exponent = length/period = copies)
                    #               i.e. at least 2 full copies of the motif; min tract = 2k bp
                    #               Note: -exp must be > 1.0 per mreps source (ARG_FLOAT check)
                    # -allowsmall → allow tracts shorter than the default random-occurrence threshold
                    # -minperiod k -maxperiod k → single-period detection per run
                    #add a remove step to not have to many file laying around
                    # might be smarter to loop over the k in the bash script to avoid having submit to many targets in the workflow.
                    gwf.target(
                        f"mreps_k{k}_{contig}_{sample}_diploid",
                        inputs=[contig_fa],
                        outputs=[f"output/tandem_repeats/{species}/{sample}/{contig}_k{k}_mreps.txt"]
                        ) << f"""
                        mkdir -p output/tandem_repeats/{species}/{sample}
                        {mreps} -minperiod {k} -maxperiod {k} -exp 2 -allowsmall -fasta {contig_fa} \
                        > output/tandem_repeats/{species}/{sample}/{contig}_k{k}_mreps.txt
                        """

                    # Step 1b: parse mreps output, canonicalise motifs, add phase
                    # Output TSV: contig, start, end, canonical_motif, k, copies, phase, raw_motif, species, sample
                    gwf.target(
                        f"parse_k{k}_{contig}_{sample}_diploid",
                        inputs=[f"output/tandem_repeats/{species}/{sample}/{contig}_k{k}_mreps.txt"],
                        outputs=[f"output/tandem_repeats/{species}/{sample}/{contig}_k{k}_trs_raw.tsv"]
                        ) << f"""
                        python {scripts}/parse_mreps_canonical.py \
                        output/tandem_repeats/{species}/{sample}/{contig}_k{k}_mreps.txt \
                        {contig} {k} {species} {sample} \
                        > output/tandem_repeats/{species}/{sample}/{contig}_k{k}_trs_raw.tsv
                        """

                    # Step 2: apply impurity models M0–M3 to produce merged tract annotations
                    # M0 = exact baseline (no merging), M1–M3 progressively relax gap/mismatch thresholds
                    # Output: one TSV per model containing merged tract coordinates and copy counts
                    gwf.target(
                        f"impurity_k{k}_{contig}_{sample}_diploid",
                        inputs=[f"output/tandem_repeats/{species}/{sample}/{contig}_k{k}_trs_raw.tsv",
                                assembly],
                        outputs=[f"output/tandem_repeats/{species}/{sample}/{contig}_k{k}_M{m}_merged.tsv"
                                 for m in range(4)]
                        ) << f"""
                        python {scripts}/impurity_models.py \
                        output/tandem_repeats/{species}/{sample}/{contig}_k{k}_trs_raw.tsv \
                        {contig} {k} {assembly} \
                        output/tandem_repeats/{species}/{sample}/{contig}_k{k}
                        """