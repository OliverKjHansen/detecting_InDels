from os.path import expanduser
from gwf import *

sample_dict = {
            "marmoset":["marmoset","lambiek"
            ]
                }

periods = [1, 2, 3, 4, 5, 6] #maybe move to config file later

assembly_type = "diploid"
inpath_fa = expanduser("~/mutationalscanning/DerivedData/fasta/HiFi")
scripts = expanduser("~/mutationalscanning/Workspaces/oliver/endjoining/scripts")
mreps = expanduser("~/mutationalscanning/Workspaces/oliver/microsatellites/development/mreps/mreps")
trf = expanduser("~/mutationalscanning/Workspaces/oliver/microsatellites/development/TRF/build/src/trf")


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

        # One GWF target per contig handles all k for both M0 (mreps) and M1 (TRF).
        # The k-loop runs inside bash; per-k intermediates live in a TMPDIR and are
        # deleted on exit. Only {contig}_M0.bed and {contig}_M1.bed land on disk.
        periods_str = " ".join(str(k) for k in periods)
        m0_k_files  = " ".join(f"$WORK/{k}_M0.bed" for k in periods)
        m1_k_files  = " ".join(f"$WORK/{k}_M1.bed" for k in periods)
        large_contigs = [c for c, l in contig_dict.items() if int(l) > 10**6]

        for contig in large_contigs:
            contig_fa = f"fasplit/{species}/{sample}/{contig}.fa"
            out_dir   = f"output/tandem_repeats/{species}/{sample}"

            gwf.target(
                f"tandem_repeats_{contig}_{sample}_diploid",
                inputs=[contig_fa, assembly],
                outputs=[f"{out_dir}/{contig}_M0.bed",
                         f"{out_dir}/{contig}_M1.bed"]
            ) << f"""
            mkdir -p {out_dir}
            WORK=$(mktemp -d)
            trap 'rm -rf "$WORK"' EXIT

            for k in {periods_str}; do

                # ---- M0: mreps exact/pure repeats ----
                # -exp 2:      at least 2 full copies (exponent = length/period)
                # -allowsmall: allow tracts shorter than random-occurrence threshold
                {mreps} -minperiod $k -maxperiod $k -exp 2 -allowsmall -fasta {contig_fa} \
                    > $WORK/${{k}}_mreps.txt
                python {scripts}/parse_mreps_canonical.py \
                    $WORK/${{k}}_mreps.txt {contig} $k {species} {sample} \
                    > $WORK/${{k}}_trs_raw.tsv
                python {scripts}/impurity_models.py \
                    $WORK/${{k}}_trs_raw.tsv {contig} $k {assembly} $WORK/$k

                # ---- M1: TRF standard impurity-tolerant detection ----
                # Parameters: match=2 mismatch=7 indel=7 PM=80 PI=10 minscore=50
                # -ngs → machine-readable stdout, -h → suppress HTML
                # awk -v passes bash $k into awk as a numeric variable
                {trf} {contig_fa} 2 7 7 80 10 50 $k -ngs -h \
                    | awk -v k=$k '/^@/{{print}} !/^@/{{if ($3==k) print}}' \
                    > $WORK/${{k}}_trf.txt
                python {scripts}/parse_trf.py \
                    $WORK/${{k}}_trf.txt $k {species} {sample} \
                    > $WORK/${{k}}_M1.bed

            done

            # Concatenate all k — header retained from the first file only
            awk 'FNR==1 && NR==1 {{print}} FNR>1 {{print}}' \
                {m0_k_files} > {out_dir}/{contig}_M0.bed
            awk 'FNR==1 && NR==1 {{print}} FNR>1 {{print}}' \
                {m1_k_files} > {out_dir}/{contig}_M1.bed

            rm -rf $WORK
            """

        # ----------------------------------------------------------------
        # Concatenate all per-contig BED files into one per model per sample
        # ----------------------------------------------------------------
        m0_inputs = [f"output/tandem_repeats/{species}/{sample}/{c}_M0.bed" for c in large_contigs]
        m1_inputs = [f"output/tandem_repeats/{species}/{sample}/{c}_M1.bed" for c in large_contigs]

        gwf.target(
            f"cat_M0_{sample}_diploid",
            inputs=m0_inputs,
            outputs=[f"output/tandem_repeats/{species}/{sample}/{sample}_M0.bed"]
        ) << f"""
        awk 'FNR==1 && NR==1 {{print}} FNR>1 {{print}}' \
        {" ".join(m0_inputs)} \
        > output/tandem_repeats/{species}/{sample}/{sample}_M0.bed
        """

        gwf.target(
            f"cat_M1_{sample}_diploid",
            inputs=m1_inputs,
            outputs=[f"output/tandem_repeats/{species}/{sample}/{sample}_M1.bed"]
        ) << f"""
        awk 'FNR==1 && NR==1 {{print}} FNR>1 {{print}}' \
        {" ".join(m1_inputs)} \
        > output/tandem_repeats/{species}/{sample}/{sample}_M1.bed
        """

