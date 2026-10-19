# Owl: how false positives are filtered

Notes on Owl v0.4.0 (Kronenberg et al., *PLoS Comput Biol* 2026, e1014423), from the code in `owl/` and the supplementary methods. Code links point to [owl.rs](owl/src/bin/owl.rs) and [lib.rs](owl/src/lib.rs).

## The pipeline in brief

1. **Alignment and phasing (before Owl).** HiFi reads are aligned to GRCh38 (no-alt) with pbmm2. DeepVariant calls SNVs, and HiPhase phases them and writes an HP (haplotype) and PS (phase block) tag on each read. No de novo assembly is involved.
2. **`owl profile`.** For each of 146,562 microsatellite markers, every read that spans the marker is taken. The read's own sequence over the marker ±5 bp is cut out, and a wrap-around dynamic programming (DP) alignment measures the repeat tract length in that read. The DP is a local alignment against the motif repeated any number of times.
3. **Per haplotype.** CV = 100 × SD / mean of the read lengths.
4. **`owl score`.** A haplotype is *high* if CV ≥ 5. The MSI score is the % of called haplotypes that are high. The paper calls a sample MSI at ≥ 10%.

## Why the design itself suppresses false positives

In normal tissue, all reads from one haplotype at a locus come from a single allele, so their lengths should agree up to sequencing error. Owl therefore measures instability as **spread within a haplotype**, not as difference from the reference or from a matched normal. Several things follow:

- **Germline heterozygosity is removed by splitting reads by HP.** Two different allele lengths are the biggest confounder: about 28% length variation when haplotypes are pooled, against about 1% when they are phased.
- **The threshold sits above HiFi noise.** It was set on a normal genome (NA12878): a gamma fit (shape 3.4, rate 1.41) gives P(CV > 5) ≈ 0.045.
- **No matched normal is needed at runtime.** Normals were used once, to calibrate the threshold and to clean the marker set.
- **Isolated artifacts barely move the score.** Each haplotype is reduced to high or low, and the score averages over about 146k loci. Normals score 1.4–5.4% and MSI tumors 10–27%.

## False-positive sources and the filters against them

### 1. Sequencing error
*How it fakes instability:* an erroneous read gets the wrong tract length, which adds spread.

| Filter | Default | Code |
|---|---|---|
| Mean base quality over the extracted window | ≥ 25 (`--min-bq`) | [owl.rs:534](owl/src/bin/owl.rs#L534) |
| Identity of the repeat alignment | ≥ 90% (`--min-idt`) | [owl.rs:536](owl/src/bin/owl.rs#L536) |
| DP scoring (+2 / −2 / −3): a positive alignment needs > 50% matches | implicit | [lib.rs:195-231](owl/src/lib.rs#L195-L231) |
| CV threshold calibrated on normal-genome noise | CV ≥ 5 (`--cov`) | [owl.rs:876](owl/src/bin/owl.rs#L876) |

### 2. Alignment and mapping error
*How it fakes instability:* reads from the wrong locus, or truncated or misplaced reads, give lengths that don't belong to this allele.

| Filter | Default | Code |
|---|---|---|
| Read must span the marker ±5 bp | hard-coded | [owl.rs:314-327](owl/src/bin/owl.rs#L314-L327) |
| Mapping quality | ≥ 10 (`--min-mapq`) | [owl.rs:347](owl/src/bin/owl.rs#L347) |
| Gap-compressed read identity (pbmm2 `mg` tag) | 0.98 (`--min-mg`), **see caveats** | [owl.rs:332-345](owl/src/bin/owl.rs#L332-L345) |
| No soft clip > 50 bp at either read end | hard-coded | [owl.rs:186-197](owl/src/bin/owl.rs#L186-L197) |
| SAM flag mask | 0, i.e. off (`--flag`) | [owl.rs:347](owl/src/bin/owl.rs#L347) |
| Whole locus no-called if > 20% of spanning reads fail the filters above | 0.2 (`--max-filt-frac`) | [owl.rs:492-500](owl/src/bin/owl.rs#L492-L500) |
| Tract length read from the sequence by the DP, not from the CIGAR | — | [owl.rs:517-519](owl/src/bin/owl.rs#L517-L519) |
| CV instead of SD, which damps outliers that grow with repeat length (e.g. insertions placed outside the window) | — | [owl.rs:571-573](owl/src/bin/owl.rs#L571-L573) |

### 3. Repeat annotation error
*How it fakes instability:* the DP measures the wrong stretch of the read.

- **Near-pure markers only.** Markers must have ≥ 97% identity to a pure repeat (hg38 RepeatMasker), because the DP isn't built for interrupted repeats.
- **Minimum length and identity.** The measured tract must be at least one motif copy long ([owl.rs:535](owl/src/bin/owl.rs#L535)) and have ≥ 90% identity.
- **Robust measurement.** The DP doesn't care which rotation of the motif the tract starts with, and it tolerates small impurities:

```
read : TT[AGCAGCAGCGCAGCAG]TT          motif CAG
       AG|CAG|CAG|C-G|CAG|CAG          → 16 bp, 94% identity
```

### 4. Phasing error
*How it fakes instability:* reads from both haplotypes end up in one group, giving a bimodal length distribution.

| Filter | Default | Code |
|---|---|---|
| Read PS must be 0 or the most common PS at the locus | hard-coded | [owl.rs:506-512](owl/src/bin/owl.rs#L506-L512), [owl.rs:537](owl/src/bin/owl.rs#L537) |
| Unphased reads (HP=0) dropped at sites that also have phased reads | hard-coded | [owl.rs:867-871](owl/src/bin/owl.rs#L867-L871) |
| Drop all HP=0 reads | off (`--unphased`) | [owl.rs:868](owl/src/bin/owl.rs#L868) |
| Warning if no phased region was found | — | [owl.rs:637-639](owl/src/bin/owl.rs#L637-L639) |
| % phased sites reported per sample | — | score output |

### 5. Germline variation and noisy loci
*How it fakes instability:* inherited length differences, or loci that are noisy in everyone.

- **Reads are grouped by HP** before the spread is computed ([owl.rs:543](owl/src/bin/owl.rs#L543)).
- **Marker design keeps clean loci.** Markers are 20–100 bp with motifs ≤ 6 bp, on chr1–22 and chrX only, excluding alt, random and Un contigs ([MARKERS.md](owl/analysis/MARKERS.md)).
- **Population filter.** Markers that were no-called or had CV > 5 in > 25% of 131 HPRC normals were removed: 17,790 markers, leaving the 146,562 shipped markers.
- **Motif-level check (analysis scripts only, not part of Owl).** A Beta-Binomial baseline per motif, built from the 131 controls, flags motifs that are unstable beyond their normal level ([motif-analysis.R](owl-analysis/motif-analysis/motif-analysis.R)).

### 6. Too little data

| Filter | Default | Code |
|---|---|---|
| Reads per haplotype | ≥ 5 (`--min-cov`, `--min-depth`) | [owl.rs:587](owl/src/bin/owl.rs#L587), [owl.rs:867](owl/src/bin/owl.rs#L867) |
| Sample QC fails if < 50% of sites pass | 50 (`--min-pass`) | [owl.rs:1025-1029](owl/src/bin/owl.rs#L1025-L1029) |
| Paper: score is stable at ≥ 20x coverage | — | — |

## Gaps and caveats

- **`--min-mg` is effectively off.** pbmm2 writes `mg` on a 0–100 scale, but Owl's default is `0.98`, so the filter never removes a pbmm2-aligned read. Upstream `main` has the same code. Use `--min-mg 98`, and confirm the scale with `samtools view x.bam | head -1 | tr '\t' '\n' | grep '^mg'`.
- **Reads without an `mg` tag are kept,** for example in BAMs aligned with minimap2.
- **No flag filtering by default,** and only soft clips are checked, not hard clips.
- **The 20% locus rule only counts alignment-level failures.** Reads dropped for base quality, identity or PS don't count toward it.
- **Base quality is averaged over the whole window,** flanks included, not over the tract alone.
- **Insertions placed at or beyond the ±5 bp window edge are cut off.** CV damps this but doesn't remove it, and there is no outlier removal.
- **Only one tract is measured per read** (the best-scoring one).
- **The score counts haplotypes, not loci:** a phased site gives two votes. The README says "loci".
- **The calibration genome is Revio data.** NA12878 is Revio (movie `m84046_231004_192713_s3`), while HPRC Release 1 controls are Sequel II. The marker filter may therefore rest on mixed-platform data.
- **Several `profile` options share the `-m` short flag.** Use the long option names.
- **Co-high output is undocumented.** It reports, per phase block, how often both haplotypes are high compared with the rate expected by chance. It is in the code ([owl.rs:912-996](owl/src/bin/owl.rs#L912-L996)) but not described in the paper. *Interpretation:* an excess likely points to phasing or locus artifacts rather than real instability.
