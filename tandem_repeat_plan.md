# Plan: Period-Stratified Tandem Repeat Framework + Repetitive BED

**TL;DR**: `repetitiveDNA.py` runs two parallel repeat-detection pipelines per contig for each k∈{1,...,6}: **mreps** (exact/pure repeats) producing **M0**, and **TRF** (standard impurity-tolerant) producing **M1**. Both models emit one merged BED per model per sample, feeding both InDel pipelines.

---

## Phase 1 — Repeat Detection ✅

### Workflow structure

**One GWF target per contig** (`tandem_repeats_{contig}_{sample}_diploid`) handles all k=1..6 for both M0 and M1. The k-loop runs inside bash; per-k intermediates are written to `WORK=$(mktemp -d)`. Only two files land on disk per contig:

- `output/tandem_repeats/{species}/{sample}/{contig}_M0.bed`
- `output/tandem_repeats/{species}/{sample}/{contig}_M1.bed`

**Target count**: ~1 faSplit + ~# contig targets + 2 cat targets per sample.

---

### M0 pipeline: mreps (exact/pure repeats) ✅

#### Step 1a: Run mreps per k (inside bash loop, TMPDIR)
- Command: `mreps -minperiod $k -maxperiod $k -exp 2 -allowsmall -fasta {contig}.fa > $WORK/${k}_mreps.txt`
- `-exp 2` → minimum exponent ≥ 2.0 (copies = length/period); minimum tract = 2k bp
- `-allowsmall` → allow tracts shorter than random-occurrence threshold

#### Step 1b: Parse + canonicalise mreps output
Script: `endjoining/scripts/parse_mreps_canonical.py`

- Reuses `kolpakov.py` parsing logic (skip first 15 and last 3 lines, field extraction, 0-based coordinate shift)
- `canonicalize(motif)` → lexicographic minimum over all rotations only (no reverse complement — `A` and `T` homopolymers remain distinct)
- `phase = start % k`
- Output TSV columns: `contig, start, end, canonical_motif, k, copies, phase, raw_motif, species, sample`
- Written to `$WORK/${k}_trs_raw.tsv` (TMPDIR — not kept on disk)

#### Step 1c: Produce M0 annotation
Script: `endjoining/scripts/impurity_models.py` with `models = [0]`

- M0 = exact baseline: each mreps tract reported as-is, no merging
- Both purity fields set to 100.0 (exact by definition)
- Output: `$WORK/${k}_M0.bed` (TMPDIR; all k then concatenated to `{contig}_M0.bed`)

---

### M1 pipeline: TRF (standard impurity-tolerant repeats) ✅

#### Step 2a: Run TRF per k (inside bash loop, TMPDIR)
- Command: `trf {contig}.fa 2 7 7 80 10 50 $k -ngs -h | awk -v k=$k '/^@/{print} !/^@/{if ($3==k) print}' > $WORK/${k}_trf.txt`
- Standard parameters: match=2 mismatch=7 indel=7 PM=80 PI=10 minscore=50 maxperiod=k
- `-ngs` → machine-readable stdout, `-h` → suppress HTML; `awk -v k=$k` post-filters to period=k rows only

#### Step 2b: Parse TRF output ✅
Script: `endjoining/scripts/parse_trf.py`

- Converts TRF -ngs rows to the unified BED format
- `trf_purity` = TRF's own SW-alignment % matches (col 6, 1-indexed)
- `model_purity` = Hamming matches/tract_length × 100 (tiled motif vs. actual repeat sequence)
- Canonicalisation: rotation-only `min(motif[i:]+motif[:i] for i in range(k))`
- TRF 1-based inclusive → 0-based half-open (start = start_1based − 1)
- Output: `$WORK/${k}_M1.bed` (TMPDIR; all k then concatenated to `{contig}_M1.bed`)

---

## Preserved future model — M_fission (code in `impurity_models.py`, not yet active)

The fission/fusion logic is preserved in `impurity_models.py` under `model == "M_fission"`. It is **not** wired into the GWF workflow yet. When activated it would merge adjacent same-motif mreps tracts when d = 1, under either:

- **Condition 1 — point mutation (substitution):** the single gap base is a mismatch to the expected motif continuation.
- **Condition 2 — 1-bp indel (phase-shifting):** the gap base matches the motif but the following tract is phase-shifted.

For d = 0 (abutting tracts): always merge.

---

## Unified BED format ✅

All output BED files (M0 and M1) share the same 11-column tab-separated format with a header line:

| Column | Type | Description |
|--------|------|-------------|
| `contig` | str | Contig/chromosome name |
| `start` | int | 0-based start (half-open) |
| `end` | int | 0-based end (half-open) |
| `canonical_motif` | str | Lexicographic-minimum rotation of motif |
| `k` | int | Period length |
| `copies` | float | Number of copies (length / period) |
| `trf_purity` | float | TRF SW-alignment % matches (100.0 for M0) |
| `model_purity` | float | Hamming matches / tract_length × 100 (100.0 for M0) |
| `model` | str | `M0` or `M1` |
| `species` | str | Species name |
| `sample` | str | Sample name |

---

## Model summary

| Model | Source | Merging | trf_purity | model_purity | Notes |
|-------|--------|---------|------------|--------------|-------|
| M0 | mreps | never | 100.0 | 100.0 | exact/pure repeat baseline |
| M1 | TRF | built-in | TRF col 6 | Hamming | standard impurity-tolerant (SW scoring) |
| M_fission | mreps | d=0 or d=1 | — | — | 1-bp fission/fusion; **preserved, not yet active** |

**2 BED files per sample** (M0, M1); M_fission added when activated.

---

## Phase 2 — Concatenate per-contig BEDs into per-sample BEDs ✅

Two GWF `cat` targets per sample, using `awk` to keep the header from the first file only:

- `cat_M0_{sample}_diploid` → `output/tandem_repeats/{species}/{sample}/{sample}_M0.bed`
- `cat_M1_{sample}_diploid` → `output/tandem_repeats/{species}/{sample}/{sample}_M1.bed`

Inputs: all `{contig}_M0.bed` / `{contig}_M1.bed` for large contigs (>1 Mb) of that sample.

---

## Phase 3 — DRL Computation (per sample × k × model) ❌ not yet started

### Step 3: Compute DRL
New script: `endjoining/scripts/drl_compute.py`

- Read `{sample}_M0.bed` and `{sample}_M1.bed` (already contain all k)
- Group by `(canonical_motif, k, model, copies)` and count tracts
- Output: `drl/{sample}_k{k}_M{m}_drl.tsv`
  - Columns: `canonical_motif, k, model, copies, count`
- GWF target: `drl_{sample}_k{k}_M{m}`, inputs = `{sample}_M0.bed` or `{sample}_M1.bed`

---

## Phase 4 — Final Merged BED (per sample × model) ❌ not yet started

### Step 4: Build repetitive BED
Shell step inside GWF — one target per model:

- Take `{sample}_M0.bed` / `{sample}_M1.bed`
- Extract columns 1–3 → `sort -k1,1 -k2,2n | bedtools merge -i -`
- Output: `output/repetitive_bed/{sample}_M{m}_repetitive.bed`
- **2 BED files per sample** (M0, M1); M_fission added when activated

---

## Phase 5 — Plots ❌ not yet started

### Step 5: R plotting
New script: `endjoining/scripts/plot_drl.R`

- Input: all DRL TSVs for a sample (all k, all models)
- Plots:
  - DRL distribution per k, overlaying models (curves per model colour)
  - Tail behaviour (log-scale count vs. copies)
  - Maximum copies and number of long tracts per model
- Output: `output/plots/{sample}_drl_plots.pdf`

---

## Files

### Implemented ✅

| File | Purpose |
|------|---------|
| `endjoining/repetitiveDNA.py` | GWF workflow — faSplit, one target per contig (all k, both models), cat targets |
| `endjoining/scripts/parse_mreps_canonical.py` | mreps parsing + canonical motif + phase |
| `endjoining/scripts/impurity_models.py` | M0 annotation (exact baseline); M_fission preserved but inactive |
| `endjoining/scripts/parse_trf.py` | TRF -ngs parsing → M1 BED (trf_purity + Hamming model_purity) |

### Pending ❌

| File | Purpose |
|------|---------|
| `endjoining/scripts/drl_compute.py` | DRL aggregation (Phase 3) |
| `endjoining/scripts/plot_drl.R` | DRL plots (Phase 5) |

---

## Relevant Existing Code to Reuse

- `microsatellites/scripts/kolpakov.py` — mreps parsing (lines 15:-3 skip, field extraction, 0-based conversion)
- `microsatellites/ratemicrosat/workflow.py` — `faSplit` target pattern, contig dict loading from `.fai`
- `endjoining/workflow.py` — BED-building shell pattern (`cat ... | sort ... | bedtools merge`)
- `endjoining/scripts/utilities.py:load_filters_with_tree` — IntervalTree BED loading (for later InDel integration)

---

## Decisions

- One BED per model (M0, M1), not a single canonical BED
- `phase = start % k` using 0-based coordinates (intermediate only — not in final BED)
- Reverse complement excluded from canonicalisation — `A` and `T` homopolymers remain distinct
- k=1 minimum tract length ≥ 2 bp (met automatically by `-exp 2`)
- Per-k intermediates written to `mktemp -d` scratch dir; `trap 'rm -rf "$WORK"' EXIT` ensures cleanup even on job failure
- InDel calling integration is **out of scope** for this milestone

---

## Verification Checklist

1. Run `parse_mreps_canonical.py` on `development/h1tg000001l.fa`; verify k=2 on a known ATATAT run gives `canonical_motif="AT"`, `phase=0`, `copies≥4`
2. Run `impurity_models.py` on a toy case; verify M0 `trf_purity=100.0`, `model_purity=100.0`
3. Run `parse_trf.py` on a known TRF output; verify `model_purity` matches manual Hamming calculation
4. BED sanity: `awk 'NR>1 && $2>=0 && $3>$2' {bed} | wc -l` matches non-header line count
5. DRL plausibility: for k=2 in a primate genome, peak copy count expected in 3–10 range
6. Check BED covers expected high-repeat regions: `bedtools intersect -a known_repeats.bed -b output_bed -u | wc -l`

---

## Open Questions / TODOs

- **[ ] Find convincing way to determine lower bound on TRF tracts.**  
  minscore=50 with match weight 2 implies a hard minimum of ~25 bp for a perfect tract, but the right biological/statistical floor may differ. Candidate approaches: (1) observed vs. expected tract count under a null model as a function of (k, copies), finding where observed ≫ expected; (2) HiFi copy-number resolution — minimum tract length where per-read copy-count σ < 0.5; (3) slippage literature lower bound (~3–4 copies). Revisit after running M0/M1 on real data.

- **[ ] Verify that excluding reverse complement from canonical motif is correct.**  
  Currently `canonicalize()` uses rotation-only (`min(motif[i:]+motif[:i] for i in range(k))`), so e.g. `AC` and `GT` (reverse complement of `AC`) are treated as distinct motifs. This was chosen to keep `A` and `T` homopolymers distinct. But for motifs with k≥2 it is not obvious this is right: a polymerase stalling on an `AC`-repeat on one strand sees a `GT`-repeat on the other — these are the same locus. Check whether downstream InDel calling or DRL analysis requires strand-collapsed (RC-inclusive) canonical forms, and update `canonicalize()` in both `parse_mreps_canonical.py` and `parse_trf.py` if so.

---

## Out of Scope (this milestone)

- InDel calling integration (BEDs feed both pipelines in the next milestone)
- Multi-species comparative analysis
- Large indel analysis
