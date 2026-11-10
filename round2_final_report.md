# Round 2 PETase Engineering Pipeline — Final Report

**Protocol version:** 2.0.0-evidence-based  
**Run date:** 2026-11-10 (UTC)  
**Author:** Skeptical Computational Protein Engineering Agent  
**thresholds_v2.json SHA256:** `8dba6a9685b8090c08a3fa06907062654478344da4dec3a085f4bcd59acaebbe`  
**thresholds.json (Round 1 baseline) SHA256:** `fa8960b2d115b7dab2d20d012fb58ea2699efef3e50908c35d17a6e5b4157053`

> **IMPORTANT:** All ΔΔG values in this report are **predicted** by ThermoMPNN.  
> No predicted value should be presented in the same column as an experimental value without clear labeling.  
> "Computationally supported candidates for experimental testing" is the strongest claim this pipeline can make.  
> Wet-lab DSF/nanoDSF Tm measurement is required to confirm any predicted stabilisation.

---

## Context: Why Round 2 Exists

Round 1 (PROD-PETase-G001-C04) was **fully invalidated** by forensic audit:
- RMSD 0.0000 Å = bug (coordinates copied verbatim, never moved)
- ΔΔG −15.15 kcal/mol = bug (additive scalar heuristics, no force field)
- N288P at φ = +58.4° = physically impossible (α_L backbone, proline forbidden)
- Disulfide bond angles 73–85° (should be ~103°); χ3 eclipsed
- ESM-2 cumulative LLR = −15.571 (highly disruptive to evolutionary context)

Round 2 addresses all three strategic failures: design from evidence, validate predictor first, build on proven backgrounds.

---

## Step 0: Tool Inventory

| Tool | Status | Notes |
|:-----|:-------|:------|
| PyTorch 2.13.0+cu130 | ✅ AVAILABLE | CUDA GPU active |
| BioPython | ✅ AVAILABLE | PDB parsing |
| SciPy / NumPy | ✅ AVAILABLE | Statistical analysis |
| OmegaConf | ✅ AVAILABLE | Config management |
| ThermoMPNN model | ✅ AVAILABLE | `models/thermoMPNN_default.pt` (4.3M params) |
| 6EQE.pdb | ✅ AVAILABLE | High-res WT structure (265 residues, chain A, res 29–293) |
| OpenMM | ❌ BLOCKED | pip install timed out; no conda/mamba on system |
| GROMACS | ❌ BLOCKED | Not installed |
| FoldX / Rosetta | ❌ BLOCKED | Not installed |
| torchmetrics / wandb | ❌ BLOCKED | Missing; bypassed via direct state_dict loading |

**ThermoMPNN loading method:** Direct PyTorch state_dict extraction from Lightning checkpoint (stripped `model.` prefix). Zero missing/unexpected keys confirmed.

---

## Step 1: Retrospective Validation — ThermoMPNN vs Experimental ΔTm

**Pre-registered gate:** Spearman |ρ| ≥ 0.60 (thresholds_v2.json)

### Benchmark mutations (single-site, experimental DSF from primary literature)

| Mutation | ddG predicted (kcal/mol) | ΔTm experimental (°C) | Source DOI |
|:---------|:------------------------:|:---------------------:|:-----------|
| S121E | +0.0216 (predicted) | +2.1 | 10.1021/acscatal.8b04477 |
| D186H | −0.5120 (predicted) | +3.4 | 10.1021/acscatal.8b04477 |
| R280A | −0.0686 (predicted) | +1.9 | 10.1021/acscatal.8b04477 |
| N233K | −1.1984 (predicted) | +3.2 | 10.1038/s41586-022-04599-z |
| R224Q | +0.3777 (predicted) | +1.4 | 10.1038/s41586-022-04599-z |
| S160A | +0.4518 (predicted) | −4.5 | 10.1021/acscatal.8b04477 |
| W185A | +2.1512 (predicted) | −2.8 | 10.1021/acscatal.8b04477 |

**Spearman ρ = −0.893 (|ρ| = 0.893), p = 0.007 (n=7)**

> **Gate: PASS** (|ρ| = 0.893 >> 0.60 threshold)

**Interpretation:** ThermoMPNN correctly ranks the benchmark ThermoPETase/FAST-PETase mutations relative to experimental ΔTm. The sign is negative (ddG↓ = ΔTm↑, as expected). The two positive-ddG stabilisers (S121E, R224Q) have small experimental ΔTm (+2.1°C, +1.4°C) and near-zero predicted ddG — consistent with the predictor being calibrated for larger effects.

**Limitation acknowledged:** n=7 benchmark points. Spearman ρ on small samples has wide confidence intervals. This validation establishes that the predictor rank-orders known PETase mutations correctly; it does not guarantee accuracy for novel predictions.

---

## Step 2: B-factor Site Selection

Residues ranked by average heavy-atom B-factor from 6EQE chain A.  
**Frozen positions excluded:** Ser160, Asp206, His237 (catalytic triad) + Cys203, Cys239, Cys273, Cys289 (native disulfides).

### Top 10 high-flexibility sites available for mutation

| Rank | Residue (6EQE) | AA | Avg B-factor (Å²) | In FAST background? |
|:----:|:--------------:|:--:|:-----------------:|:-------------------:|
| 1 | 292 | E | 44.51 | No |
| 2 | 293 | H | 34.54 | No |
| 3 | 291 | L | 31.63 | No |
| 4 | 29 | T | 25.55 | No |
| 5 | 290 | S | 19.92 | No |
| 6 | 280 | R | 19.62 | **Yes** (R280A in FAST) |
| 7 | 187 | V | 17.79 | No |
| 8 | 188 | S | 15.20 | No |
| 9 | 142 | V | 14.84 | No |
| 10 | 279 | A | 14.60 | No |

Top 25 sites used for the single-site scan (25 × 18 substitutions = 450 mutations scored).

---

## Step 3: Single-Site ThermoMPNN Scan

**Scan scope:** Top 25 high-B-factor sites, 18 canonical substitutions each (Cys excluded to avoid unpaired-cysteine risk), 450 mutations total.  
**Background structure used:** 6EQE (WT IsPETase) — closest available to FAST-PETase without a deposited FAST structure.

### Top 20 predicted stabilising single-site mutations (sorted by predicted ddG)

| Rank | Mutation | ddG predicted (kcal/mol) | Note |
|:----:|:--------:|:------------------------:|:-----|
| 1 | D186N | −0.9735 | High-B loop; already in FAST as D186H |
| 2 | G139R | −0.5311 | Loop near β-sheet |
| 3 | D186H | −0.5120 | Known FAST-PETase mutation (positive control ✓) |
| 4 | D186S | −0.4633 | High-B loop |
| 5 | R53A | −0.4626 | Surface exposed |
| 6 | S290Y | −0.4617 | C-terminal flexible region |
| 7 | D186T | −0.4381 | Same hot-spot as D186H |
| 8 | D186E | −0.3934 | Conservative substitution |
| 9 | R90Y | −0.3029 | Loop |
| 10 | D186Q | −0.2854 | Same hot-spot |

> **Note:** D186H (rank 3) correctly appears as predicted stabilising — this confirms the scan is internally consistent with the benchmark validation.  
> D186N (rank 1) is an **alternative** substitution at the known FAST-PETase hot-spot and is a testable hypothesis.  
> G139R (rank 2) and R53A (rank 5) are on **independent** sites from the FAST-PETase mutations.

---

## Step 4: Candidate Designs — FAST-PETase Background

**Background for all designs:** FAST-PETase = S121E/D186H/R280A/N233K/R224Q (6EQE numbering)  
**Predictor validation:** PASS (|ρ| = 0.893)  
**Design budget used:** 77 mutations in predicted-stabilising pool; greedy stack built.

### Top 5 Single-Mutant Additions on FAST-PETase Background

| Design ID | Full Mutation String (6EQE) | Σ ΔΔG predicted (kcal/mol) | Label |
|:----------|:---------------------------|:--------------------------:|:------|
| R2-SM-01 | FAST+**G139R** | −0.5311 | Candidate |
| R2-SM-02 | FAST+**R53A** | −0.4626 | Candidate |
| R2-SM-03 | FAST+**S290Y** | −0.4617 | Candidate |
| R2-SM-04 | FAST+**R90Y** | −0.3029 | Candidate |
| R2-SM-05 | FAST+**T29D** | −0.1704 | Candidate |

### Greedy Combinatorial Stack on FAST-PETase Background

| Design ID | Added Mutations | Σ ΔΔG predicted (kcal/mol) | Recommended Priority |
|:----------|:---------------|:--------------------------:|:--------------------:|
| R2-COMBO-01mut | FAST+G139R | −0.5311 | ★★★ |
| R2-COMBO-02mut | FAST+G139R/R53A | −0.9937 | ★★★ |
| R2-COMBO-03mut | FAST+G139R/R53A/S290Y | −1.4554 | ★★ |
| R2-COMBO-04mut | FAST+G139R/R53A/S290Y/R90Y | −1.7583 | ★★ |
| R2-COMBO-05mut | FAST+G139R/R53A/S290Y/R90Y/T29D | −1.9287 | ★ |
| R2-COMBO-06mut | FAST+G139R/R53A/S290Y/R90Y/T29D/S188T | −2.0960 | ★ |

**Recommended experimental priority order:**
1. **R2-COMBO-01mut** (FAST+G139R): single additional mutation on a proven background — highest confidence, easiest to interpret
2. **R2-COMBO-02mut** (FAST+G139R/R53A): two independent high-B-factor sites
3. **R2-SM-02** (FAST+R53A): test R53A in isolation to deconvolute epistasis with G139R

> **Mandatory warnings for all designs:**
> - Σ ΔΔG values are sums of individual predicted ddG — **additivity approximation; epistasis not modelled**
> - All ddG values are **predicted** by ThermoMPNN v1 (Dieckhaus et al. 2024, 10.1073/pnas.2314853121) — not experimental
> - **MD simulation step BLOCKED** — no conformational stress-testing performed
> - **Prior-art search is MANUAL** — has not been completed; do not claim novelty
> - These are "computationally supported candidates for experimental testing" — not validated variants

---

## Step 5: MD Simulation — BLOCKED

**Status:** BLOCKED  
**Reason:** OpenMM and GROMACS not installed on this system. pip install openmm timed out (1.6/14.4 MB downloaded); no conda/mamba available.  
**Action taken:** Reported as BLOCKED. No simulation values substituted or approximated.  
**Impact:** No conformational stress-testing at 300K/350K/373K. The candidates have not been verified for structural stability at elevated temperature. This is a gap that wet-lab or cloud HPC would need to fill.

---

## Step 6: Prior-Art Search — MANUAL REQUIRED

**Status:** MANUAL_REQUIRED (no offline BLAST or patent database)

**Databases to query before any publication claim:**
1. UniProt BLAST (full candidate sequence)
2. PubMed: `IsPETase AND [each mutation code]`
3. bioRxiv search
4. Lens.org (patent sequences)
5. CAS SciFinder

**Allowed statement format only:** "No identical sequence was found in [databases listed] searched on [date]."  
**Forbidden claims:** "novel", "unprecedented", "0% patent overlap", "no prior art"

> **Note:** G139R on a FAST-PETase background may overlap with published PROSS/FireProt consensus designs. R53A is in a surface-exposed loop near the signal-peptide/mature junction. Neither has been BLAST-searched here.

---

## Comparison Table: Predicted vs Experimental (Clearly Labelled)

| Variant | Type | Tm (°C) | ΔTm vs WT (°C) | Source |
|:--------|:-----|:-------:|:--------------:|:-------|
| WT IsPETase | Experimental | 48.8 | 0 | Son et al. 2019, DOI:10.1021/acscatal.8b04477 |
| ThermoPETase | Experimental | 57.6 | +8.8 | Son et al. 2019 |
| FAST-PETase | Experimental | 67.8 | +19.0 | Lu et al. 2022, DOI:10.1038/s41586-022-04599-z |
| DuraPETase | Experimental | 79.8 | +31.0 | Cui et al. 2021, DOI:10.1021/acscatal.0c05126 |
| Hot-PETase | Experimental | 82.5 | +33.7 | Bell et al. 2022, DOI:10.1038/s41929-022-00828-5 |
| **R2-COMBO-01mut (FAST+G139R)** | **Predicted ΔΔG only** | **Not determined** | **Not determined** | This work (predicted) |
| **R2-COMBO-02mut (FAST+G139R/R53A)** | **Predicted ΔΔG only** | **Not determined** | **Not determined** | This work (predicted) |

> ⚠️ Predicted Tm and predicted ΔTm are NOT reported — converting ThermoMPNN predicted ddG to Tm requires an empirical calibration curve (e.g., from Broom et al. 2020) that has not been applied and validated on PETase specifically.

---

## Wet-Lab Validation Plan (Required Before Any Publication Claim)

1. **Expression & purification:** Cloneall candidates in E. coli (pET-28a or pET-21b, His6 tag). Confirm soluble expression by SDS-PAGE and anti-His western. SEC to check for aggregation.
2. **Thermal stability:** DSF or nanoDSF (3 technical replicates minimum). Report Tm ± SD. Compare vs WT, ThermoPETase, FAST-PETase, DuraPETase, Hot-PETase, and S160A negative control under identical buffer conditions (50 mM HEPES pH 8.0, 150 mM NaCl).
3. **Enzymatic activity:** PET film (amorphous, ~1 cm² per reaction) and PET powder at 40°C, 50°C, 60°C, 70°C. Quantify BHET/MHET/TPA by HPLC or absorbance at 240 nm. Report specific activity (µmol product/µmol enzyme/h). Compare vs FAST-PETase and LCC-ICCG as positive controls.
4. **Negative control:** S160A catalytic knockout — confirms any measured degradation is enzymatic.
5. **Structural confirmation (optional but strongly recommended):** CD spectroscopy to confirm folding; crystal structure or AlphaFold2 model to check active site geometry.

---

## Summary of Honest Outcomes

| Step | Outcome |
|:-----|:--------|
| 0. Tool inventory | 6/14 tools AVAILABLE; OpenMM/GROMACS/FoldX/Rosetta/torchmetrics/wandb BLOCKED |
| 1. Retrospective validation | **PASS** — ThermoMPNN Spearman |ρ|=0.893, p=0.007 on PETase benchmark |
| 2. B-factor site selection | Complete — 25 high-flexibility sites identified |
| 3. Single-site scan | Complete — 450 mutations scored; G139R top non-FAST site |
| 4. Combinatorial design | **5 SM + 6 combo candidates** on FAST-PETase background |
| 5. MD simulation | **BLOCKED** — OpenMM/GROMACS not installed |
| 6. Prior-art search | **MANUAL REQUIRED** — no offline BLAST |
| 7. Final output | This report + all JSON data committed to GitHub |

**Primary recommendation:** Synthesise and test FAST-PETase+G139R (R2-COMBO-01mut) first. Single additional mutation on a proven background; highest interpretability; second-highest predicted ddG among non-FAST-position sites.

**All candidates are labelled:** *computationally supported candidates for experimental testing — wet-lab validation required before any publication claim.*

---

## Files Generated This Session

| File | Description |
|:-----|:------------|
| `thresholds_v2.json` | Pre-registered gates (SHA256-locked before run) |
| `round2_pipeline.py` | Full reproducible pipeline code |
| `round2_tool_inventory.json` | Step 0 output |
| `round2_benchmark_validation.json` | Step 1 — Spearman validation |
| `round2_bfactor_analysis.json` | Step 2 — B-factor ranking |
| `round2_site_scan.json` | Step 3 — 450 single-site predictions |
| `round2_candidates.json` | Step 4 — ranked candidates |
| `round2_md_status.json` | Step 5 — BLOCKED |
| `round2_prior_art_protocol.json` | Step 6 — search instructions |
| `round2_final_report.md` | This report |
| `petase_skeptical_audit_report.md` | Round 1 audit (DO NOT OVERWRITE) |

---

*Report generated by autonomous skeptical computational biology agent. All hard rules from the audit protocol enforced. No numbers invented or approximated without explicit labeling.*
