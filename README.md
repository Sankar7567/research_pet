# Adversarial Self-Auditing PETase Engineering Pipeline & Forensic Audit

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python: 3.11](https://img.shields.io/badge/Python-3.11-blue.svg)](https://www.python.org/)
[![OpenMM: 8.6.1](https://img.shields.io/badge/OpenMM-8.6.1-green.svg)](https://openmm.org/)
[![Tectonic: 0.17.0](https://img.shields.io/badge/Tectonic-0.17.0-purple.svg)](https://tectonic-typesetting.github.io/)

This repository documents an adversarial forensic self-audit and pre-registered negative-result case study of an autonomous AI-agent protein design pipeline tasked with engineering hyper-thermostable variants of *Ideonella sakaiensis* PETase (IsPETase).

---

## ⚠️ Formal Retraction Notice (Round 1 Claims)

In early autonomous iterations, the generative pipeline claimed the design of an engineered candidate designated **`PROD-PETase-G001-C04`** with:
- Melting temperature: $T_m = 88.70^\circ\text{C}$ ($\Delta T_m = +39.88^\circ\text{C}$)
- Folding free energy change: $\Delta\Delta G = -15.15\text{ kcal/mol}$
- Catalytic triad RMSD: $0.0000\text{ \AA}$

**ALL THREE HEADLINE CLAIMS ARE RETRACTED AS SOFTWARE HALLUCINATIONS.**  
Our forensic audit demonstrated that:
1. **$0.0000\text{ \AA}$ RMSD**: An artifact of the PDB coordinate generator copying template active-site atom records verbatim without energetic relaxation (actual relaxed shift is $0.0082\text{ \AA}$; overall backbone RMSD is $0.0984\text{ \AA}$).
2. **$\Delta\Delta G = -15.15\text{ kcal/mol}$**: Derived by summing crude hardcoded heuristic scalars ($-1.5\text{ kcal/mol}$ per proline, $-2.5$ per disulfide) rather than evaluating an atomic force field or physical potential.
3. **Steric Impossibilities**: Introduced unclosed proline rings ($C_\delta\text{--}N = 3.14\text{ \AA}$ vs $1.47\text{ \AA}$), an $\alpha_L$ positive-$\phi$ proline ($N288P$, $\phi = +58.4^\circ$, $\psi = +33.1^\circ$) causing a lethal $1.08\text{ \AA}$ backbone collision with Gly287, and distorted disulfides ($73^\circ\text{--}85^\circ$ angles).

All fabricated legacy files have been prefixed with `RETRACTED_` and preserved exclusively for audit traceability. See [`petase_skeptical_audit_report.md`](petase_skeptical_audit_report.md).

---

## 🔬 Retrospective Validation on Verified Primary Literature

Following the retraction of Round 1 claims, we established a pre-registered protocol requiring retrospective validation of candidate scoring models (ThermoMPNN and ESM-2) against primary peer-reviewed literature data:
- **On benchmark variants** ($n=5$ mutations present in the FAST/Thermo backgrounds: S121E, D186H, R280A, N233K, R224Q): ThermoMPNN displayed strong rank correlation:
  $$\text{Spearman } \rho = -0.8000 \quad (95\%\text{ CI: } [-1.0000, \, -0.1111], \, B=10,000)$$
- **On primary held-out stability mutants** ($n=5$ non-knockout single mutants: P181A, W185A, F229Y, K95A, R132N):
  $$\text{Spearman } \rho = -0.6000 \quad (95\%\text{ CI: } [-1.0000, \, +1.0000], \, B=10,000)$$
- **On held-out set with functional knockout controls** ($n=7$, adding S160A catalytic knockout and C203A/C239A disulfide knockout):
  $$\text{Spearman } \rho = -0.0714 \quad (95\%\text{ CI: } [-0.7647, \, +1.0000], \, B=10,000)$$
- **On combined literature dataset** ($n=12$):
  $$\text{Spearman } \rho = -0.4685 \quad (95\%\text{ CI: } [-0.8750, \, +0.1643], \, B=10,000)$$

### Statistical Power Constraint & Inconclusive Interpretation
Analytical power calculations indicate that sample sizes of $n=5\text{--}7$ only have $80\%$ power at $\alpha=0.05$ to detect extreme effects ($|\rho| \geq 0.886\text{--}0.963$). Because confidence intervals span zero, the retrospective validation is **statistically inconclusive** on novel single mutants. In strict accordance with the pre-registered decision rule, candidate promotion was halted, locking the outcome to **`AUDIT-ONLY`**.

---

## 📜 Hash-Timestamped Pre-Registration Ledgers

All operational constraints, pass/fail gates, and decision trees were version-controlled and hashed before executing analyses:
- **Round 1 Audit**: [`thresholds.json`](thresholds.json)  
  `SHA256: fa8960b2d115b7dab2d20d012fb58ea2699efef3e50908c35d17a6e5b4157053`
- **Round 2 Baseline**: [`thresholds_v2.json`](thresholds_v2.json)  
  `SHA256: 8dba6a9685b8090c08a3fa06907062654478344da4dec3a085f4bcd59acaebbe`
- **Phase B/C Self-Audit Protocol**: [`thresholds_v3.json`](thresholds_v3.json)  
  `SHA256: 32cb1429033bba8e2f9ed011eb89ee69fdd5f9faa4bd26c37d4a35a6260ff9ee`  
  *Git Commit Timestamp:* `2026-11-10T13:01:59Z` (Commit `1fb0754`)

---

## 📁 Repository Structure

```text
├── manuscript.pdf                   # Submission-ready bioRxiv preprint (Root directory)
├── manuscript.tex                   # Complete LaTeX source (11pt, single column, lineno)
├── data_verification_log.csv        # Primary literature verification log (16 variants)
├── claim_audit.csv                  # 35-point claim-to-file traceability ledger
├── phase_c_validation_results.json  # Full statistical output and bootstrap resamples
├── run_phase_c_validation.py        # Bootstrap statistical evaluation engine
├── generate_manuscript_figures.py   # Publication vector figure generator
├── figures/
│   ├── figure1_pipeline_architecture.pdf
│   ├── figure2_audit_failure_modes.pdf
│   ├── figure3_retrospective_validation.pdf
│   └── figure4_power_and_bootstrap.pdf
├── cover_letters.md                 # Tailored cover letters for candidate journals
├── environment.yml                  # Pinned Conda/Python environment specification
├── thresholds_v3.json               # Pre-registered protocol & gates
├── round2_pipeline.py               # Round 2 evidence pipeline
├── petase_skeptical_audit_report.md # Round 1 forensic audit report
├── round2_final_report.md           # Round 2 design evaluation report
└── LICENSE                          # MIT License
```

---

## 🛠️ Environment Setup & Reproduction

### 1. Install dependencies
```bash
# Using uv (recommended)
uv venv env_openmm --python 3.11
uv pip install --system-certs openmm mdtraj pdbfixer numpy scipy biopython torch transformers omegaconf
```

### 2. Run retrospective validation
```bash
python3 run_phase_c_validation.py
```

### 3. Generate manuscript figures
```bash
python3 generate_manuscript_figures.py
```

### 4. Compile the preprint PDF
```bash
tectonic manuscript.tex
```
Output compiles directly to `manuscript.pdf` in the workspace root.

---

## 🌐 Zenodo DOI Archival Plan

Prior to formal bioRxiv preprint submission, this repository will be archived on [Zenodo](https://zenodo.org/) via GitHub integration to mint a permanent, immutable Digital Object Identifier (DOI) for all scripts, PDB coordinates, figures, and raw validation datasets.
