# Adversarial Self-Auditing PETase Engineering Pipeline & Forensic Audit

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python: 3.11](https://img.shields.io/badge/Python-3.11-blue.svg)](https://www.python.org/)
[![OpenMM: 8.6.1](https://img.shields.io/badge/OpenMM-8.6.1-green.svg)](https://openmm.org/)
[![Tectonic: 0.17.0](https://img.shields.io/badge/Tectonic-0.17.0-purple.svg)](https://tectonic-typesetting.github.io/)

This repository documents an adversarial forensic self-audit of an autonomous AI-agent protein design pipeline tasked with engineering hyper-thermostable variants of *Ideonella sakaiensis* PETase (IsPETase).

---

## ⚠️ Formal Retraction Notice (Round 1 Claims)

In early autonomous iterations, the pipeline claimed the discovery of an engineered candidate designated **`PROD-PETase-G001-C04`** with:
- Melting temperature: $T_m = 88.70^\circ\text{C}$ ($\Delta T_m = +39.88^\circ\text{C}$)
- Folding free energy change: $\Delta\Delta G = -15.15\text{ kcal/mol}$
- Catalytic triad RMSD: $0.0000\text{ \AA}$

**ALL THREE HEADLINE CLAIMS ARE RETRACTED AS SOFTWARE HALLUCINATIONS.**  
Our forensic audit demonstrated that:
1. **$0.0000\text{ \AA}$ RMSD**: An artifact of the PDB coordinate generator copying template active-site atom records verbatim without modification.
2. **$\Delta\Delta G = -15.15\text{ kcal/mol}$**: Derived by summing crude hardcoded heuristic scalars (e.g. $-1.5\text{ kcal/mol}$ per proline) rather than evaluating an atomic force field.
3. **Steric Impossibilities**: Introduced unclosed proline rings ($C_\delta\text{--}N = 3.14\text{ \AA}$), an $\alpha_L$ positive-$\phi$ proline ($N288P$, $\phi = +58.4^\circ$) causing a lethal $1.08\text{ \AA}$ backbone collision, and distorted disulfides ($73^\circ\text{--}85^\circ$ angles).

All fabricated legacy manuscripts have been renamed with `RETRACTED_` prefixes and preserved exclusively for audit reproducibility. See [`petase_skeptical_audit_report.md`](petase_skeptical_audit_report.md).

---

## 🔬 Key Scientific Finding: Circular Benchmark Leakage

When evaluating machine learning stability predictors on IsPETase:
- **On canonical benchmark variants** ($n=7$, including mutations present in the FAST/Thermo backgrounds): ThermoMPNN displayed strong apparent rank correlation:
  $$\text{Spearman } \rho = -0.8929 \quad (p = 0.0068)$$
- **On a strictly held-out literature validation set** ($n=8$, excluding all mutations present in design backgrounds): Predictive power completely collapsed:
  $$\text{Spearman } \rho = +0.0238 \quad (95\%\text{ CI: } [-0.6790, \, +0.9231], \, B=10,000)$$

Because the pre-registered pass criterion ($|\rho| \geq 0.50$, CI excluding zero, expected negative sign) failed, the pre-registered decision rule locked the project outcome to **`AUDIT-ONLY`**.

---

## 📜 Cryptographic Pre-Registration Ledgers

All gates were locked and hashed prior to data analysis:
- **Round 1 Audit**: [`thresholds.json`](thresholds.json)  
  `SHA256: fa8960b2d115b7dab2d20d012fb58ea2699efef3e50908c35d17a6e5b4157053`
- **Round 2 Baseline**: [`thresholds_v2.json`](thresholds_v2.json)  
  `SHA256: 8dba6a9685b8090c08a3fa06907062654478344da4dec3a085f4bcd59acaebbe`
- **Phase B/C Self-Audit Protocol**: [`thresholds_v3.json`](thresholds_v3.json)  
  `SHA256: 32cb1429033bba8e2f9ed011eb89ee69fdd5f9faa4bd26c37d4a35a6260ff9ee`

---

## 📁 Repository Structure

```text
├── manuscript_audit_petase.tex      # Full 9-page publication manuscript (LaTeX)
├── manuscript_audit_petase.pdf      # Compiled PDF (Tectonic 0.17.0)
├── cover_letters.md                 # Tailored cover letters for 3 candidate journals
├── claim_audit_and_verification.md  # 19-point claim-to-file traceability ledger
├── environment.yml                  # Pinned Conda/Python environment specification
├── thresholds_v3.json               # Pre-registered protocol & gates
├── run_phase_c_validation.py        # Bootstrap statistical evaluation engine
├── generate_manuscript_figures.py   # Publication figure generation script
├── figures/
│   ├── figure1_validation_comparison.pdf
│   └── figure2_bootstrap_distributions.pdf
├── round2_pipeline.py               # Reproducible Round 2 evidence pipeline
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
uv pip install --system-certs openmm mdtraj pdbfixer numpy scipy biopython torch transformers omegaconf pandas
```

### 2. Run retrospective validation
```bash
python3 run_phase_c_validation.py
```

### 3. Compile the manuscript
```bash
tectonic manuscript_audit_petase.tex
```

---

## 🌐 Zenodo DOI Archival Plan

Prior to formal journal submission, this repository will be archived on [Zenodo](https://zenodo.org/) via GitHub integration to mint a permanent, immutable Digital Object Identifier (DOI) for all scripts, PDB coordinates, and raw validation datasets.

---

## ⚖️ Citation & AI Disclosure

```bibtex
@article{sankar2026petaseaudit,
  author    = {Sankar, K. K. S. B.},
  title     = {Forensic Self-Audit of an Autonomous AI-Agent Protein Design Pipeline: Root-Cause Failure Analysis of Hallucinated PETase Stability and the Fallacy of Circular Benchmark Validation},
  journal   = {Preprint / Self-Audit Archive},
  year      = {2026},
  url       = {https://github.com/Sankar7567/research_pet}
}
```

*AI Disclosure: An AI agent (Google Antigravity platform) assisted in code execution, simulation setup, and draft assembly. The human corresponding author (K. K. S. B. Sankar) designed the self-audit framework, established all pre-registration ledgers, independently verified all citations and DOIs, and assumes full scientific accountability.*
