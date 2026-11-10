# Phase G: Scientific Claim-Audit and Human Verification Checklist

**Audit Date:** November 10, 2026  
**Auditor:** Skeptical Computational Biology Self-Auditing Engine  
**Target Repository:** `github.com/Sankar7567/research_pet`

---

## 1. Quantitative Claim-Audit Matrix

Every quantitative claim presented in `manuscript_audit_petase.tex` (and compiled in `manuscript_audit_petase.pdf`) was programmatically checked against verifiable machine-generated logs on disk.

| # | Specific Quantitative Claim | Numerical Value | Source File on Disk | Verdict | Notes |
|:--|:---|:---:|:---|:---:|:---|
| 1 | Round 1 headline $T_m$ fabrication | $88.70^\circ\text{C}$ ($\Delta T_m = +39.88^\circ\text{C}$) | `RETRACTED_manuscript_nature_catalysis.tex` | **RETRACTED** | Confirmed heuristic ungrounded claim |
| 2 | Round 1 headline $\Delta\Delta G$ fabrication | $-15.15\text{ kcal/mol}$ | `compute_thermodynamics.py`, `RETRACTED_manuscript_nature_catalysis.tex` | **RETRACTED** | Confirmed linear heuristic sum |
| 3 | Round 1 active site RMSD bug | $0.0000\text{ \AA}$ | `run_structural_validation.py`, `audit_metrics.json` | **RETRACTED** | Confirmed template identity artifact |
| 4 | N288P wild-type $\phi$ dihedral angle | $+58.4^\circ$ | `audit_metrics.json`, `audit_pipeline_reproduce.py` | **VERIFIED** | Computed from `6EQE.pdb` |
| 5 | N288P steric clash with Ala287(C) | $1.08\text{ \AA}$ | `audit_metrics.json` | **VERIFIED** | Distance below $2.00\text{ \AA}$ VdW cutoff |
| 6 | Unclosed proline $C_\delta\text{--}N$ distance | $3.14\text{ \AA}$ | `candidate_petase_structure.pdb`, `audit_metrics.json` | **VERIFIED** | Physical bond length is $1.47\text{ \AA}$ |
| 7 | Disulfide $C_\alpha\text{--}C_\alpha$ distances | $7.18\text{ \AA}$ and $7.55\text{ \AA}$ | `audit_metrics.json` | **VERIFIED** | Ideal range is $5.0\text{--}6.5\text{ \AA}$ |
| 8 | Disulfide $C_\beta\text{--}S_\gamma\text{--}S_\gamma$ angles | $73.2^\circ$ and $85.1^\circ$ | `audit_metrics.json` | **VERIFIED** | Ideal angle is $103^\circ \pm 5^\circ$ |
| 9 | Round 1 ESM-2 cumulative LLR | $-15.571$ | `esm2_audit_results.json` | **VERIFIED** | Evaluated on `facebook/esm2_t6_8M` |
| 10 | S103C single-mutation ESM-2 LLR | $-7.813$ | `esm2_audit_results.json` | **VERIFIED** | Severe evolutionary penalty |
| 11 | Benchmark correlation (seen backgrounds) | $\rho = -0.8929, p = 0.0068$ ($n=7$) | `round2_benchmark_validation.json` | **VERIFIED** | Tested on 7 benchmark variants |
| 12 | Held-out validation sample size | $n = 8$ | `phase_c_validation_results.json` | **VERIFIED** | Primary literature curation |
| 13 | ThermoMPNN held-out Spearman $\rho$ | $+0.0238$ | `phase_c_validation_results.json` | **VERIFIED** | $10,000$ bootstrap resamples |
| 14 | ThermoMPNN $95\%$ bootstrap CI | $[-0.6790, +0.9231]$ | `phase_c_validation_results.json` | **VERIFIED** | Percentile interval spans zero |
| 15 | ESM-2 LLR held-out Spearman $\rho$ | $-0.0238$ | `phase_c_validation_results.json` | **VERIFIED** | $10,000$ bootstrap resamples |
| 16 | ESM-2 LLR $95\%$ bootstrap CI | $[-0.8481, +0.8987]$ | `phase_c_validation_results.json` | **VERIFIED** | Spans zero symmetrically |
| 17 | Pre-registered gate outcome | `AUDIT-ONLY` | `phase_c_validation_results.json` | **VERIFIED** | Failed $\rho \leq -0.50$ requirement |
| 18 | `thresholds_v3.json` SHA256 hash | `32cb1429033bba8e2f9ed011eb89ee69...` | `thresholds_v3.json` | **VERIFIED** | Verified by `sha256sum` |
| 19 | Compiled Manuscript PDF Size | 9 pages, 213,856 bytes | `manuscript_audit_petase.pdf` | **VERIFIED** | Compiled with Tectonic 0.17.0 |

---

## 2. Mandatory Human Verification Checklist Before Submission

> **CRITICAL ACTIONS REQUIRED BY HUMAN RESEARCHER (K. K. S. B. Sankar):**

- [ ] **1. Primary Source DOI Verification**:
  Every citation in the manuscript was verified against live DOI resolvers. Manually confirm that all volume, page, and year numbers match your institutional library holdings:
  - Son et al. 2019 (*ACS Catalysis*): DOI `10.1021/acscatal.9b00568`
  - Joo et al. 2018 (*Nature Communications*): DOI `10.1038/s41467-018-02881-1`
  - Austin et al. 2018 (*PNAS*): DOI `10.1073/pnas.1718804115`
  - Cui et al. 2021 (*ACS Catalysis*): DOI `10.1021/acscatal.0c05126`
  - Lu et al. 2022 (*Nature*): DOI `10.1038/s41586-022-04599-z`
  - Bell et al. 2022 (*Nature Catalysis*): DOI `10.1038/s41929-022-00828-5`
  - Meng et al. 2021 (*Int. J. Biol. Macromol.*): DOI `10.1016/j.ijbiomac.2021.03.058`
  - Ding et al. 2024 (*Biochem. Eng. J.*): DOI `10.1016/j.bej.2024.109420`
  - Brott et al. 2022 (*Eng. Life Sci.*): DOI `10.1002/elsc.202100105`
  - Dieckhaus et al. 2024 (*PNAS*): DOI `10.1073/pnas.2314853121`
  - Lin et al. 2023 (*Science*): DOI `10.1126/science.ade2574`

- [ ] **2. Institutional Affiliation & Approvals**:
  Confirm that your affiliation details (*Department of Computer Science and Engineering, Srinivasa Institute of Engineering and Technology, Cheyyeru, Amalapuram, Andhra Pradesh 533216, India*) and official institutional email comply with institutional publication policies.

- [ ] **3. Prior-Art & Legal Status**:
  Ensure that no patent or commercial exclusivity claims are made. The manuscript explicitly frames all variants as public-domain audit benchmarks, and notes that sequence queries in Lens.org/Google Patents do not constitute formal legal freedom-to-operate (FTO).

- [ ] **4. Journal Policy & Pre-Print Posting**:
  Check target journal pre-print policies. We recommend depositing the compiled manuscript (`manuscript_audit_petase.pdf`) to **bioRxiv** or **arXiv (q-bio.BM)** prior to formal peer-review submission.

- [ ] **5. Faculty / Peer Review**:
  Request an independent internal review by a biochemistry or bioinformatics faculty colleague at your institute to verify that the negative-result audit framing and tone meet professional standards.
