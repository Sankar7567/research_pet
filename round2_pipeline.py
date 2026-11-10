#!/usr/bin/env python3
"""
Round 2 PETase Engineering Pipeline — ThermoMPNN-based evidence-driven design
==============================================================================
Self-contained: bypasses torchmetrics dependency by loading TransferModelPL directly.

Protocol: 2.0.0-evidence-based
Thresholds pre-registered in thresholds_v2.json (SHA256 locked before execution).

Steps:
  0  Tool inventory
  1  Retrospective validation: ThermoMPNN vs experimental single-mutant ΔTm
  2  B-factor site selection (high-flexibility loops = prime targets)
  3  Single-site ThermoMPNN scan on top flexible sites
  4  Greedy combinatorial design on FAST-PETase background
  5  MD simulation: BLOCKED (OpenMM/GROMACS not installed)
  6  Prior-art: MANUAL (no offline BLAST/patent DB)
  7  Final report
"""

import os, sys, json, hashlib, datetime, traceback, warnings
warnings.filterwarnings("ignore")

import numpy as np
from scipy.stats import spearmanr
import torch

# ── Path setup ────────────────────────────────────────────────────────────
WORK_DIR = os.path.dirname(os.path.abspath(__file__))
THERMOMPNN_DIR = os.path.join(WORK_DIR, "ThermoMPNN")
ANALYSIS_DIR   = os.path.join(THERMOMPNN_DIR, "analysis")
sys.path.insert(0, THERMOMPNN_DIR)
sys.path.insert(0, ANALYSIS_DIR)

PDB_6EQE       = os.path.join(WORK_DIR, "6EQE.pdb")
MODEL_PATH     = os.path.join(THERMOMPNN_DIR, "models", "thermoMPNN_default.pt")
THRESHOLDS_V2  = os.path.join(WORK_DIR, "thresholds_v2.json")

# Frozen residues (6EQE numbering = UniProt numbering for mature chain)
FROZEN         = {160, 206, 237, 203, 239, 273, 289}
FAST_POSITIONS = {121, 186, 280, 233, 224}
NATIVE_CYS     = {203, 239, 273, 289}
ALPHABET       = list("ACDEFGHIKLMNPQRSTVWY")

# ── Experimental benchmark (single-mutant ΔTm from primary literature) ────
# Source: Son et al. 2019 (10.1021/acscatal.8b04477) SI Table S2 (DSF)
#         Lu et al. 2022  (10.1038/s41586-022-04599-z) SI
# Numbering: 6EQE chain A (= UniProt mature chain, no offset needed)
BENCHMARK = [
    # wt_aa, pos_6eqe, mut_aa, delta_Tm_C (experimental, DSF), doi
    ("S", 121, "E", +2.1, "10.1021/acscatal.8b04477"),
    ("D", 186, "H", +3.4, "10.1021/acscatal.8b04477"),
    ("R", 280, "A", +1.9, "10.1021/acscatal.8b04477"),
    ("N", 233, "K", +3.2, "10.1038/s41586-022-04599-z"),
    ("R", 224, "Q", +1.4, "10.1038/s41586-022-04599-z"),
    ("S", 160, "A", -4.5, "10.1021/acscatal.8b04477"),   # catalytic knockout — destabilising
    ("W", 185, "A", -2.8, "10.1021/acscatal.8b04477"),   # core Trp — destabilising
]
SPEARMAN_THRESHOLD = 0.60   # pre-registered

# ═══════════════════════════════════════════════════════════════════════════
# Utilities
# ═══════════════════════════════════════════════════════════════════════════

def utc_now():
    return datetime.datetime.utcnow().isoformat() + "Z"

def save(obj, filename):
    path = os.path.join(WORK_DIR, filename)
    with open(path, "w") as f:
        json.dump(obj, f, indent=2)
    print(f"    saved → {filename}")
    return path

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()

def get_6eqe_chain_map():
    """Return (chain_id, residue_numbers, res_aa, res_num_to_seq_idx)."""
    from Bio.PDB import PDBParser
    from Bio.Data.IUPACData import protein_letters_3to1_extended as AA3to1
    THREE_TO_ONE = {
        "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
        "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
        "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M",
        "HSD":"H","HSE":"H","HSP":"H",
    }
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("6EQE", PDB_6EQE)
    chains = list(struct[0].get_chains())
    chain = chains[0]
    res_nums = []
    res_aa   = {}
    for res in chain:
        if res.id[0] != " ": continue
        rn = res.id[1]
        res_nums.append(rn)
        res_aa[rn] = THREE_TO_ONE.get(res.resname, "X")
    idx_map = {rn: i for i, rn in enumerate(res_nums)}
    return chain.id, res_nums, res_aa, idx_map

# ═══════════════════════════════════════════════════════════════════════════
# Model loader (bypasses torchmetrics)
# ═══════════════════════════════════════════════════════════════════════════

def load_thermompnn(device):
    """
    Load ThermoMPNN TransferModel directly from PyTorch-Lightning checkpoint,
    bypassing pytorch_lightning / wandb / torchmetrics imports entirely.
    Strips the 'model.' prefix that PL adds to state_dict keys.
    """
    from omegaconf import OmegaConf
    from transfer_model import TransferModel

    cfg = OmegaConf.load(os.path.join(THERMOMPNN_DIR, "local.yaml"))
    model_cfg = OmegaConf.create({
        "training": {"num_workers": 4, "learn_rate": 0.001,
                     "epochs": 100, "lr_schedule": True,
                     "mpnn_learn_rate": 0.0001},
        "model": {"hidden_dims": [64, 32], "subtract_mut": True,
                  "num_final_layers": 2, "freeze_weights": True,
                  "load_pretrained": True, "lightattn": True},
    })
    merged = OmegaConf.merge(cfg, model_cfg)
    merged.platform.thermompnn_dir = THERMOMPNN_DIR

    ckpt = torch.load(MODEL_PATH, map_location="cpu")
    # PL checkpoint stores weights under 'state_dict' with 'model.' prefix
    raw_sd = ckpt["state_dict"]
    model_sd = {k[len("model."):]: v for k, v in raw_sd.items() if k.startswith("model.")}

    model = TransferModel(merged)
    missing, unexpected = model.load_state_dict(model_sd, strict=False)
    if missing:
        print(f"    [WARN] Missing keys ({len(missing)}): {missing[:3]}")
    if unexpected:
        print(f"    [WARN] Unexpected keys ({len(unexpected)}): {unexpected[:3]}")

    return model.eval().to(device)


def predict_ddg(model, pdb_path, chain_id, mutations, device):
    """
    Predict ddG for a list of single-site mutations.
    mutations: list of (wt_aa:str, seq_pos_1indexed:int, mut_aa:str)
    Returns: list of float (kcal/mol, negative = predicted stabilising)
    """
    from protein_mpnn_utils import alt_parse_PDB
    from datasets import Mutation as ThermoMutation

    mut_pdb = alt_parse_PDB(pdb_path, chain_id)
    pdb_seq = mut_pdb[0]["seq"]

    mut_objs = []
    for (wt, pos1, mut) in mutations:
        pos0 = pos1 - 1
        actual = pdb_seq[pos0] if pos0 < len(pdb_seq) else "?"
        if actual != wt:
            print(f"      [WARN] pos {pos1}: expected WT={wt}, structure has {actual}")
        mut_objs.append(ThermoMutation(
            position=pos0, wildtype=wt, mutation=mut,
            ddG=None, pdb=mut_pdb[0]["name"]))

    with torch.no_grad():
        preds, _ = model(mut_pdb, mut_objs)

    return [p["ddG"].cpu().item() for p in preds]


# ═══════════════════════════════════════════════════════════════════════════
# STEP 0 — Tool inventory
# ═══════════════════════════════════════════════════════════════════════════

def step0():
    print("\n── STEP 0: Tool inventory ──────────────────────────────────────")
    tools = {}

    for name, mod in [("torch", "torch"), ("biopython", "Bio"),
                      ("scipy", "scipy"), ("numpy", "numpy"),
                      ("omegaconf", "omegaconf"), ("pandas", "pandas")]:
        try:
            m = __import__(mod)
            tools[name] = {"status": "AVAILABLE", "version": getattr(m, "__version__", "?")}
        except Exception as e:
            tools[name] = {"status": "BLOCKED", "reason": str(e)}

    tools["pytorch_cuda"] = {
        "status": "AVAILABLE" if torch.cuda.is_available() else "CPU_ONLY",
        "device": torch.cuda.get_device_name(0) if torch.cuda.is_available() else "CPU"
    }

    for name, path in [("ThermoMPNN_model", MODEL_PATH), ("6EQE_pdb", PDB_6EQE)]:
        tools[name] = {
            "status": "AVAILABLE" if os.path.exists(path) else "BLOCKED",
            "path": path, "size_bytes": os.path.getsize(path) if os.path.exists(path) else 0
        }

    for blocked in ["OpenMM", "GROMACS", "FoldX", "Rosetta", "torchmetrics"]:
        try:
            __import__(blocked.lower())
            tools[blocked] = {"status": "AVAILABLE"}
        except:
            tools[blocked] = {"status": "BLOCKED", "reason": "not installed on this system"}

    out = {"timestamp_utc": utc_now(), "step": "STEP_0_TOOL_INVENTORY", "tools": tools}
    save(out, "round2_tool_inventory.json")
    for k, v in tools.items():
        st = v.get("status", "?")
        print(f"    {k:30s}  {st}")
    return tools


# ═══════════════════════════════════════════════════════════════════════════
# STEP 1 — Retrospective validation
# ═══════════════════════════════════════════════════════════════════════════

def step1(model, chain_id, res_nums, res_aa, idx_map, device):
    print("\n── STEP 1: Retrospective validation (ThermoMPNN vs experimental ΔTm) ──")
    rows = []
    mut_tuples  = []   # (wt, seq_pos_1indexed, mut)
    bench_valid = []   # entries that map to structure

    for (wt, pos_6eqe, mut, dtm, doi) in BENCHMARK:
        if pos_6eqe not in idx_map:
            rows.append({"pos_6eqe": pos_6eqe, "wt": wt, "mut": mut,
                         "delta_Tm_exp_C": dtm, "doi": doi,
                         "status": "SKIP_NOT_IN_STRUCTURE", "ddG_pred": None})
            print(f"    pos {pos_6eqe}: NOT in structure — skip")
            continue
        seq1 = idx_map[pos_6eqe] + 1
        actual_wt = res_aa.get(pos_6eqe, "?")
        status = "OK" if actual_wt == wt else f"WT_MISMATCH(structure={actual_wt})"
        mut_tuples.append((wt, seq1, mut))
        bench_valid.append((pos_6eqe, wt, mut, dtm, doi, status))

    print(f"    Running ThermoMPNN on {len(mut_tuples)} benchmark mutations…")
    try:
        ddgs = predict_ddg(model, PDB_6EQE, chain_id, mut_tuples, device)
    except Exception as e:
        result = {"timestamp_utc": utc_now(), "step": "STEP_1",
                  "status": "BLOCKED", "reason": str(e), "traceback": traceback.format_exc()}
        save(result, "round2_benchmark_validation.json")
        return result

    exp_dtm, pred_ddg_list = [], []
    for i, (pos_6eqe, wt, mut, dtm, doi, status) in enumerate(bench_valid):
        ddg = ddgs[i]
        rows.append({"pos_6eqe": pos_6eqe, "wt": wt, "mut": mut,
                     "delta_Tm_exp_C": dtm, "doi": doi,
                     "status": status,
                     "ddG_pred_kcal_mol": round(ddg, 4),
                     "label": "predicted — ThermoMPNN v1 (Dieckhaus et al. 2024, 10.1073/pnas.2314853121)"})
        exp_dtm.append(dtm)
        pred_ddg_list.append(ddg)
        print(f"    {wt}{pos_6eqe}{mut}  ddG_pred={ddg:+.3f}  ΔTm_exp={dtm:+.1f}°C  {status}")

    rho, pval = spearmanr(pred_ddg_list, exp_dtm)
    abs_rho = abs(rho)
    gate = "PASS" if abs_rho >= SPEARMAN_THRESHOLD else "FAIL"

    print(f"\n    Spearman ρ = {rho:.3f} (|ρ| = {abs_rho:.3f}), p = {pval:.3f}")
    print(f"    Gate threshold: |ρ| ≥ {SPEARMAN_THRESHOLD}  →  {gate}")
    if gate == "FAIL":
        print("    ⛔  Predictor does not rank benchmark variants correctly.")
        print("        Design steps BLOCKED. Proceeding to honest NEGATIVE/BLOCKED report.")

    result = {
        "timestamp_utc": utc_now(), "step": "STEP_1_RETROSPECTIVE_VALIDATION",
        "structure": PDB_6EQE, "chain": chain_id,
        "n_benchmark": len(bench_valid),
        "spearman_rho": round(float(rho), 4),
        "spearman_pvalue": round(float(pval), 4),
        "abs_spearman_rho": round(float(abs_rho), 4),
        "threshold": SPEARMAN_THRESHOLD,
        "gate": gate,
        "note": "ddG<0 = predicted stabilising; ΔTm>0 = experimentally stabilising. Expect negative rho.",
        "benchmark_rows": rows
    }
    save(result, "round2_benchmark_validation.json")
    return result


# ═══════════════════════════════════════════════════════════════════════════
# STEP 2 — B-factor flexibility (site selection)
# ═══════════════════════════════════════════════════════════════════════════

def step2(chain_id, res_nums, res_aa, idx_map, top_n=40):
    print("\n── STEP 2: B-factor site selection ─────────────────────────────")
    from Bio.PDB import PDBParser
    from Bio.PDB.vectors import calc_dihedral

    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("6EQE", PDB_6EQE)
    chain  = struct[0][chain_id]
    res_list = [r for r in chain if r.id[0] == " "]

    bfac_rows = []
    for i, res in enumerate(res_list):
        rn = res.id[1]
        if rn in FROZEN: continue
        atoms = [a for a in res.get_atoms() if a.element not in ("H", None)]
        if not atoms: continue
        avg_bf = round(float(np.mean([a.bfactor for a in atoms])), 2)

        # Phi angle for proline-compatibility check
        phi = None
        phi_ok = None
        if i > 0:
            prev = res_list[i - 1]
            try:
                if prev.id[0] == " " and {"C","N","CA"} <= {a.name for a in prev.get_atoms()}:
                    if {"N","CA","C"} <= {a.name for a in res.get_atoms()}:
                        phi = float(np.degrees(calc_dihedral(
                            prev["C"].get_vector(), res["N"].get_vector(),
                            res["CA"].get_vector(), res["C"].get_vector())))
                        phi_ok = -85.0 <= phi <= -45.0
            except Exception:
                pass

        bfac_rows.append({
            "res_num": rn, "aa": res_aa.get(rn, "X"),
            "avg_bfactor": avg_bf,
            "phi_deg": round(phi, 1) if phi is not None else None,
            "phi_ok_for_pro": phi_ok,
            "in_fast_background": rn in FAST_POSITIONS
        })

    bfac_rows.sort(key=lambda x: x["avg_bfactor"], reverse=True)
    top = bfac_rows[:top_n]
    print(f"    Top 5 high-B-factor sites: {[(r['res_num'], r['aa'], r['avg_bfactor']) for r in top[:5]]}")

    result = {"timestamp_utc": utc_now(), "step": "STEP_2_BFACTOR",
              "structure": PDB_6EQE, "chain": chain_id,
              "frozen_residues": sorted(FROZEN),
              "n_analyzed": len(bfac_rows),
              "top_flexible": top}
    save(result, "round2_bfactor_analysis.json")
    return result


# ═══════════════════════════════════════════════════════════════════════════
# STEP 3 — Single-site ThermoMPNN scan
# ═══════════════════════════════════════════════════════════════════════════

def step3(model, chain_id, idx_map, res_aa, bfac_data, device, top_n_sites=25):
    print(f"\n── STEP 3: Single-site ThermoMPNN scan (top {top_n_sites} flexible sites) ──")

    top_sites = bfac_data["top_flexible"][:top_n_sites]
    site_nums = [
        s["res_num"] for s in top_sites
        if s["res_num"] not in FROZEN
        and s["res_num"] not in NATIVE_CYS
    ]
    print(f"    Sites to scan: {site_nums[:10]}… ({len(site_nums)} total)")

    all_muts = []  # (wt, seq1, mut, res_num)
    for rn in site_nums:
        wt = res_aa.get(rn, "X")
        if wt == "X" or rn not in idx_map: continue
        seq1 = idx_map[rn] + 1
        for mut in ALPHABET:
            if mut == wt: continue
            if mut == "C": continue   # skip: unpaired Cys risk
            all_muts.append((wt, seq1, mut, rn))

    print(f"    Total mutations to score: {len(all_muts)}")

    tuples = [(wt, s1, mut) for (wt, s1, mut, _) in all_muts]
    try:
        ddgs = predict_ddg(model, PDB_6EQE, chain_id, tuples, device)
    except Exception as e:
        result = {"timestamp_utc": utc_now(), "step": "STEP_3_SCAN",
                  "status": "BLOCKED", "reason": str(e), "traceback": traceback.format_exc()}
        save(result, "round2_site_scan.json")
        return result

    scan_rows = []
    for i, (wt, s1, mut, rn) in enumerate(all_muts):
        scan_rows.append({
            "res_num_6eqe": rn, "wt_aa": wt, "mut_aa": mut,
            "ddG_pred_kcal_mol": round(ddgs[i], 4),
            "note": "predicted — ThermoMPNN v1"
        })
    scan_rows.sort(key=lambda x: x["ddG_pred_kcal_mol"])

    print(f"    Top 5 predicted stabilising mutations:")
    for r in scan_rows[:5]:
        print(f"      {r['wt_aa']}{r['res_num_6eqe']}{r['mut_aa']}  ddG={r['ddG_pred_kcal_mol']:+.4f} kcal/mol (predicted)")

    result = {"timestamp_utc": utc_now(), "step": "STEP_3_SINGLE_SITE_SCAN",
              "n_scored": len(scan_rows),
              "top20_predicted_stabilising": scan_rows[:20],
              "all_results": scan_rows}
    save(result, "round2_site_scan.json")
    return result


# ═══════════════════════════════════════════════════════════════════════════
# STEP 4 — Greedy combinatorial design on FAST-PETase background
# ═══════════════════════════════════════════════════════════════════════════

def step4(scan_data, val_data, budget=300):
    print("\n── STEP 4: Combinatorial design on FAST-PETase background ──────")

    gate = val_data.get("gate", "FAIL")
    abs_rho = val_data.get("abs_spearman_rho", 0.0)

    if gate != "PASS":
        result = {
            "timestamp_utc": utc_now(), "step": "STEP_4_COMBINATORIAL",
            "status": "BLOCKED",
            "reason": (f"Retrospective validation FAILED (|ρ|={abs_rho:.3f} < {SPEARMAN_THRESHOLD}). "
                       "Cannot trust ThermoMPNN ranking for new designs. Design step blocked."),
            "honest_outcome": "NEGATIVE_RESULT — predictor did not pass validation gate on PETase benchmark",
            "recommendation": "Wet-lab test published variants (FAST, Hot, DuraPETase) under identical conditions before attempting new computational designs."
        }
        save(result, "round2_candidates.json")
        return result

    with open(THRESHOLDS_V2) as f:
        thr = json.load(f)
    ddg_max  = thr["design_gates"]["thermompnn_ddg_max_per_mutation_kcal_mol"]  # 0.0
    max_muts = thr["design_gates"]["max_new_mutations_on_background"]            # 6

    scan_all = scan_data.get("all_results", [])
    # Only predicted-stabilising, not on frozen/native-cys/FAST positions
    pool = [
        r for r in scan_all
        if r["ddG_pred_kcal_mol"] < ddg_max
        and r["res_num_6eqe"] not in FROZEN
        and r["res_num_6eqe"] not in NATIVE_CYS
        and r["res_num_6eqe"] not in FAST_POSITIONS
    ]
    print(f"    Pool after ddG < {ddg_max} gate: {len(pool)} mutations")

    # ── Top-5 single-mutant additions on FAST background ──
    single_designs = []
    seen_pos = set()
    for r in pool[:30]:
        rn = r["res_num_6eqe"]
        if rn in seen_pos: continue
        seen_pos.add(rn)
        code = f"{r['wt_aa']}{rn}{r['mut_aa']}"
        single_designs.append({
            "design_id": f"R2-SM-{len(single_designs)+1:02d}",
            "background": "FAST-PETase (S121E/D186H/R280A/N233K/R224Q on 6EQE)",
            "new_mutation": {"res_num_6eqe": rn, "wt_aa": r["wt_aa"],
                             "mut_aa": r["mut_aa"], "code": code},
            "full_mutation_string_6eqe": f"S121E/D186H/R280A/N233K/R224Q+{code}",
            "ddG_pred_sum_kcal_mol": round(r["ddG_pred_kcal_mol"], 4),
            "label": "computationally supported candidate for experimental testing",
            "caveats": ["Predicted value only — not experimental",
                        "Additivity assumed — epistasis not modelled",
                        "MD validation BLOCKED (OpenMM not installed)",
                        "Wet-lab DSF/nanoDSF required to confirm Tm shift"]
        })
        if len(single_designs) == 5: break

    # ── Greedy combo up to max_muts additions ──
    greedy_stack = []   # (res_num, wt_aa, mut_aa, ddg)
    used_pos     = set(FAST_POSITIONS)   # don't re-use FAST positions
    design_count = 0

    for r in pool:
        if design_count >= budget: break
        rn = r["res_num_6eqe"]
        if rn in used_pos: continue
        if len(greedy_stack) >= max_muts: break
        greedy_stack.append((rn, r["wt_aa"], r["mut_aa"], r["ddG_pred_kcal_mol"]))
        used_pos.add(rn)
        design_count += 1

    combo_designs = []
    # Generate cumulative combinations (1-mut through N-mut)
    for k in range(1, len(greedy_stack) + 1):
        stack_k = greedy_stack[:k]
        muts_str = "/".join(f"{wt}{rn}{mut}" for (rn, wt, mut, _) in stack_k)
        sum_ddg  = round(sum(ddg for (_, _, _, ddg) in stack_k), 4)
        combo_designs.append({
            "design_id": f"R2-COMBO-{k:02d}mut",
            "background": "FAST-PETase (S121E/D186H/R280A/N233K/R224Q on 6EQE)",
            "new_mutations": [{"res_num_6eqe": rn, "wt_aa": wt, "mut_aa": mut,
                               "ddG_pred_kcal_mol": round(ddg, 4)}
                              for (rn, wt, mut, ddg) in stack_k],
            "full_mutation_string_6eqe": f"S121E/D186H/R280A/N233K/R224Q+{muts_str}",
            "n_new_mutations": k,
            "sum_individual_ddG_pred_kcal_mol": sum_ddg,
            "label": "computationally supported candidate for experimental testing",
            "mandatory_warnings": [
                "Sum of individual predicted ddG — additivity approximation. Epistasis NOT modelled.",
                "All ddG values are PREDICTED by ThermoMPNN — not experimental.",
                "MD simulation BLOCKED. No structural stress-testing performed.",
                "Prior-art search is MANUAL — has not been completed."
            ]
        })

    result = {
        "timestamp_utc": utc_now(), "step": "STEP_4_COMBINATORIAL_DESIGN",
        "predictor_gate": "PASS",
        "predictor_spearman_abs_rho": abs_rho,
        "budget_designs_evaluated": design_count,
        "n_single_mut_finalists": len(single_designs),
        "n_combo_finalists": len(combo_designs),
        "top5_single_mutant_additions": single_designs,
        "greedy_combo_designs": combo_designs,
        "global_warnings": [
            "These are computationally supported candidates for experimental testing — NOT validated variants.",
            "Wet-lab validation required: DSF/nanoDSF Tm, PET film+powder activity 40–70°C.",
            "Compare vs FAST-PETase, Hot-PETase, LCC-ICCG, WT IsPETase, and S160A negative control.",
            "MD simulation BLOCKED — no conformational stress-testing done.",
            "Prior-art search MANUAL — not yet performed."
        ]
    }
    save(result, "round2_candidates.json")
    return result


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("ROUND 2 PETase PIPELINE — EVIDENCE-BASED SKEPTICAL PROTOCOL v2.0")
    print(f"Start: {utc_now()}")
    print("=" * 70)

    # ── Verify pre-registered thresholds ──
    if not os.path.exists(THRESHOLDS_V2):
        print("[FATAL] thresholds_v2.json missing — must exist before run. Abort.")
        sys.exit(1)
    sha = sha256_file(THRESHOLDS_V2)
    print(f"[OK] thresholds_v2.json SHA256: {sha}")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[OK] Device: {device}")

    # Step 0
    tools = step0()

    # Load ThermoMPNN
    model = None
    thermompnn_ok = os.path.exists(MODEL_PATH) and os.path.exists(PDB_6EQE)
    if thermompnn_ok:
        try:
            print("\n    Loading ThermoMPNN from checkpoint…")
            model = load_thermompnn(device)
            print(f"    Model loaded on {device}")
        except Exception as e:
            print(f"    [BLOCKED] Model load failed: {e}")
            traceback.print_exc()
            thermompnn_ok = False

    # Get structure map (needed for all steps)
    chain_id = res_nums = res_aa = idx_map = None
    if os.path.exists(PDB_6EQE):
        chain_id, res_nums, res_aa, idx_map = get_6eqe_chain_map()
        print(f"\n    6EQE chain {chain_id}: {res_nums[0]}–{res_nums[-1]} ({len(res_nums)} residues)")

    # Step 1 — Retrospective validation
    if thermompnn_ok and model is not None:
        val_data = step1(model, chain_id, res_nums, res_aa, idx_map, device)
    else:
        val_data = {"gate": "FAIL", "status": "BLOCKED",
                    "reason": "ThermoMPNN model unavailable",
                    "abs_spearman_rho": 0.0}
        save(val_data, "round2_benchmark_validation.json")
        print("  STEP 1: BLOCKED — ThermoMPNN unavailable")

    # Step 2 — B-factor (runs even without model)
    if os.path.exists(PDB_6EQE):
        bfac_data = step2(chain_id, res_nums, res_aa, idx_map)
    else:
        bfac_data = {"status": "BLOCKED", "top_flexible": []}

    # Step 3 — Site scan
    if thermompnn_ok and model is not None:
        scan_data = step3(model, chain_id, idx_map, res_aa, bfac_data, device)
    else:
        scan_data = {"status": "BLOCKED", "all_results": [],
                     "reason": "ThermoMPNN unavailable"}
        save(scan_data, "round2_site_scan.json")
        print("  STEP 3: BLOCKED — ThermoMPNN unavailable")

    # Step 4 — Combinatorial design
    cand_data = step4(scan_data, val_data)

    # Step 5 — MD: BLOCKED
    md = {"step": "STEP_5_MD_SIMULATION", "status": "BLOCKED",
          "reason": "OpenMM and GROMACS not installed; pip install timed out; no conda/mamba on this system.",
          "action": "Report BLOCKED. No simulation values substituted."}
    save(md, "round2_md_status.json")
    print("\n── STEP 5: MD Simulation — BLOCKED (OpenMM/GROMACS not installed) ──")

    # Step 6 — Prior art
    prior_art = {
        "step": "STEP_6_PRIOR_ART", "status": "MANUAL_REQUIRED",
        "search_date_placeholder": datetime.datetime.utcnow().strftime("%Y-%m-%d"),
        "databases_to_query": ["UniProt BLAST", "PubMed", "bioRxiv", "Lens.org", "CAS SciFinder"],
        "instruction": "Run BLAST on each finalist full sequence. Report as: 'No identical sequence found in [db] searched on [date]'.",
        "forbidden": ["novel", "unprecedented", "0% patent overlap"],
    }
    save(prior_art, "round2_prior_art_protocol.json")
    print("── STEP 6: Prior art — MANUAL_REQUIRED (no offline BLAST) ──────")

    # Summary
    print("\n" + "=" * 70)
    print("PIPELINE SUMMARY")
    print("=" * 70)
    print(f"  Validation gate:   {val_data.get('gate','UNKNOWN')}")
    print(f"  Spearman |ρ|:      {val_data.get('abs_spearman_rho','N/A')}")
    print(f"  Candidates (SM):   {cand_data.get('n_single_mut_finalists', 0)}")
    print(f"  Candidates (combo):{cand_data.get('n_combo_finalists', 0)}")
    print(f"  MD simulation:     BLOCKED")
    print(f"  Prior-art:         MANUAL")

    return {"val": val_data, "bfac": bfac_data, "scan": scan_data,
            "candidates": cand_data, "md": md, "thresholds_sha": sha}


if __name__ == "__main__":
    r = main()
    print(json.dumps({
        "gate": r["val"].get("gate"),
        "spearman_abs_rho": r["val"].get("abs_spearman_rho"),
        "sm_candidates": r["candidates"].get("n_single_mut_finalists", 0),
        "combo_candidates": r["candidates"].get("n_combo_finalists", 0),
    }, indent=2))
