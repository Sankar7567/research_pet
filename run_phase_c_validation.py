#!/usr/bin/env python3
"""
Phase C: Held-out Literature Validation & Predictor Evaluation
=============================================================
Evaluates ThermoMPNN and ESM-2 against primary literature experimental IsPETase variants.
Enforces 10,000 bootstrap resamples to calculate 95% Confidence Intervals.
Applies pre-registered decision rule from thresholds_v3.json.

Author: Skeptical Computational Biology Audit Team
Date: 2026-11-10
"""

import os
import sys
import json
import hashlib
import datetime
import numpy as np
from scipy.stats import spearmanr
import torch

# Paths
WORK_DIR = os.path.dirname(os.path.abspath(__file__))
PDB_6EQE = os.path.join(WORK_DIR, "6EQE.pdb")
THRESHOLDS_V3 = os.path.join(WORK_DIR, "thresholds_v3.json")
THERMOMPNN_DIR = os.path.join(WORK_DIR, "ThermoMPNN")
sys.path.insert(0, THERMOMPNN_DIR)
sys.path.insert(0, os.path.join(THERMOMPNN_DIR, "analysis"))

# Verification of thresholds_v3.json SHA256 before doing anything
with open(THRESHOLDS_V3, "rb") as f:
    v3_sha = hashlib.sha256(f.read()).hexdigest()
print(f"[PRE-REGISTRATION] thresholds_v3.json loaded, SHA256: {v3_sha}")

# ═══════════════════════════════════════════════════════════════════════════
# 1. HELD-OUT VALIDATION SET (Strictly excluding FAST/Thermo/Dura/Hot backgrounds)
# ═══════════════════════════════════════════════════════════════════════════
# Rule: Exclude positions/mutations in FAST (121, 186, 224, 233, 280),
# Dura (117, 119, 140, 159, 165, 168, 180, 188, 214, 280), Thermo (121, 186, 280).
# All data from primary peer-reviewed papers with exact table/figure citations.

HELD_OUT_LITERATURE_VARIANTS = [
    {
        "variant_id": "Son2019_P181A",
        "mutations_6eqe": [{"wt": "P", "pos": 181, "mut": "A"}],
        "delta_Tm_exp_C": 0.45,
        "Tm_exp_C": 49.25,
        "method": "DSF (Protein Thermal Shift)",
        "doi": "10.1021/acscatal.9b00568",
        "source": "Son et al. 2019, ACS Catalysis, Table S2",
        "category": "held_out"
    },
    {
        "variant_id": "Son2019_W185A",
        "mutations_6eqe": [{"wt": "W", "pos": 185, "mut": "A"}],
        "delta_Tm_exp_C": -2.80,
        "Tm_exp_C": 46.00,
        "method": "DSF (Protein Thermal Shift)",
        "doi": "10.1021/acscatal.9b00568",
        "source": "Son et al. 2019, ACS Catalysis, Table S2",
        "category": "held_out"
    },
    {
        "variant_id": "Son2019_S160A",
        "mutations_6eqe": [{"wt": "S", "pos": 160, "mut": "A"}],
        "delta_Tm_exp_C": -4.50,
        "Tm_exp_C": 44.30,
        "method": "DSF (Protein Thermal Shift)",
        "doi": "10.1021/acscatal.9b00568",
        "source": "Son et al. 2019, ACS Catalysis, Table S2",
        "category": "held_out_control"
    },
    {
        "variant_id": "Joo2018_C203A_C239A",
        "mutations_6eqe": [
            {"wt": "C", "pos": 203, "mut": "A"},
            {"wt": "C", "pos": 239, "mut": "A"}
        ],
        "delta_Tm_exp_C": -13.00,
        "Tm_exp_C": 35.80,
        "method": "Protein Thermal Shift (DSF)",
        "doi": "10.1038/s41467-018-02881-1",
        "source": "Joo et al. 2018, Nature Communications, Fig. 2c",
        "category": "held_out_disulfide_knockout"
    },
    {
        "variant_id": "Meng2021_F229Y",
        "mutations_6eqe": [{"wt": "F", "pos": 229, "mut": "Y"}],
        "delta_Tm_exp_C": 3.60,
        "Tm_exp_C": 54.40,
        "method": "DSF (SYPRO Orange)",
        "doi": "10.1016/j.ijbiomac.2021.03.058",
        "source": "Meng et al. 2021, Int. J. Biol. Macromol., Table 1",
        "category": "held_out"
    },
    {
        "variant_id": "Ding2024_K95A",
        "mutations_6eqe": [{"wt": "K", "pos": 95, "mut": "A"}],
        "delta_Tm_exp_C": 1.20,
        "Tm_exp_C": 44.70,
        "method": "DSF (Prometheus nanoDSF)",
        "doi": "10.1016/j.bej.2024.109420",
        "source": "Ding et al. 2024, Biochem. Eng. J., Table 1",
        "category": "held_out"
    },
    {
        "variant_id": "Ding2024_R132N",
        "mutations_6eqe": [{"wt": "R", "pos": 132, "mut": "N"}],
        "delta_Tm_exp_C": -0.80,
        "Tm_exp_C": 42.70,
        "method": "DSF (Prometheus nanoDSF)",
        "doi": "10.1016/j.bej.2024.109420",
        "source": "Ding et al. 2024, Biochem. Eng. J., Table 1",
        "category": "held_out"
    },
    {
        "variant_id": "Brott2022_F201I",
        "mutations_6eqe": [{"wt": "F", "pos": 201, "mut": "I"}],
        "delta_Tm_exp_C": 2.80,
        "Tm_exp_C": 51.60,
        "method": "DSF (nanoDSF)",
        "doi": "10.1002/elsc.202100105",
        "source": "Brott et al. 2022, Eng. Life Sci., Table 1",
        "category": "held_out"
    }
]

# Benchmark set including background mutations (for comparison)
BACKGROUND_BENCHMARK_VARIANTS = [
    {
        "variant_id": "Son2019_S121E",
        "mutations_6eqe": [{"wt": "S", "pos": 121, "mut": "E"}],
        "delta_Tm_exp_C": 2.10,
        "doi": "10.1021/acscatal.9b00568",
        "source": "Son et al. 2019, ACS Catalysis, Table S2"
    },
    {
        "variant_id": "Son2019_D186H",
        "mutations_6eqe": [{"wt": "D", "pos": 186, "mut": "H"}],
        "delta_Tm_exp_C": 3.40,
        "doi": "10.1021/acscatal.9b00568",
        "source": "Son et al. 2019, ACS Catalysis, Table S2"
    },
    {
        "variant_id": "Son2019_R280A",
        "mutations_6eqe": [{"wt": "R", "pos": 280, "mut": "A"}],
        "delta_Tm_exp_C": 1.90,
        "doi": "10.1021/acscatal.9b00568",
        "source": "Son et al. 2019, ACS Catalysis, Table S2"
    },
    {
        "variant_id": "Lu2022_N233K",
        "mutations_6eqe": [{"wt": "N", "pos": 233, "mut": "K"}],
        "delta_Tm_exp_C": 3.20,
        "doi": "10.1038/s41586-022-04599-z",
        "source": "Lu et al. 2022, Nature, Extended Data Fig. 3"
    },
    {
        "variant_id": "Lu2022_R224Q",
        "mutations_6eqe": [{"wt": "R", "pos": 224, "mut": "Q"}],
        "delta_Tm_exp_C": 1.40,
        "doi": "10.1038/s41586-022-04599-z",
        "source": "Lu et al. 2022, Nature, Extended Data Fig. 3"
    }
]

# ═══════════════════════════════════════════════════════════════════════════
# 2. MODEL LOADERS
# ═══════════════════════════════════════════════════════════════════════════

def load_thermompnn(device):
    from omegaconf import OmegaConf
    from transfer_model import TransferModel

    cfg = OmegaConf.load(os.path.join(THERMOMPNN_DIR, "local.yaml"))
    model_cfg = OmegaConf.create({
        "training": {"num_workers": 4, "learn_rate": 0.001, "epochs": 100, "lr_schedule": True, "mpnn_learn_rate": 0.0001},
        "model": {"hidden_dims": [64, 32], "subtract_mut": True, "num_final_layers": 2, "freeze_weights": True, "load_pretrained": True, "lightattn": True},
    })
    merged = OmegaConf.merge(cfg, model_cfg)
    merged.platform.thermompnn_dir = THERMOMPNN_DIR

    ckpt_path = os.path.join(THERMOMPNN_DIR, "models", "thermoMPNN_default.pt")
    ckpt = torch.load(ckpt_path, map_location="cpu")
    raw_sd = ckpt["state_dict"]
    model_sd = {k[len("model."):]: v for k, v in raw_sd.items() if k.startswith("model.")}
    model = TransferModel(merged)
    model.load_state_dict(model_sd, strict=False)
    return model.eval().to(device)


def load_esm2(device):
    from transformers import AutoTokenizer, EsmForMaskedLM
    model_name = "facebook/esm2_t6_8M_UR50D"
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = EsmForMaskedLM.from_pretrained(model_name)
    return tokenizer, model.eval().to(device)


# ═══════════════════════════════════════════════════════════════════════════
# 3. PREDICTOR SCORING FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════

def score_thermompnn_variant(model, mut_pdb, mutations_list, idx_map, device):
    from datasets import Mutation as ThermoMutation
    mut_objs = []
    for m in mutations_list:
        pos_seq = idx_map[m["pos"]]
        mut_objs.append(ThermoMutation(
            position=pos_seq,
            wildtype=m["wt"],
            mutation=m["mut"],
            ddG=None,
            pdb=mut_pdb[0]["name"]
        ))
    with torch.no_grad():
        preds, _ = model(mut_pdb, mut_objs)
    # Return sum of predicted ddGs (additive approximation)
    return float(sum(p["ddG"].cpu().item() for p in preds))


def score_esm2_variant(tokenizer, model, wt_sequence_290, mutations_list, device):
    """
    Log-Likelihood Ratio (LLR): log P(mut | context) - log P(wt | context)
    Positive LLR = favorable in evolutionary context.
    Negative LLR = deleterious.
    """
    seq_chars = list(wt_sequence_290)
    total_llr = 0.0

    for m in mutations_list:
        pos_1based = m["pos"]  # 1-based in UniProt / 6EQE
        pos_0based = pos_1based - 1

        # Check WT
        assert seq_chars[pos_0based] == m["wt"], f"WT mismatch at {pos_1based}: expected {m['wt']}, got {seq_chars[pos_0based]}"

        # Masked sequence
        masked_chars = list(seq_chars)
        masked_chars[pos_0based] = tokenizer.mask_token
        masked_str = "".join(masked_chars)

        tokens = tokenizer(masked_str, return_tensors="pt").to(device)
        with torch.no_grad():
            logits = model(**tokens).logits

        mask_token_idx = (tokens.input_ids == tokenizer.mask_token_id).nonzero(as_tuple=True)[1][0]
        probs = torch.softmax(logits[0, mask_token_idx], dim=-1)

        wt_tok = tokenizer.convert_tokens_to_ids(m["wt"])
        mut_tok = tokenizer.convert_tokens_to_ids(m["mut"])

        wt_prob = probs[wt_tok].item()
        mut_prob = probs[mut_tok].item()

        llr = np.log(mut_prob + 1e-12) - np.log(wt_prob + 1e-12)
        total_llr += llr

    return float(total_llr)


# ═══════════════════════════════════════════════════════════════════════════
# 4. BOOTSTRAP STATISTICAL ENGINE
# ═══════════════════════════════════════════════════════════════════════════

def bootstrap_spearman(x, y, n_resamples=10000, alpha=0.05, seed=42):
    rng = np.random.default_rng(seed)
    n = len(x)
    observed_rho, _ = spearmanr(x, y)

    boot_rhos = []
    for _ in range(n_resamples):
        indices = rng.choice(n, size=n, replace=True)
        # Avoid zero variance
        if len(np.unique(x[indices])) < 2 or len(np.unique(y[indices])) < 2:
            continue
        r, _ = spearmanr(x[indices], y[indices])
        if not np.isnan(r):
            boot_rhos.append(r)

    boot_rhos = np.array(boot_rhos)
    low_ci = float(np.percentile(boot_rhos, (alpha / 2.0) * 100))
    high_ci = float(np.percentile(boot_rhos, (1.0 - alpha / 2.0) * 100))

    return {
        "observed_rho": float(observed_rho),
        "ci_lower_95": low_ci,
        "ci_upper_95": high_ci,
        "ci_excludes_zero": bool(low_ci > 0.0 or high_ci < 0.0),
        "expected_negative_sign": bool(observed_rho < 0.0),
        "n_samples": n,
        "n_valid_bootstraps": len(boot_rhos)
    }


# ═══════════════════════════════════════════════════════════════════════════
# MAIN EXECUTION
# ═══════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("PHASE C: HELD-OUT VALIDATION & SKEPTICAL AUDIT")
    print("=" * 70)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Compute device: {device}")

    # Load 6EQE structure
    from Bio.PDB import PDBParser
    from protein_mpnn_utils import alt_parse_PDB

    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("6EQE", PDB_6EQE)
    chain = list(struct[0].get_chains())[0]
    res_nums = [r.id[1] for r in chain if r.id[0] == " "]
    idx_map = {rn: i for i, rn in enumerate(res_nums)}

    mut_pdb = alt_parse_PDB(PDB_6EQE, chain.id)
    wt_seq_structure = mut_pdb[0]["seq"]

    # Full WT sequence (290 aa, UniProt A0A0K8P6T7)
    wt_seq_290 = (
        "MNFPRASRLMQAAVLGGLMAVSAAATAQTNPYARGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPR"
        "LASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLI"
        "FACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCS"
    )

    # 1. Load models
    print("\n[STEP 1] Loading models...")
    thermompnn = load_thermompnn(device)
    print("  ThermoMPNN v1: LOADED")
    esm_tokenizer, esm_model = load_esm2(device)
    print("  ESM-2 8M: LOADED")

    # 2. Score Held-out Validation Set
    print(f"\n[STEP 2] Scoring Held-out Validation Set (n={len(HELD_OUT_LITERATURE_VARIANTS)})...")
    held_out_results = []

    for v in HELD_OUT_LITERATURE_VARIANTS:
        # Check mapping to structure
        in_struct = all(m["pos"] in idx_map for m in v["mutations_6eqe"])
        if not in_struct:
            print(f"  [WARN] {v['variant_id']} not in 6EQE structure!")
            continue

        pred_ddg = score_thermompnn_variant(thermompnn, mut_pdb, v["mutations_6eqe"], idx_map, device)
        pred_llr = score_esm2_variant(esm_tokenizer, esm_model, wt_seq_290, v["mutations_6eqe"], device)

        row = dict(v)
        row["ddG_pred_thermompnn_kcal_mol"] = round(pred_ddg, 4)
        row["llr_pred_esm2"] = round(pred_llr, 4)
        held_out_results.append(row)

        print(f"  {v['variant_id']:25s} | exp ΔTm: {v['delta_Tm_exp_C']:+6.2f}°C | ThermoMPNN ddG: {pred_ddg:+6.3f} kcal/mol | ESM-2 LLR: {pred_llr:+6.3f}")

    # 3. Statistical Analysis on Held-out Set
    exp_dtms = np.array([r["delta_Tm_exp_C"] for r in held_out_results])
    pred_ddgs = np.array([r["ddG_pred_thermompnn_kcal_mol"] for r in held_out_results])
    pred_llrs = np.array([r["llr_pred_esm2"] for r in held_out_results])

    print("\n[STEP 3] Running 10,000-replicate bootstrap on held-out validation set...")
    stats_thermompnn = bootstrap_spearman(pred_ddgs, exp_dtms, n_resamples=10000, seed=42)
    stats_esm2 = bootstrap_spearman(pred_llrs, exp_dtms, n_resamples=10000, seed=42)

    print(f"  ThermoMPNN vs Exp ΔTm:")
    print(f"    Observed Spearman rho: {stats_thermompnn['observed_rho']:+.4f}")
    print(f"    95% Bootstrap CI:     [{stats_thermompnn['ci_lower_95']:+.4f}, {stats_thermompnn['ci_upper_95']:+.4f}]")
    print(f"    CI excludes zero:      {stats_thermompnn['ci_excludes_zero']}")

    print(f"  ESM-2 LLR vs Exp ΔTm:")
    print(f"    Observed Spearman rho: {stats_esm2['observed_rho']:+.4f} (expected positive)")
    print(f"    95% Bootstrap CI:     [{stats_esm2['ci_lower_95']:+.4f}, {stats_esm2['ci_upper_95']:+.4f}]")
    print(f"    CI excludes zero:      {stats_esm2['ci_excludes_zero']}")

    # 4. Evaluate Pre-registered Pass Criterion
    # Pre-registered pass: |rho| >= 0.50 AND CI excludes 0 AND expected sign (rho < 0 for ThermoMPNN)
    pass_thermompnn = (
        stats_thermompnn["observed_rho"] <= -0.50 and
        stats_thermompnn["ci_upper_95"] < 0.0
    )

    print(f"\n[STEP 4] Pre-Registered Decision Gate Evaluation:")
    print(f"  Requirement 1: |rho| >= 0.50 -> Observed |rho| = {abs(stats_thermompnn['observed_rho']):.4f} ({'PASS' if abs(stats_thermompnn['observed_rho']) >= 0.50 else 'FAIL'})")
    print(f"  Requirement 2: Expected negative sign -> {'PASS' if stats_thermompnn['expected_negative_sign'] else 'FAIL'}")
    print(f"  Requirement 3: 95% CI excludes zero -> {'PASS' if stats_thermompnn['ci_excludes_zero'] else 'FAIL (CI spans zero)'}")

    # DECISION RULE
    if pass_thermompnn:
        project_outcome = "DESIGN"
        outcome_reason = "Retrospective validation passed pre-registered gates on held-out literature set."
    else:
        project_outcome = "AUDIT-ONLY"
        outcome_reason = (
            f"Retrospective validation failed pre-registered gate on held-out set "
            f"(rho = {stats_thermompnn['observed_rho']:+.3f}, 95% CI [{stats_thermompnn['ci_lower_95']:+.3f}, {stats_thermompnn['ci_upper_95']:+.3f}], "
            f"n={len(held_out_results)}). Because sample size is small and confidence interval spans zero / does not satisfy pre-registered thresholds, "
            f"claims of computationally validated design are NOT supported. Manuscript OUTCOME is locked to AUDIT-ONLY."
        )

    print(f"\n>>> PRE-REGISTERED PROJECT OUTCOME: {project_outcome} <<<")
    print(f"Reason: {outcome_reason}")

    # Save full results
    output_data = {
        "timestamp_utc": datetime.datetime.utcnow().isoformat() + "Z",
        "thresholds_v3_sha256": v3_sha,
        "actual_n_held_out": len(held_out_results),
        "held_out_variants": held_out_results,
        "statistical_evaluation": {
            "thermompnn": stats_thermompnn,
            "esm2": stats_esm2
        },
        "pre_registered_pass_criteria": {
            "threshold_min_abs_rho": 0.50,
            "ci_must_exclude_zero": True,
            "expected_sign": "negative",
            "passed": pass_thermompnn
        },
        "project_outcome": project_outcome,
        "outcome_reason": outcome_reason
    }

    out_file = os.path.join(WORK_DIR, "phase_c_validation_results.json")
    with open(out_file, "w") as f:
        json.dump(output_data, f, indent=2)
    print(f"\nResults written to: {out_file}")

    return output_data


if __name__ == "__main__":
    main()
