#!/usr/bin/env python3
"""
Generate publication-quality figures for the forensic audit manuscript.
Outputs both PNG (300 DPI) and PDF vector formats.
"""

import os
import json
import numpy as np
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

WORK_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(WORK_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Set style
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.titlesize": 13,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight"
})

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 1: Circular Benchmark Leakage vs. Held-Out Generalization
# ═══════════════════════════════════════════════════════════════════════════

def generate_figure1():
    # Load Phase C data
    with open(os.path.join(WORK_DIR, "phase_c_validation_results.json")) as f:
        data = json.load(f)

    held_out = data["held_out_variants"]
    exp_dtm_held = [v["delta_Tm_exp_C"] for v in held_out]
    pred_ddg_held = [v["ddG_pred_thermompnn_kcal_mol"] for v in held_out]
    labels_held = [v["variant_id"].replace("Son2019_", "").replace("Ding2024_", "").replace("Meng2021_", "").replace("Joo2018_", "").replace("Brott2022_", "") for v in held_out]

    # Benchmark variants (seen in backgrounds)
    bench_variants = [
        ("S121E", 2.10, 0.0216),
        ("D186H", 3.40, -0.5120),
        ("R280A", 1.90, -0.0686),
        ("N233K", 3.20, -1.1984),
        ("R224Q", 1.40, 0.3777),
        ("S160A", -4.50, 0.4518),
        ("W185A", -2.80, 2.1512),
    ]
    labels_bench = [b[0] for b in bench_variants]
    exp_dtm_bench = [b[1] for b in bench_variants]
    pred_ddg_bench = [b[2] for b in bench_variants]

    fig = plt.figure(figsize=(10, 4.2))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1], wspace=0.3)

    # Panel A: Benchmark Set (Seen in backgrounds)
    ax1 = fig.add_subplot(gs[0])
    ax1.scatter(exp_dtm_bench, pred_ddg_bench, color="#1f77b4", s=60, edgecolors="k", zorder=3)
    for i, txt in enumerate(labels_bench):
        offset = (5, 5) if i % 2 == 0 else (-25, -15)
        ax1.annotate(txt, (exp_dtm_bench[i], pred_ddg_bench[i]), textcoords="offset points", xytext=offset, fontsize=8)

    # Trendline
    m1, b1 = np.polyfit(exp_dtm_bench, pred_ddg_bench, 1)
    x_vals1 = np.linspace(min(exp_dtm_bench)-1, max(exp_dtm_bench)+1, 100)
    ax1.plot(x_vals1, m1*x_vals1 + b1, "--", color="#1f77b4", alpha=0.7, label="Linear fit")

    ax1.axhline(0, color="gray", linestyle=":", alpha=0.6)
    ax1.axvline(0, color="gray", linestyle=":", alpha=0.6)
    ax1.set_xlabel("Experimental $\\Delta T_m$ (°C)")
    ax1.set_ylabel("ThermoMPNN Predicted $\\Delta\\Delta G$ (kcal/mol)")
    ax1.set_title("A. Benchmark Variants (Seen Backgrounds)\nSpearman $\\rho = -0.893$ ($p = 0.007$)", fontweight="bold")
    ax1.grid(True, linestyle="--", alpha=0.4)

    # Panel B: Held-Out Variants (Strict Exclusion)
    ax2 = fig.add_subplot(gs[1])
    ax2.scatter(exp_dtm_held, pred_ddg_held, color="#d62728", s=60, edgecolors="k", zorder=3)
    for i, txt in enumerate(labels_held):
        offset = (5, 5) if i % 2 == 0 else (-35, -12)
        ax2.annotate(txt, (exp_dtm_held[i], pred_ddg_held[i]), textcoords="offset points", xytext=offset, fontsize=8)

    # Trendline
    m2, b2 = np.polyfit(exp_dtm_held, pred_ddg_held, 1)
    x_vals2 = np.linspace(min(exp_dtm_held)-1, max(exp_dtm_held)+1, 100)
    ax2.plot(x_vals2, m2*x_vals2 + b2, "--", color="#d62728", alpha=0.7, label="Linear fit")

    ax2.axhline(0, color="gray", linestyle=":", alpha=0.6)
    ax2.axvline(0, color="gray", linestyle=":", alpha=0.6)
    ax2.set_xlabel("Experimental $\\Delta T_m$ (°C)")
    ax2.set_ylabel("ThermoMPNN Predicted $\\Delta\\Delta G$ (kcal/mol)")
    ax2.set_title("B. Held-Out Literature Variants ($n=8$)\nSpearman $\\rho = +0.024$ (95% CI: [-0.68, +0.92])", fontweight="bold")
    ax2.grid(True, linestyle="--", alpha=0.4)

    fig.tight_layout()
    png_path = os.path.join(FIG_DIR, "figure1_validation_comparison.png")
    pdf_path = os.path.join(FIG_DIR, "figure1_validation_comparison.pdf")
    fig.savefig(png_path)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"Generated Figure 1: {png_path} and {pdf_path}")


# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 2: Bootstrap Distribution of Correlation Coefficients
# ═══════════════════════════════════════════════════════════════════════════

def generate_figure2():
    with open(os.path.join(WORK_DIR, "phase_c_validation_results.json")) as f:
        data = json.load(f)

    held_out = data["held_out_variants"]
    exp_dtm = np.array([v["delta_Tm_exp_C"] for v in held_out])
    pred_ddg = np.array([v["ddG_pred_thermompnn_kcal_mol"] for v in held_out])
    pred_llr = np.array([v["llr_pred_esm2"] for v in held_out])

    rng = np.random.default_rng(42)
    n = len(exp_dtm)
    boot_rhos_ddg = []
    boot_rhos_llr = []

    for _ in range(10000):
        idx = rng.choice(n, size=n, replace=True)
        if len(np.unique(pred_ddg[idx])) > 1 and len(np.unique(exp_dtm[idx])) > 1:
            r1, _ = spearmanr(pred_ddg[idx], exp_dtm[idx])
            if not np.isnan(r1): boot_rhos_ddg.append(r1)
        if len(np.unique(pred_llr[idx])) > 1 and len(np.unique(exp_dtm[idx])) > 1:
            r2, _ = spearmanr(pred_llr[idx], exp_dtm[idx])
            if not np.isnan(r2): boot_rhos_llr.append(r2)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8))

    # ThermoMPNN
    ax1.hist(boot_rhos_ddg, bins=40, color="#2ca02c", alpha=0.7, edgecolor="k", density=True)
    ax1.axvline(0, color="k", linestyle="--", linewidth=1.5, label="$\\rho = 0$")
    ax1.axvline(-0.50, color="r", linestyle=":", linewidth=2, label="Pre-registered threshold ($\\rho \\leq -0.50$)")
    ax1.set_xlabel("Bootstrap Spearman $\\rho$")
    ax1.set_ylabel("Density")
    ax1.set_title("A. ThermoMPNN Bootstrap Distribution\n(10,000 resamples, 95% CI: [-0.68, +0.92])", fontweight="bold")
    ax1.legend(loc="upper right", fontsize=8)
    ax1.grid(True, linestyle="--", alpha=0.3)

    # ESM-2 LLR
    ax2.hist(boot_rhos_llr, bins=40, color="#9467bd", alpha=0.7, edgecolor="k", density=True)
    ax2.axvline(0, color="k", linestyle="--", linewidth=1.5, label="$\\rho = 0$")
    ax2.set_xlabel("Bootstrap Spearman $\\rho$")
    ax2.set_ylabel("Density")
    ax2.set_title("B. ESM-2 LLR Bootstrap Distribution\n(10,000 resamples, 95% CI: [-0.85, +0.90])", fontweight="bold")
    ax2.legend(loc="upper right", fontsize=8)
    ax2.grid(True, linestyle="--", alpha=0.3)

    fig.tight_layout()
    png_path = os.path.join(FIG_DIR, "figure2_bootstrap_distributions.png")
    pdf_path = os.path.join(FIG_DIR, "figure2_bootstrap_distributions.pdf")
    fig.savefig(png_path)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"Generated Figure 2: {png_path} and {pdf_path}")


if __name__ == "__main__":
    generate_figure1()
    generate_figure2()
