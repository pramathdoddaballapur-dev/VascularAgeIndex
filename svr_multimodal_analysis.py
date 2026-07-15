"""
svr_multimodal_analysis.py
--------------------------
Full multimodal biological age clock — Manning Lab
Integrates 2-photon collagen architecture with PA vascular mechanics.

Feature set (primary multimodal model, manuscript v8):
    vm_circ   = von Mises PDF at 90°  (circumferential orientation)
    vm_axial  = von Mises PDF at 0°   (axial orientation)
    COLLAGEN STRAIGHTNESS
    Distensibility (mmHg-1)
    PA Circumferential Stiffness (MPa)
    PWV  (m/s)

Performance (normotensive aging cohort, 6–24mo):
    n=36  |  LOO R²=0.834  |  MAE=2.26 mo
    Female-only: n=18, R²=0.960  |  Male-only: n=18, R²=0.658

Preprocessing per LOO fold:
    1. Winsorise at 5/95th percentile (fit on training fold only)
    2. Yeo-Johnson power transform
    3. StandardScaler
    → no data leakage

SVR hyperparameters: kernel=rbf, C=200, ε=1.5, γ=0.05

Figures produced:
    Figure 1 – Stepwise model comparison (2-photon only → add distensibility
               → add circ. stiffness → add PWV)
    Figure 2 – Predicted vs actual age (LOO-CV), coloured by sex
    Figure 3 – Residuals by age group (box + strip)
    Figure 4 – Sex-stratified analysis (per-sex LOO + R² comparison)
    Figure 5 – PCHIP parameter trajectories (θ, κ, straightness, distensibility)

Usage:
    python svr_multimodal_analysis.py
    python svr_multimodal_analysis.py Aging_PA_merged_updated.xlsx

Column names expected in Excel (must match exactly):
    Sample Name, Age, Sex
    ORIENTATION-θ, ORIENTATION-κ, COLLAGEN STRAIGHTNESS
    Distensibility (mmHg-1), PA Circumferential Stiffness (MPa), PWV
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.special import i0 as bessel_i0
from scipy.interpolate import PchipInterpolator
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler, PowerTransformer
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import r2_score
from sklearn.model_selection import LeaveOneOut


# ── Colours ───────────────────────────────────────────────────────────────────
C_MALE    = "#185FA5"
C_FEMALE  = "#993556"
C_2PHOTON = "#0F6E56"   # 2-photon-only baseline (green)
C_MULTI   = "#993556"   # multimodal (pink)
C_LINEAR  = "#185FA5"
C_GBM     = "#BA7517"
GREY      = "#888780"

# ── Constants ─────────────────────────────────────────────────────────────────
AGE_MAP = {"3mo": 6, "6mo": 6, "12mo": 12, "18mo": 18, "24mo": 24}
REF_CIRC  = np.deg2rad(90)
REF_AXIAL = np.deg2rad(0)

# SVR hyperparameters (tuned on normotensive cohort — do not retune on new data)
BEST_C     = 200
BEST_EPS   = 1.5
BEST_GAMMA = 0.05

# Feature sets
PHOTON_FEATURES = ["vm_circ", "vm_axial", "COLLAGEN STRAIGHTNESS"]
MECH_FEATURES   = ["Distensibility (mmHg-1)",
                    "PA Circumferential Stiffness (MPa)",
                    "PWV"]
MULTIMODAL_FEATURES = PHOTON_FEATURES + MECH_FEATURES   # primary model


# ── Von Mises PDF ─────────────────────────────────────────────────────────────
def vm_pdf(x_rad, theta_deg, kappa):
    theta_rad = np.deg2rad(theta_deg)
    return np.exp(kappa * np.cos(x_rad - theta_rad)) / (2 * np.pi * bessel_i0(kappa))


# ── Data loading ──────────────────────────────────────────────────────────────
def load_data(path: str) -> pd.DataFrame:
    df = pd.read_excel(path)
    df = df[df["Sample Name"].notna() & (df["Sample Name"] != "Average")]
    df = df[df["Age"].notna()]
    df["Age_num"] = df["Age"].map(AGE_MAP)
    df = df[df["Age_num"].notna()].copy()   # drops 3mo and 27mo

    # Compute von Mises features
    theta = df["ORIENTATION-θ"]
    kappa = df["ORIENTATION-κ"]
    valid = theta.notna() & kappa.notna()
    df.loc[valid, "vm_circ"]  = vm_pdf(REF_CIRC,  theta[valid], kappa[valid])
    df.loc[valid, "vm_axial"] = vm_pdf(REF_AXIAL, theta[valid], kappa[valid])

    return df


# ── LOO helpers ───────────────────────────────────────────────────────────────
def loo_predict_svr(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    LOO predictions with per-fold preprocessing:
      Winsorise (5/95th pct) → Yeo-Johnson → StandardScaler
    All transformers fit on training fold only — no leakage.
    """
    preds = np.empty(len(y))
    for train_idx, test_idx in LeaveOneOut().split(X):
        X_tr, X_te = X[train_idx], X[test_idx]
        # Winsorise
        lo = np.percentile(X_tr, 5,  axis=0)
        hi = np.percentile(X_tr, 95, axis=0)
        X_tr = np.clip(X_tr, lo, hi)
        X_te = np.clip(X_te, lo, hi)
        # Yeo-Johnson
        pt = PowerTransformer(method='yeo-johnson')
        X_tr = pt.fit_transform(X_tr)
        X_te = pt.transform(X_te)
        # StandardScaler
        sc = StandardScaler()
        X_tr = sc.fit_transform(X_tr)
        X_te = sc.transform(X_te)
        # Fit and predict
        m = SVR(kernel="rbf", C=BEST_C, epsilon=BEST_EPS, gamma=BEST_GAMMA)
        m.fit(X_tr, y[train_idx])
        preds[test_idx] = m.predict(X_te)
    return preds


def loo_r2_mae(preds, y):
    ss_res = np.sum((y - preds) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2  = 1 - ss_res / ss_tot
    mae = np.mean(np.abs(y - preds))
    return r2, mae


# ── Figure 1: stepwise model comparison ──────────────────────────────────────
def plot_stepwise_comparison(df: pd.DataFrame, ax: plt.Axes):
    """
    Shows LOO R² as mechanical features are progressively added to the
    2-photon baseline. Mirrors Figure 3C logic from the manuscript.
    """
    configs = [
        ("2-photon only\n(vm_circ + vm_axial\n+ Collagen)", PHOTON_FEATURES),
        ("+ Distensibility",
         PHOTON_FEATURES + ["Distensibility (mmHg-1)"]),
        ("+ Circ. Stiffness",
         PHOTON_FEATURES + ["Distensibility (mmHg-1)",
                             "PA Circumferential Stiffness (MPa)"]),
        ("+ PWV\n[Primary multimodal]",
         MULTIMODAL_FEATURES),
    ]

    r2_vals, mae_vals, n_vals = [], [], []
    for _, feats in configs:
        sub = df[["Age_num"] + feats].dropna()
        X   = sub[feats].values
        y   = sub["Age_num"].values
        p   = loo_predict_svr(X, y)
        r2, mae = loo_r2_mae(p, y)
        r2_vals.append(r2)
        mae_vals.append(mae)
        n_vals.append(len(y))

    labels = [c[0] for c in configs]
    colors = [C_2PHOTON, "#5A9E8C", "#2E7D9A", C_MULTI]
    x = np.arange(len(configs))

    bars = ax.bar(x, r2_vals, color=colors, alpha=0.88, edgecolor="white",
                  zorder=3, width=0.55)
    for bar, r2v, maev, nv in zip(bars, r2_vals, mae_vals, n_vals):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.01,
                f"R²={r2v:.3f}\nMAE={maev:.2f}mo\nn={nv}",
                ha="center", va="bottom", fontsize=8, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("LOO-CV R²", fontsize=10)
    ax.set_ylim(0, 1.0)
    ax.set_title(
        "Stepwise model improvement — 2-photon → multimodal\n"
        f"SVR RBF  |  C={BEST_C}, ε={BEST_EPS}, γ={BEST_GAMMA}  |  "
        "Winsorise 5/95 + Yeo-Johnson + StandardScaler",
        fontsize=10, fontweight="bold"
    )
    ax.axhline(0, color=GREY, lw=0.8)
    ax.grid(axis="y", color="#dddddd", zorder=0)

    # Annotate delta R² arrows
    for i in range(1, len(r2_vals)):
        delta = r2_vals[i] - r2_vals[i - 1]
        sign  = "+" if delta >= 0 else ""
        ax.annotate(
            f"{sign}{delta:.3f}",
            xy=(i, r2_vals[i] + 0.04),
            ha="center", fontsize=8, color="#333333",
            fontweight="bold"
        )


# ── Figure 2: predicted vs actual (multimodal) ───────────────────────────────
def plot_predicted_vs_actual(df: pd.DataFrame, ax: plt.Axes):
    sub  = df[["Age_num", "Sex"] + MULTIMODAL_FEATURES].dropna()
    X    = sub[MULTIMODAL_FEATURES].values
    y    = sub["Age_num"].values
    sex  = sub["Sex"].values

    preds     = loo_predict_svr(X, y)
    r2, mae   = loo_r2_mae(preds, y)

    for label, color in [("Male", C_MALE), ("Female", C_FEMALE)]:
        mask = sex == label
        ax.scatter(y[mask], preds[mask], c=color, alpha=0.72, s=58,
                   edgecolors="white", linewidths=0.4, label=label, zorder=3)

    lims = [4, 26]
    ax.plot(lims, lims, ls="--", lw=1.5, color=GREY, zorder=2, label="Identity")
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel("Actual age (months)", fontsize=9)
    ax.set_ylabel("Predicted age (months)", fontsize=9)
    ax.set_title(
        f"Multimodal SVR — predicted vs actual age  (LOO-CV)\n"
        f"vm_circ + vm_axial + Collagen + Distensibility + PA Circ. Stiffness + PWV\n"
        f"LOO R² = {r2:.3f}  |  MAE = {mae:.2f} mo  |  n = {len(y)}",
        fontsize=9, fontweight="bold"
    )
    ax.legend(fontsize=8)
    ax.grid(color="#eeeeee", zorder=0)
    ax.set_xticks([6, 12, 18, 24])
    return preds, y, sex


# ── Figure 3: residuals ───────────────────────────────────────────────────────
def plot_residuals(preds, y, sex, ax: plt.Axes):
    residuals  = preds - y
    age_groups = sorted(np.unique(y).astype(int))
    rng        = np.random.default_rng(42)

    for gi, age in enumerate(age_groups):
        mask = y == age
        res  = residuals[mask]
        sx   = sex[mask]

        q1, med, q3 = np.percentile(res, [25, 50, 75])
        bw = 0.22
        ax.broken_barh([(gi - bw, 2 * bw)], (q1, q3 - q1),
                       facecolors="none", edgecolors=C_LINEAR,
                       linewidth=1.2, zorder=3)
        ax.hlines(med, gi - bw, gi + bw, colors=C_LINEAR, linewidth=2, zorder=4)
        ax.vlines(gi, res.min(), q1,     colors=C_LINEAR, linewidth=0.9, zorder=3)
        ax.vlines(gi, q3, res.max(),     colors=C_LINEAR, linewidth=0.9, zorder=3)
        ax.hlines([res.min(), res.max()],
                  gi - 0.1, gi + 0.1,   colors=C_LINEAR, linewidth=0.9, zorder=3)

        for label, color in [("Male", C_MALE), ("Female", C_FEMALE)]:
            pts = res[sx == label]
            jx  = gi + rng.uniform(-0.28, 0.28, size=len(pts))
            ax.scatter(jx, pts, c=color, alpha=0.6, s=28,
                       edgecolors="white", linewidths=0.3, zorder=5)

    ax.axhline(0, ls="--", lw=1.5, color=GREY, zorder=2)
    ax.set_xticks(range(len(age_groups)))
    ax.set_xticklabels([f"{g}mo" for g in age_groups], fontsize=9)
    ax.set_xlabel("Actual age group", fontsize=9)
    ax.set_ylabel("Residual  (predicted − actual, months)", fontsize=9)
    ax.set_title("Residuals by age group — Multimodal SVR (LOO-CV)",
                 fontsize=10, fontweight="bold")
    patches = [mpatches.Patch(color=C_MALE, label="Male"),
               mpatches.Patch(color=C_FEMALE, label="Female")]
    ax.legend(handles=patches, fontsize=8)
    ax.grid(axis="y", color="#eeeeee", zorder=0)


# ── Figure 4: sex-stratified ──────────────────────────────────────────────────
def plot_sex_stratified(df: pd.DataFrame,
                        ax_male, ax_female, ax_r2):
    sub = df[["Age_num", "Sex"] + MULTIMODAL_FEATURES].dropna()

    sex_data = {}
    for sex_label, color, ax in [
        ("Male",   C_MALE,   ax_male),
        ("Female", C_FEMALE, ax_female)
    ]:
        s = sub[sub["Sex"] == sex_label]
        if len(s) < 5:
            ax.set_title(f"{sex_label}: insufficient data", fontsize=10)
            continue
        X = s[MULTIMODAL_FEATURES].values
        y = s["Age_num"].values

        preds   = loo_predict_svr(X, y)
        r2, mae = loo_r2_mae(preds, y)
        sex_data[sex_label] = {"y": y, "preds": preds, "r2": r2, "mae": mae}

        lims = [4, 26]
        ax.scatter(y, preds, c=color, alpha=0.75, s=60,
                   edgecolors="white", linewidths=0.5, zorder=3)
        ax.plot(lims, lims, ls="--", lw=1.5, color=GREY, zorder=2)
        ax.set_xlim(lims); ax.set_ylim(lims)
        ax.set_xticks([6, 12, 18, 24])
        ax.set_xlabel("Actual age (months)", fontsize=9)
        ax.set_ylabel("Predicted age (months)", fontsize=9)
        ax.set_title(
            f"{sex_label}  (n={len(y)})\n"
            f"LOO R² = {r2:.3f}  |  MAE = {mae:.2f} mo",
            fontsize=10, fontweight="bold", color=color
        )
        ax.grid(color="#eeeeee", zorder=0)

    # R² bar chart
    X_all = sub[MULTIMODAL_FEATURES].values
    y_all = sub["Age_num"].values
    p_all = loo_predict_svr(X_all, y_all)
    r2_all, mae_all = loo_r2_mae(p_all, y_all)

    labels   = ["Combined\n(M+F)", "Male only", "Female only"]
    r2_vals  = [r2_all,
                sex_data.get("Male",   {}).get("r2",  np.nan),
                sex_data.get("Female", {}).get("r2",  np.nan)]
    mae_vals = [mae_all,
                sex_data.get("Male",   {}).get("mae", np.nan),
                sex_data.get("Female", {}).get("mae", np.nan)]
    colors   = [GREY, C_MALE, C_FEMALE]

    bars = ax_r2.bar(labels, r2_vals, color=colors, alpha=0.85,
                     edgecolor="white", zorder=3, width=0.5)
    for bar, r2v, maev in zip(bars, r2_vals, mae_vals):
        if not np.isnan(r2v):
            ax_r2.text(bar.get_x() + bar.get_width() / 2,
                       bar.get_height() + 0.01,
                       f"R²={r2v:.3f}\nMAE={maev:.2f}mo",
                       ha="center", va="bottom", fontsize=8, fontweight="bold")

    ax_r2.set_ylim(0, 1.05)
    ax_r2.set_ylabel("LOO-CV R²", fontsize=9)
    ax_r2.set_title("Model performance by sex", fontsize=10, fontweight="bold")
    ax_r2.axhline(0, color=GREY, lw=0.8)
    ax_r2.grid(axis="y", color="#dddddd", zorder=0)


# ── Figure 5: PCHIP parameter trajectories ───────────────────────────────────
def plot_pchip_trajectories(df: pd.DataFrame):
    """
    PCHIP splines of raw parameter values vs age, by sex.
    Panels: θ, κ, collagen straightness, distensibility
    (mirrors Figure 6 panels A–D in manuscript v8)
    """
    panels = [
        ("ORIENTATION-θ",           "Fiber orientation  θ  (°)",         "θ  (degrees)"),
        ("ORIENTATION-κ",           "Fiber concentration  κ",             "κ"),
        ("COLLAGEN STRAIGHTNESS",   "Collagen straightness",              "Straightness"),
        ("Distensibility (mmHg-1)", "PA distensibility  (mmHg⁻¹)",       "Distensibility (mmHg⁻¹)"),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(20, 5), constrained_layout=True)

    for ax, (col, title, ylabel) in zip(axes, panels):
        sub = df[["Age_num", "Sex", col]].dropna()
        age_nodes = np.array([6, 12, 18, 24], dtype=float)

        for sex_label, color in [("Male", C_MALE), ("Female", C_FEMALE)]:
            mask = sub["Sex"] == sex_label
            true_ages  = sub.loc[mask, "Age_num"].values
            param_vals = sub.loc[mask, col].values

            means, sds = [], []
            for age in age_nodes:
                bin_mask = true_ages == age
                means.append(param_vals[bin_mask].mean() if bin_mask.sum() > 0 else np.nan)
                sds.append(param_vals[bin_mask].std()    if bin_mask.sum() > 1 else 0.0)

            means = np.array(means, dtype=float)
            sds   = np.array(sds,   dtype=float)
            valid = ~np.isnan(means)
            if valid.sum() < 2:
                continue

            spl    = PchipInterpolator(age_nodes[valid], means[valid])
            spl_sd = PchipInterpolator(age_nodes[valid], sds[valid])
            x_plot = np.linspace(age_nodes[valid].min(), age_nodes[valid].max(), 300)
            y_c    = spl(x_plot)
            y_sd   = np.clip(spl_sd(x_plot), 0, None)

            ax.plot(x_plot, y_c, color=color, lw=2.5, label=sex_label, zorder=3)
            ax.fill_between(x_plot, y_c - y_sd, y_c + y_sd,
                            color=color, alpha=0.15, zorder=2)
            ax.scatter(age_nodes[valid], means[valid],
                       color=color, s=50, zorder=4, edgecolors="white", linewidths=0.5)

        ax.set_xlabel("Age (months)", fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_xticks([6, 12, 18, 24])
        ax.legend(fontsize=9)
        ax.grid(color="#eeeeee", zorder=0)

    fig.suptitle(
        "PCHIP parameter trajectories  ·  Male vs Female  ·  6–24mo\n"
        "Mean ± 1 SD per age group  ·  PCHIP spline interpolation",
        fontsize=12, fontweight="bold"
    )
    return fig


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "Aging_PA_merged_updated.xlsx"
    print(f"Loading data: {path}")
    df = load_data(path)

    n_photon = df[["Age_num"] + PHOTON_FEATURES].dropna().shape[0]
    n_multi  = df[["Age_num"] + MULTIMODAL_FEATURES].dropna().shape[0]
    print(f"  2-photon-only subset:  n={n_photon}")
    print(f"  Multimodal subset:     n={n_multi}")
    print(f"  Sex: {df.groupby('Sex').size().to_dict()}")
    print(f"  Timepoints: {sorted(df['Age_num'].dropna().unique().astype(int))}")

    # ── Figure 1: stepwise comparison ─────────────────────────────────────────
    print("\nFigure 1: stepwise model comparison …")
    fig1, ax1 = plt.subplots(figsize=(11, 6), constrained_layout=True)
    plot_stepwise_comparison(df, ax1)
    fig1.suptitle("Vascular Aging — Multimodal SVR  (6–24mo)",
                  fontsize=13, fontweight="bold")
    fig1.savefig("fig1_stepwise_comparison.png", dpi=150, bbox_inches="tight")
    print("  → fig1_stepwise_comparison.png")

    # ── Figure 2: predicted vs actual ─────────────────────────────────────────
    print("Figure 2: predicted vs actual …")
    fig2, ax2 = plt.subplots(figsize=(7, 6), constrained_layout=True)
    preds, y_true, sex = plot_predicted_vs_actual(df, ax2)
    fig2.savefig("fig2_predicted_vs_actual.png", dpi=150, bbox_inches="tight")
    print("  → fig2_predicted_vs_actual.png")

    # ── Figure 3: residuals ────────────────────────────────────────────────────
    print("Figure 3: residuals …")
    fig3, ax3 = plt.subplots(figsize=(9, 5), constrained_layout=True)
    plot_residuals(preds, y_true, sex, ax3)
    fig3.savefig("fig3_residuals.png", dpi=150, bbox_inches="tight")
    print("  → fig3_residuals.png")

    # ── Figure 4: sex-stratified ───────────────────────────────────────────────
    print("Figure 4: sex-stratified …")
    fig4 = plt.figure(figsize=(14, 10), constrained_layout=True)
    gs   = fig4.add_gridspec(2, 2)
    ax_m  = fig4.add_subplot(gs[0, 0])
    ax_f  = fig4.add_subplot(gs[0, 1])
    ax_r2 = fig4.add_subplot(gs[1, 0])
    plot_sex_stratified(df, ax_m, ax_f, ax_r2)
    fig4.suptitle("Sex-Stratified Multimodal SVR  (6–24mo)",
                  fontsize=13, fontweight="bold")
    fig4.savefig("fig4_sex_stratified.png", dpi=150, bbox_inches="tight")
    print("  → fig4_sex_stratified.png")

    # ── Figure 5: PCHIP trajectories ──────────────────────────────────────────
    print("Figure 5: PCHIP parameter trajectories …")
    fig5 = plot_pchip_trajectories(df)
    fig5.savefig("fig5_pchip_trajectories.png", dpi=150, bbox_inches="tight")
    print("  → fig5_pchip_trajectories.png")

    print("\nDone.")
    plt.show()


if __name__ == "__main__":
    main()
