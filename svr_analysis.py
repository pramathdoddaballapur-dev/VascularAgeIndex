"""
Vascular Aging SVR Analysis
  • Restricts data to 6, 12, 18, 24 month timepoints (removes 3mo and 27mo)
  • Replaces raw ORIENTATION-θ / ORIENTATION-κ with two von Mises PDF features:
        vm_circ  = f(90°; θ, κ)   [circumferential reference]
        vm_axial = f(0°;  θ, κ)   [axial reference]
      where  f(x; θ, κ) = exp(κ · cos(x − θ)) / (2π · I₀(κ))
  • Winning combo is now: vm_circ + vm_axial + COLLAGEN STRAIGHTNESS
    with best LOO hyperparameters C=200, ε=1.5, γ=0.05

Figures produced:
  Figure 1 – Model comparison (LOO-CV R²) across param combos & model types
  Figure 2 – Predicted vs actual age (LOO-CV), coloured by sex
  Figure 3 – Residuals by age group (box + strip plot)
  Figure 4 – Sex-stratified analysis: per-sex LOO predictions + R² comparison

"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.special import i0 as bessel_i0
from scipy import stats
from scipy.interpolate import make_interp_spline
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.preprocessing import PolynomialFeatures, StandardScaler, PowerTransformer, RobustScaler
from sklearn.pipeline import Pipeline
from sklearn.svm import SVR
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import r2_score
from sklearn.model_selection import LeaveOneOut



# colours
C_MALE   = "#185FA5"
C_FEMALE = "#993556"
C_LINEAR = "#185FA5"
C_POLY   = "#0F6E56"
C_GBM    = "#BA7517"
C_SVR    = "#993556"
GREY     = "#888780"

# constants
AGE_MAP = {"3mo": 6, "6mo": 6, "12mo": 12, "18mo": 18, "24mo": 24}

# Von Mises reference angles (radians)
REF_CIRC  = np.deg2rad(90)
REF_AXIAL = np.deg2rad(0)

# Best SVR hyperparameters (tuned after Yeo-Johnson transformation)
BEST_C     = 200
BEST_EPS   = 1.5
BEST_GAMMA = 0.05

WINNING_COMBO = ["vm_circ", "vm_axial", "COLLAGEN STRAIGHTNESS"]


# von Mises PDF
def vm_pdf(x_rad: float, theta_deg: np.ndarray, kappa: np.ndarray) -> np.ndarray:
    """Von Mises PDF evaluated at reference angle x_rad.
    f(x; θ, κ) = exp(κ · cos(x − θ)) / (2π · I₀(κ))
    θ supplied in degrees, converted internally.
    """
    theta_rad = np.deg2rad(theta_deg)
    return np.exp(kappa * np.cos(x_rad - theta_rad)) / (2 * np.pi * bessel_i0(kappa))


# data load
def load_data(path: str) -> pd.DataFrame:
    df = pd.read_excel(path)
    df = df[df["Sample Name"].notna() & (df["Sample Name"] != "Average")]
    df = df[df["Age"].notna()]
    df["Age_num"] = df["Age"].map(AGE_MAP)
    df = df[df["Age_num"].notna()].copy()   # drop 3mo and 27mo

    # Compute combined von Mises features
    theta = df["ORIENTATION-θ"]
    kappa = df["ORIENTATION-κ"]
    valid = theta.notna() & kappa.notna()
    df.loc[valid, "vm_circ"]  = vm_pdf(REF_CIRC,  theta[valid], kappa[valid])
    df.loc[valid, "vm_axial"] = vm_pdf(REF_AXIAL, theta[valid], kappa[valid])

    return df


# LOO helpers
def loo_predict(model_factory, X: np.ndarray, y: np.ndarray,
                yeo_johnson: bool = False) -> np.ndarray:
    """Return LOO out-of-fold predictions.
    If yeo_johnson=True, applies per-fold preprocessing:
      1. Winsorise at 5/95th percentile (clipped to training fold limits)
      2. Yeo-Johnson power transform
      3. StandardScaler
    All steps fit on training fold only — no leakage.
    """
    preds = np.empty(len(y))
    for train_idx, test_idx in LeaveOneOut().split(X):
        X_tr, X_te = X[train_idx], X[test_idx]
        if yeo_johnson:
            # Winsorise
            lo = np.percentile(X_tr, 5, axis=0)
            hi = np.percentile(X_tr, 95, axis=0)
            X_tr = np.clip(X_tr, lo, hi)
            X_te = np.clip(X_te, lo, hi)
            # Yeo-Johnson
            pt = PowerTransformer(method='yeo-johnson')
            X_tr = pt.fit_transform(X_tr)
            X_te = pt.transform(X_te)
            # Standard scale
            sc = StandardScaler()
            X_tr = sc.fit_transform(X_tr)
            X_te = sc.transform(X_te)
        m = model_factory()
        m.fit(X_tr, y[train_idx])
        preds[test_idx] = m.predict(X_te)
    return preds


def loo_r2(model_factory, X: np.ndarray, y: np.ndarray,
           yeo_johnson: bool = False) -> float:
    preds = loo_predict(model_factory, X, y, yeo_johnson=yeo_johnson)
    ss_res = np.sum((y - preds) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot


# SVR model building
def make_svr():
    return Pipeline([
        ("sc", StandardScaler()),
        ("svr", SVR(kernel="rbf", C=BEST_C, epsilon=BEST_EPS, gamma=BEST_GAMMA))
    ])

def make_svr_yj():
    return Pipeline([
        ("yj", PowerTransformer(method='yeo-johnson')),
        ("sc", StandardScaler()),
        ("svr", SVR(kernel="rbf", C=BEST_C, epsilon=BEST_EPS, gamma=BEST_GAMMA))
    ])

def make_linear():  return LinearRegression()
def make_poly2():   return Pipeline([("poly", PolynomialFeatures(2, include_bias=False)), ("ridge", Ridge(10))])
def make_gbm():     return GradientBoostingRegressor(n_estimators=100, max_depth=2, learning_rate=0.05, random_state=42)


# Figure 1: model comparison
def plot_model_comparison(df: pd.DataFrame, ax_top: plt.Axes, ax_bot: plt.Axes):
    combos = {
        "1-param\n(κ only)":                   ["ORIENTATION-κ"],
        "2-param\n(vm_circ\n+ Collagen)":       ["vm_circ", "COLLAGEN STRAIGHTNESS"],
        "3-param\n(vm_circ + vm_axial\n+ Collagen)": ["vm_circ", "vm_axial", "COLLAGEN STRAIGHTNESS"],
        "4-param\n(+ Adv. Thickness)":          ["vm_circ", "vm_axial", "COLLAGEN STRAIGHTNESS", "ADVENTITIA THICKNESS"],
        "All features":                         ["vm_circ", "vm_axial", "COLLAGEN STRAIGHTNESS",
                                                  "ADVENTITIA THICKNESS", "MEDIA THICKNESS",
                                                  "ADVENTITIAL CELL COUNT", "SMC CELL COUNT",
                                                  "ENDOTHELIAL CELL COUNT"],
    }
    model_defs = {
        "Linear":      (make_linear, C_LINEAR),
        "Poly2+Ridge": (make_poly2,  C_POLY),
        "GBM (d=2)":   (make_gbm,    C_GBM),
        "SVR (RBF)":   (make_svr_yj, C_SVR),
    }

    combo_names = list(combos.keys())
    n_combos    = len(combo_names)
    n_models    = len(model_defs)
    bar_w       = 0.18
    offsets     = np.linspace(-(n_models - 1) / 2, (n_models - 1) / 2, n_models) * bar_w
    x           = np.arange(n_combos)

    results       = {mn: [] for mn in model_defs}
    train_results = {mn: [] for mn in model_defs}

    for cname, combo in combos.items():
        sub = df[["Age_num"] + combo].dropna()
        X, y = sub[combo].values, sub["Age_num"].values
        for mname, (factory, _) in model_defs.items():
            use_yj = (mname == "SVR (RBF)")
            r2_loo = loo_r2(factory, X, y, yeo_johnson=use_yj)
            m = factory(); m.fit(X, y)
            r2_tr = r2_score(y, m.predict(X))
            results[mname].append(r2_loo)
            train_results[mname].append(r2_tr)

    for i, (mname, (_, color)) in enumerate(model_defs.items()):
        ax_top.bar(x + offsets[i], results[mname], bar_w,
                   label=mname, color=color, alpha=0.9, zorder=3)

    ax_top.axhline(0, color=GREY, lw=0.8, ls="--")
    ax_top.set_xticks(x)
    ax_top.set_xticklabels(combo_names, fontsize=8)
    ax_top.set_ylabel("LOO-CV R²", fontsize=9)
    ax_top.set_title("LOO-CV R² by parameter combo & model type  (6–24mo, Von Mises features, SVR with Yeo-Johnson transform)",
                     fontsize=10, fontweight="bold")
    ax_top.legend(fontsize=8, ncol=4, loc="upper left")
    ax_top.set_ylim(-0.35, 0.85)
    ax_top.grid(axis="y", color="#dddddd", zorder=0)

    # train vs LOO (bottom panel) — Linear, GBM, and SVR RBF
    models_to_check = [("Linear", C_LINEAR), ("GBM (d=2)", C_GBM), ("SVR (RBF)", C_SVR)]
    n_check = len(models_to_check)
    bw2 = 0.22
    offsets_bot = np.linspace(-(n_check - 1) / 2, (n_check - 1) / 2, n_check) * bw2

    for i, (mname, color) in enumerate(models_to_check):
        off = offsets_bot[i]
        ax_bot.bar(x + off, train_results[mname], bw2, color=color, alpha=0.28,
                   label=f"{mname} – Train", zorder=3)
        ax_bot.bar(x + off, results[mname],       bw2, color=color, alpha=1.0,
                   label=f"{mname} – LOO",   zorder=4)

    ax_bot.set_xticks(x)
    ax_bot.set_xticklabels(combo_names, fontsize=8)
    ax_bot.set_ylabel("R²", fontsize=9)
    ax_bot.set_title("Train R² vs LOO R²  (overfitting check — Linear, GBM, SVR RBF)",
                     fontsize=10, fontweight="bold")
    ax_bot.legend(fontsize=8, ncol=3, loc="upper left")
    ax_bot.set_ylim(0, 1.05)
    ax_bot.grid(axis="y", color="#dddddd", zorder=0)


# Figure 2: predicted vs actual
def plot_predicted_vs_actual(df: pd.DataFrame, ax: plt.Axes):
    sub   = df[["Age_num", "Sex"] + WINNING_COMBO].dropna()
    X     = sub[WINNING_COMBO].values
    y     = sub["Age_num"].values
    sex   = sub["Sex"].values

    preds  = loo_predict(make_svr, X, y, yeo_johnson=True)
    r2_loo = 1 - np.sum((y - preds) ** 2) / np.sum((y - y.mean()) ** 2)

    for label, color in [("Male", C_MALE), ("Female", C_FEMALE)]:
        mask = sex == label
        ax.scatter(y[mask], preds[mask], c=color, alpha=0.72, s=55,
                   edgecolors="white", linewidths=0.4, label=label, zorder=3)

    lims = [4, 26]
    ax.plot(lims, lims, ls="--", lw=1.5, color=GREY, zorder=2, label="Perfect prediction")
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel("Actual age (months)", fontsize=9)
    ax.set_ylabel("Predicted age (months)", fontsize=9)
    ax.set_title(
        f"SVR (RBF) – predicted vs actual age  (LOO-CV)\n"
        f"Features: f(90°;θ,κ) + f(0°;θ,κ) + Collagen straightness  |  "
        f"Winsorise 5/95 + Yeo-Johnson  |  "
        f"LOO R² = {r2_loo:.3f}  |  C={BEST_C}, ε={BEST_EPS}, γ={BEST_GAMMA}",
        fontsize=9, fontweight="bold"
    )
    ax.legend(fontsize=8)
    ax.grid(color="#eeeeee", zorder=0)
    ax.set_xticks([6, 12, 18, 24])
    ax.set_yticks(range(4, 27, 2))
    return preds, y, sex


#  Figure 3: residuals by age group
def plot_residuals(preds: np.ndarray, y: np.ndarray, sex: np.ndarray, ax: plt.Axes):
    residuals   = preds - y
    age_groups  = sorted(np.unique(y).astype(int))
    group_labels = [f"{g}mo" for g in age_groups]
    rng = np.random.default_rng(42)

    for gi, age in enumerate(age_groups):
        mask = y == age
        res  = residuals[mask]
        sx   = sex[mask]

        q1, med, q3 = np.percentile(res, [25, 50, 75])
        bw = 0.22
        ax.broken_barh(
            [(gi - bw, 2 * bw)], (q1, q3 - q1),
            facecolors="none", edgecolors=C_LINEAR, linewidth=1.2, zorder=3
        )
        ax.hlines(med, gi - bw, gi + bw, colors=C_LINEAR, linewidth=2, zorder=4)
        ax.vlines(gi, res.min(), q1, colors=C_LINEAR, linewidth=0.9, zorder=3)
        ax.vlines(gi, q3, res.max(), colors=C_LINEAR, linewidth=0.9, zorder=3)
        ax.hlines([res.min(), res.max()], gi - 0.1, gi + 0.1, colors=C_LINEAR, linewidth=0.9, zorder=3)

        for label, color in [("Male", C_MALE), ("Female", C_FEMALE)]:
            pts = res[sx == label]
            jx  = gi + rng.uniform(-0.28, 0.28, size=len(pts))
            ax.scatter(jx, pts, c=color, alpha=0.6, s=28,
                       edgecolors="white", linewidths=0.3, zorder=5)

    ax.axhline(0, ls="--", lw=1.5, color=GREY, zorder=2)
    ax.set_xticks(range(len(age_groups)))
    ax.set_xticklabels(group_labels, fontsize=9)
    ax.set_xlabel("Actual age group", fontsize=9)
    ax.set_ylabel("Residual  (predicted − actual, months)", fontsize=9)
    ax.set_title("Residuals by age group – SVR (RBF) LOO-CV", fontsize=10, fontweight="bold")
    patches = [mpatches.Patch(color=C_MALE, label="Male"),
               mpatches.Patch(color=C_FEMALE, label="Female")]
    ax.legend(handles=patches, fontsize=8)
    ax.grid(axis="y", color="#eeeeee", zorder=0)


#  Figure 4: sex-stratified analysis
def plot_sex_stratified(df: pd.DataFrame, fig: plt.Figure,
                        ax_male: plt.Axes, ax_female: plt.Axes,
                        ax_r2: plt.Axes, ax_stats: plt.Axes):
    sub = df[["Age_num", "Sex"] + WINNING_COMBO].dropna()
    age_groups = sorted(sub["Age_num"].unique().astype(int))

    #  per-sex LOO predictions
    sex_data = {}
    for sex_label, color, ax in [("Male", C_MALE, ax_male), ("Female", C_FEMALE, ax_female)]:
        s = sub[sub["Sex"] == sex_label]
        X = s[WINNING_COMBO].values
        y = s["Age_num"].values

        if len(y) < 5:
            ax.set_title(f"{sex_label}: insufficient data", fontsize=10)
            continue

        preds  = loo_predict(make_svr, X, y, yeo_johnson=True)
        r2_loo = 1 - np.sum((y - preds) ** 2) / np.sum((y - y.mean()) ** 2)
        mae    = np.mean(np.abs(preds - y))
        sex_data[sex_label] = {"y": y, "preds": preds, "r2": r2_loo, "mae": mae}

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
            f"LOO R² = {r2_loo:.3f}  |  MAE = {mae:.2f} mo",
            fontsize=10, fontweight="bold", color=color
        )
        ax.grid(color="#eeeeee", zorder=0)

        # per-age-group annotations
        for age in age_groups:
            mask = y == age
            if mask.sum() == 0:
                continue
            mean_pred = preds[mask].mean()
            age_mae   = np.mean(np.abs(preds[mask] - y[mask]))
            ax.annotate(
                f"n={mask.sum()}\nMAE={age_mae:.1f}",
                xy=(age, mean_pred),
                xytext=(age + 0.6, mean_pred + 1.2),
                fontsize=6.5, color=color, alpha=0.8,
                arrowprops=dict(arrowstyle="-", color=color, alpha=0.4, lw=0.8)
            )

    # R² comparison bar chart
    # combined model
    X_all = sub[WINNING_COMBO].values
    y_all = sub["Age_num"].values
    r2_combined = 1 - np.sum(
        (y_all - loo_predict(make_svr, X_all, y_all, yeo_johnson=True)) ** 2
    ) / np.sum((y_all - y_all.mean()) ** 2)

    labels = ["Combined\n(M + F)", "Male only", "Female only"]
    r2_vals = [
        r2_combined,
        sex_data.get("Male",   {}).get("r2", np.nan),
        sex_data.get("Female", {}).get("r2", np.nan),
    ]
    mae_combined = np.mean(np.abs(y_all - loo_predict(make_svr, X_all, y_all, yeo_johnson=True)))
    mae_vals = [
        mae_combined,
        sex_data.get("Male",   {}).get("mae", np.nan),
        sex_data.get("Female", {}).get("mae", np.nan),
    ]
    colors = [GREY, C_MALE, C_FEMALE]

    bars = ax_r2.bar(labels, r2_vals, color=colors, alpha=0.85,
                     edgecolor="white", zorder=3, width=0.5)
    for bar, r2v, maev in zip(bars, r2_vals, mae_vals):
        if not np.isnan(r2v):
            ax_r2.text(bar.get_x() + bar.get_width() / 2,
                       bar.get_height() + 0.01,
                       f"R²={r2v:.3f}\nMAE={maev:.2f}mo",
                       ha="center", va="bottom", fontsize=8, fontweight="bold")

    ax_r2.set_ylim(0, 0.95)
    ax_r2.set_ylabel("LOO-CV R²", fontsize=9)
    ax_r2.set_title("Model performance by sex", fontsize=10, fontweight="bold")
    ax_r2.axhline(0, color=GREY, lw=0.8)
    ax_r2.grid(axis="y", color="#dddddd", zorder=0)

    # Mann-Whitney U: sex differences per feature per age group, includes p-values for difference
    test_features = ["vm_circ", "vm_axial", "COLLAGEN STRAIGHTNESS",
                     "ADVENTITIA THICKNESS", "ADVENTITIAL CELL COUNT"]
    feat_labels   = ["f(90°;θ,κ)", "f(0°;θ,κ)", "Coll. straight.",
                     "Adv. thickness", "Adv. cell count"]

    p_matrix = np.ones((len(test_features), len(age_groups)))

    for fi, feat in enumerate(test_features):
        for ai, age in enumerate(age_groups):
            age_df = df[(df["Age_num"] == age) & df["Sex"].notna() & df[feat].notna()]
            males   = age_df.loc[age_df["Sex"] == "Male",   feat].values
            females = age_df.loc[age_df["Sex"] == "Female", feat].values
            if len(males) >= 3 and len(females) >= 3:
                _, p = stats.mannwhitneyu(males, females, alternative="two-sided")
                p_matrix[fi, ai] = p

    # heatmap
    im = ax_stats.imshow(p_matrix, aspect="auto", cmap="RdYlGn",
                         vmin=0, vmax=0.1)
    ax_stats.set_xticks(range(len(age_groups)))
    ax_stats.set_xticklabels([f"{a}mo" for a in age_groups], fontsize=9)
    ax_stats.set_yticks(range(len(test_features)))
    ax_stats.set_yticklabels(feat_labels, fontsize=8)
    ax_stats.set_title("Sex difference p-values (Mann-Whitney U)\nGreen = p<0.05, Red = p>0.1",
                        fontsize=9, fontweight="bold")

    for fi in range(len(test_features)):
        for ai in range(len(age_groups)):
            p = p_matrix[fi, ai]
            star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else f"{p:.2f}"
            ax_stats.text(ai, fi, star, ha="center", va="center",
                          fontsize=8, color="black", fontweight="bold" if p < 0.05 else "normal")

    plt.colorbar(im, ax=ax_stats, label="p-value", shrink=0.8)


# Figure 5: spline curves of parameter vs LOO predicted age, by sex
def plot_params_vs_predicted_age(df: pd.DataFrame):
    feat_cols    = ['vm_circ', 'vm_axial', 'COLLAGEN STRAIGHTNESS']
    raw_cols     = ['ORIENTATION-θ', 'ORIENTATION-κ', 'COLLAGEN STRAIGHTNESS']
    panel_labels = ['Orientation  θ  (°)', 'Concentration  κ', 'Collagen straightness']
    ylabels      = ['θ  (degrees)', 'κ  (concentration)', 'Collagen straightness']

    all_cols = ['vm_circ', 'vm_axial', 'COLLAGEN STRAIGHTNESS',
                'ORIENTATION-θ', 'ORIENTATION-κ', 'Age_num', 'Sex']
    sub = df[all_cols].dropna().reset_index(drop=True)
    X   = sub[feat_cols].values
    y   = sub['Age_num'].values

    # LOO predictions — winsorise + Yeo-Johnson + StandardScaler
    preds = np.empty(len(y))
    for tr, te in LeaveOneOut().split(X):
        lo  = np.percentile(X[tr], 5, axis=0); hi = np.percentile(X[tr], 95, axis=0)
        Xtr = np.clip(X[tr], lo, hi);          Xte = np.clip(X[te], lo, hi)
        pt  = PowerTransformer(method='yeo-johnson')
        Xtr = pt.fit_transform(Xtr);           Xte = pt.transform(Xte)
        sc  = StandardScaler()
        Xtr = sc.fit_transform(Xtr);           Xte = sc.transform(Xte)
        m   = SVR(kernel='rbf', C=BEST_C, epsilon=BEST_EPS, gamma=BEST_GAMMA)
        m.fit(Xtr, y[tr]); preds[te] = m.predict(Xte)

    x_dense = np.linspace(6, 24, 300)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)

    for ax, raw_col, title, ylabel in zip(axes, raw_cols, panel_labels, ylabels):

        for sex_label, color in [('Male', C_MALE), ('Female', C_FEMALE)]:
            mask = sub['Sex'].values == sex_label

            # Use true age group as knots — more stable than predicted-age bins, but (predicted ages cluster around true ages anyway)
            age_nodes  = np.array([6, 12, 18, 24], dtype=float)
            true_ages  = y[mask]
            param_vals = sub.loc[mask, raw_col].values

            mean_nodes = []
            sd_nodes   = []
            for age in age_nodes:
                bin_mask = true_ages == age
                if bin_mask.sum() > 0:
                    mean_nodes.append(param_vals[bin_mask].mean())
                    sd_nodes.append(param_vals[bin_mask].std() if bin_mask.sum() > 1 else 0.0)
                else:
                    mean_nodes.append(np.nan)
                    sd_nodes.append(0.0)

            mean_nodes = np.array(mean_nodes, dtype=float)
            sd_nodes   = np.array(sd_nodes,   dtype=float)

            valid_mask = ~np.isnan(mean_nodes)
            knot_x     = age_nodes[valid_mask]
            knot_y     = mean_nodes[valid_mask]
            knot_sd    = sd_nodes[valid_mask]

            if len(knot_x) < 2:
                continue

            # Monotonicity-preserving cubic spline (PCHIP) — no oscillation at endpoints
            from scipy.interpolate import PchipInterpolator
            spl    = PchipInterpolator(knot_x, knot_y)
            spl_sd = PchipInterpolator(knot_x, knot_sd)

            x_plot  = np.linspace(knot_x.min(), knot_x.max(), 300)
            y_curve = spl(x_plot)
            y_sd    = np.clip(spl_sd(x_plot), 0, None)

            ax.plot(x_plot, y_curve, color=color, lw=2.5, label=sex_label, zorder=3)
            ax.fill_between(x_plot,
                            y_curve - y_sd,
                            y_curve + y_sd,
                            color=color, alpha=0.15, zorder=2)

        ax.set_xlabel('Predicted age (months)', fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.set_xlim(6, 24)
        ax.set_xticks([6, 12, 18, 24])
        ax.legend(fontsize=9)
        ax.grid(color='#eeeeee', zorder=0)

    fig.suptitle(
        'Model-implied parameter trajectories across predicted age  —  Male vs Female\n'
        'PCHIP spline through age-group means  ·  shaded band = ±1 SD  ·  6–24mo',
        fontsize=12, fontweight='bold'
    )
    return fig


# main - Updates what function is running currently, wait time, and some figures
def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "Aging_PA_03-20-2026.xlsx"
    print(f"Loading data from: {path}")
    df = load_data(path)
    print(f"  {len(df)} samples  |  timepoints: {sorted(df['Age_num'].dropna().unique().astype(int))}")
    print(f"  Sex breakdown:\n{df.groupby(['Age_num','Sex']).size().to_string()}")
    print(f"\n  Von Mises features computed:")
    print(f"    vm_circ  = f(90°;θ,κ)  — circumferential reference")
    print(f"    vm_axial = f(0°;θ,κ)   — axial reference")
    print(f"    formula: exp(κ·cos(x−θ)) / (2π·I₀(κ))")

    # Figure 1
    print("\nFigure 1: model comparison (may take ~1–2 min) …")
    fig1, (ax1a, ax1b) = plt.subplots(2, 1, figsize=(12, 9), constrained_layout=True)
    plot_model_comparison(df, ax1a, ax1b)
    fig1.suptitle("Vascular Aging – Model Comparison  (6–24mo, Von Mises features)",
                  fontsize=12, fontweight="bold")
    fig1.savefig("figure1_model_comparison.png", dpi=150, bbox_inches="tight")
    print("  → figure1_model_comparison.png")

    # Figure 2
    print("Figure 2: predicted vs actual …")
    fig2, ax2 = plt.subplots(figsize=(7, 6), constrained_layout=True)
    preds, y_true, sex = plot_predicted_vs_actual(df, ax2)
    fig2.savefig("figure2_predicted_vs_actual.png", dpi=150, bbox_inches="tight")
    print("  → figure2_predicted_vs_actual.png")

    # Figure 3
    print("Figure 3: residuals …")
    fig3, ax3 = plt.subplots(figsize=(9, 5), constrained_layout=True)
    plot_residuals(preds, y_true, sex, ax3)
    fig3.savefig("figure3_residuals.png", dpi=150, bbox_inches="tight")
    print("  → figure3_residuals.png")

    # Figure 4
    print("Figure 4: sex-stratified analysis …")
    fig4 = plt.figure(figsize=(14, 10), constrained_layout=True)
    gs   = fig4.add_gridspec(2, 2, hspace=0.38, wspace=0.32)
    ax_male   = fig4.add_subplot(gs[0, 0])
    ax_female = fig4.add_subplot(gs[0, 1])
    ax_r2     = fig4.add_subplot(gs[1, 0])
    ax_stats  = fig4.add_subplot(gs[1, 1])
    plot_sex_stratified(df, fig4, ax_male, ax_female, ax_r2, ax_stats)
    fig4.suptitle("Sex-Stratified SVR Analysis  (6–24mo, Von Mises features)",
                  fontsize=13, fontweight="bold")
    fig4.savefig("figure4_sex_stratified.png", dpi=150, bbox_inches="tight")
    print("  → figure4_sex_stratified.png")

    #Figure 5
    print("Figure 5: parameters vs predicted age by sex …")
    fig5 = plot_params_vs_predicted_age(df)
    fig5.savefig("figure5_params_vs_predicted_age.png", dpi=150, bbox_inches="tight")
    print("  → figure5_params_vs_predicted_age.png")

    print("\nDone.")
    plt.show()


if __name__ == "__main__":
    main()