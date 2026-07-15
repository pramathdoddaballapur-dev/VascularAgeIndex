"""
save_svr_model.py
-----------------
Trains the final multimodal SVR biological age clock on ALL available data
and saves the model + preprocessing pipeline to disk.

Use this ONCE after your LOO evaluation is done.
The saved artifacts are then used by predict_age.py for inference on new samples.

Preprocessing mirrors the per-fold pipeline in svr_multimodal_analysis.py:
    1. Winsorise at 5/95th percentile (fit on full training set)
    2. Yeo-Johnson power transform
    3. StandardScaler

Outputs:
    svr_final_model.pkl     — trained SVR (RBF, C=200, ε=1.5, γ=0.05)
    svr_preprocessor.pkl    — dict with winsor bounds, PowerTransformer, StandardScaler
    svr_feature_names.pkl   — ordered feature list (for column alignment check)
"""

import sys
import numpy as np
import joblib
from scipy.special import i0 as bessel_i0
import pandas as pd
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler, PowerTransformer

# ── Constants (must match svr_multimodal_analysis.py exactly) ─────────────────
AGE_MAP = {"3mo": 6, "6mo": 6, "12mo": 12, "18mo": 18, "24mo": 24}
REF_CIRC  = np.deg2rad(90)
REF_AXIAL = np.deg2rad(0)

BEST_C     = 200
BEST_EPS   = 1.5
BEST_GAMMA = 0.05

PHOTON_FEATURES     = ["vm_circ", "vm_axial", "COLLAGEN STRAIGHTNESS"]
MECH_FEATURES       = ["Distensibility (mmHg-1)",
                        "PA Circumferential Stiffness (MPa)",
                        "PWV"]
MULTIMODAL_FEATURES = PHOTON_FEATURES + MECH_FEATURES


def vm_pdf(x_rad, theta_deg, kappa):
    theta_rad = np.deg2rad(theta_deg)
    return np.exp(kappa * np.cos(x_rad - theta_rad)) / (2 * np.pi * bessel_i0(kappa))


def load_data(path: str) -> pd.DataFrame:
    df = pd.read_excel(path)
    df = df[df["Sample Name"].notna() & (df["Sample Name"] != "Average")]
    df = df[df["Age"].notna()]
    df["Age_num"] = df["Age"].map(AGE_MAP)
    df = df[df["Age_num"].notna()].copy()

    theta = df["ORIENTATION-θ"]
    kappa = df["ORIENTATION-κ"]
    valid = theta.notna() & kappa.notna()
    df.loc[valid, "vm_circ"]  = vm_pdf(REF_CIRC,  theta[valid], kappa[valid])
    df.loc[valid, "vm_axial"] = vm_pdf(REF_AXIAL, theta[valid], kappa[valid])
    return df


def train_and_save(data_path: str):
    print(f"Loading data from: {data_path}")
    df = load_data(data_path)

    sub = df[["Age_num"] + MULTIMODAL_FEATURES].dropna()
    X = sub[MULTIMODAL_FEATURES].values
    y = sub["Age_num"].values
    print(f"Training on n={len(y)} samples, {X.shape[1]} features")

    # ── Preprocessing (fit on all data — this is the final model, not LOO) ────
    lo = np.percentile(X, 5,  axis=0)
    hi = np.percentile(X, 95, axis=0)
    X_clipped = np.clip(X, lo, hi)

    pt = PowerTransformer(method='yeo-johnson')
    X_transformed = pt.fit_transform(X_clipped)

    sc = StandardScaler()
    X_scaled = sc.fit_transform(X_transformed)

    # ── Train final SVR ────────────────────────────────────────────────────────
    model = SVR(kernel="rbf", C=BEST_C, epsilon=BEST_EPS, gamma=BEST_GAMMA)
    model.fit(X_scaled, y)

    # Quick sanity check: in-sample R² (will be higher than LOO — expected)
    preds_insample = model.predict(X_scaled)
    ss_res = np.sum((y - preds_insample) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2_insample = 1 - ss_res / ss_tot
    print(f"In-sample R² (not LOO, expected to be high): {r2_insample:.3f}")
    print("  (Use LOO R²=0.834 from svr_multimodal_analysis.py as the reported metric)")

    # ── Save ──────────────────────────────────────────────────────────────────
    preprocessor = {
        "winsor_lo": lo,
        "winsor_hi": hi,
        "power_transformer": pt,
        "scaler": sc,
    }

    joblib.dump(model,              "svr_final_model.pkl")
    joblib.dump(preprocessor,       "svr_preprocessor.pkl")
    joblib.dump(MULTIMODAL_FEATURES, "svr_feature_names.pkl")

    print("\nSaved:")
    print("  svr_final_model.pkl")
    print("  svr_preprocessor.pkl")
    print("  svr_feature_names.pkl")
    print("\nDone. Use predict_age.py to run inference on new samples.")


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "Aging_PA_merged_updated.xlsx"
    train_and_save(path)
