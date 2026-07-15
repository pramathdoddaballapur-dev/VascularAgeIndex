"""
save_svr_2photon_model.py
-------------------------
Trains the 2-photon-only SVR biological age clock on all available data
and saves the model + preprocessing pipeline.

Performance (from manuscript, LOO-CV):
    n=81  |  LOO R2=0.596  |  MAE=3.43 mo

Features: vm_circ, vm_axial, COLLAGEN STRAIGHTNESS
(vm_circ and vm_axial computed from raw ORIENTATION-theta and ORIENTATION-kappa)

Outputs:
    svr_2photon_model.pkl          -- trained SVR (RBF, C=200, eps=1.5, gamma=0.05)
    svr_2photon_preprocessor.pkl   -- winsor bounds, PowerTransformer, StandardScaler
    svr_2photon_feature_names.pkl  -- ordered feature list
"""

import sys
import numpy as np
import joblib
import pandas as pd
from scipy.special import i0 as bessel_i0
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler, PowerTransformer

# Constants (must match svr_multimodal_analysis.py exactly)
AGE_MAP    = {"3mo": 6, "6mo": 6, "12mo": 12, "18mo": 18, "24mo": 24}
REF_CIRC   = np.deg2rad(90)
REF_AXIAL  = np.deg2rad(0)

BEST_C     = 200
BEST_EPS   = 1.5
BEST_GAMMA = 0.05

PHOTON_FEATURES = ["vm_circ", "vm_axial", "COLLAGEN STRAIGHTNESS"]


def vm_pdf(x_rad, theta_deg, kappa):
    theta_rad = np.deg2rad(theta_deg)
    return np.exp(kappa * np.cos(x_rad - theta_rad)) / (2 * np.pi * bessel_i0(kappa))


def load_data(path):
    df = pd.read_excel(path)
    df = df[df["Sample Name"].notna() & (df["Sample Name"] != "Average")]
    df = df[df["Age"].notna()]
    df["Age_num"] = df["Age"].map(AGE_MAP)
    df = df[df["Age_num"].notna()].copy()

    theta = df["ORIENTATION-\u03b8"]
    kappa = df["ORIENTATION-\u03ba"]
    valid = theta.notna() & kappa.notna()
    df.loc[valid, "vm_circ"]  = vm_pdf(REF_CIRC,  theta[valid], kappa[valid])
    df.loc[valid, "vm_axial"] = vm_pdf(REF_AXIAL, theta[valid], kappa[valid])
    return df


def train_and_save(data_path):
    print(f"Loading data from: {data_path}")
    df = load_data(data_path)

    sub = df[["Age_num"] + PHOTON_FEATURES].dropna()
    X = sub[PHOTON_FEATURES].values
    y = sub["Age_num"].values
    print(f"Training on n={len(y)} samples, {X.shape[1]} features")
    print(f"Features: {PHOTON_FEATURES}")

    if len(y) != 81:
        print(f"  Warning: expected n=81 for 2-photon model, got n={len(y)}.")

    # Preprocessing
    lo = np.percentile(X, 5,  axis=0)
    hi = np.percentile(X, 95, axis=0)
    X_clipped     = np.clip(X, lo, hi)
    pt = PowerTransformer(method='yeo-johnson')
    X_transformed = pt.fit_transform(X_clipped)
    sc = StandardScaler()
    X_scaled      = sc.fit_transform(X_transformed)

    # Train
    model = SVR(kernel="rbf", C=BEST_C, epsilon=BEST_EPS, gamma=BEST_GAMMA)
    model.fit(X_scaled, y)

    # Sanity check
    preds = model.predict(X_scaled)
    r2 = 1 - np.sum((y - preds)**2) / np.sum((y - y.mean())**2)
    print(f"In-sample R2 (not LOO, expected to be high): {r2:.3f}")
    print("  (Manuscript LOO R2=0.596 from svr_multimodal_analysis.py is the reported metric)")

    # Save
    preprocessor = {
        "winsor_lo": lo,
        "winsor_hi": hi,
        "power_transformer": pt,
        "scaler": sc,
    }

    joblib.dump(model,           "svr_2photon_model.pkl")
    joblib.dump(preprocessor,    "svr_2photon_preprocessor.pkl")
    joblib.dump(PHOTON_FEATURES, "svr_2photon_feature_names.pkl")

    print("\nSaved:")
    print("  svr_2photon_model.pkl")
    print("  svr_2photon_preprocessor.pkl")
    print("  svr_2photon_feature_names.pkl")
    print("\nDone. Use predict_age.py --2photon to run inference on new samples.")


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "Aging_PA_merged_updated.xlsx"
    train_and_save(path)
