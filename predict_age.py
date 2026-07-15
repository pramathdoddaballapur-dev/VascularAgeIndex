"""
predict_age.py
--------------
Loads a saved SVR biological age clock and predicts age on new test samples.

Three models supported:
    Primary multimodal (default): 2-photon + PA mechanics (LOO R2=0.834)
    Lung model (--lung):          2-photon + Rn           (LOO R2=0.854)
    2-photon only (--2photon):    vm_circ + vm_axial + straightness (LOO R2=0.596)

Usage:
    python predict_age.py new_samples.xlsx
    python predict_age.py new_samples.xlsx --lung
    python predict_age.py new_samples.xlsx --2photon
    python predict_age.py new_samples.xlsx --2photon --sheet Sheet2
    python predict_age.py new_samples.xlsx --2photon --output predictions.csv

Required columns - primary multimodal:
    ORIENTATION-theta, ORIENTATION-kappa, COLLAGEN STRAIGHTNESS,
    Distensibility (mmHg-1), PA Circumferential Stiffness (MPa), PWV

Required columns - lung model:
    ORIENTATION-theta, ORIENTATION-kappa, COLLAGEN STRAIGHTNESS,
    Rn  (alias for: lung Resistance (Rn) (cmH2O.s/mL))

Required columns - 2-photon only:
    ORIENTATION-theta, ORIENTATION-kappa, COLLAGEN STRAIGHTNESS

Optional columns (any model):
    Sample Name, Age, Sex, Group

Prerequisites:
    Primary model:  run save_svr_model.py first
    Lung model:     run save_svr_lung_model.py first
    2-photon model: run save_svr_2photon_model.py first
"""

import argparse
import numpy as np
import pandas as pd
import joblib
from scipy.special import i0 as bessel_i0

AGE_MAP   = {"3mo": 6, "6mo": 6, "12mo": 12, "18mo": 18, "24mo": 24}
REF_CIRC  = np.deg2rad(90)
REF_AXIAL = np.deg2rad(0)

COLUMN_ALIASES = {
    "Rn": "lung Resistance (Rn) (cmH2O.s/mL)",
    "Distensibility (MPa-1)": "Distensibility (mmHg-1)"
}

# Extract numeric age from group label e.g. "Hypoxic 8" -> 8.0, "Hypoxic 5.5" -> 5.5
def parse_group_age(group_label):
    import re
    m = re.search(r"[\d]+\.?[\d]*", str(group_label))
    return float(m.group()) if m else np.nan


def vm_pdf(x_rad, theta_deg, kappa):
    theta_rad = np.deg2rad(theta_deg)
    return np.exp(kappa * np.cos(x_rad - theta_rad)) / (2 * np.pi * bessel_i0(kappa))


def load_and_prep(path, sheet=0):
    df = pd.read_excel(path, sheet_name=sheet)
    df.rename(columns=COLUMN_ALIASES, inplace=True)

    theta_col = next((c for c in df.columns if "theta" in c.lower() or "\u03b8" in c), None)
    kappa_col = next((c for c in df.columns if "kappa" in c.lower() or "\u03ba" in c), None)

    if theta_col and kappa_col:
        theta = df[theta_col]
        kappa = df[kappa_col]
        valid = theta.notna() & kappa.notna()
        df.loc[valid, "vm_circ"]  = vm_pdf(REF_CIRC,  theta[valid], kappa[valid])
        df.loc[valid, "vm_axial"] = vm_pdf(REF_AXIAL, theta[valid], kappa[valid])
    elif "vm_circ" not in df.columns:
        raise ValueError(
            "Input must have ORIENTATION-theta and ORIENTATION-kappa columns "
            "or pre-computed vm_circ and vm_axial columns."
        )

    if "Age" in df.columns:
        df["Age_num"] = df["Age"].map(AGE_MAP)

    # Parse actual age from Group label if present and Age_num not available
    if "Group" in df.columns and "Age_num" not in df.columns:
        df["Age_num"] = df["Group"].apply(parse_group_age)

    return df


def preprocess(X, preprocessor):
    lo = preprocessor["winsor_lo"]
    hi = preprocessor["winsor_hi"]
    pt = preprocessor["power_transformer"]
    sc = preprocessor["scaler"]
    return sc.transform(pt.transform(np.clip(X, lo, hi)))


def predict(data_path, model_type="multimodal", sheet=0, output_path=None):

    if model_type == "lung":
        model_file  = "svr_lung_model.pkl"
        prep_file   = "svr_lung_preprocessor.pkl"
        feat_file   = "svr_lung_feature_names.pkl"
        model_label = "Lung model (2-photon + Rn)  |  LOO R2=0.854"
    elif model_type == "2photon":
        model_file  = "svr_2photon_model.pkl"
        prep_file   = "svr_2photon_preprocessor.pkl"
        feat_file   = "svr_2photon_feature_names.pkl"
        model_label = "2-photon only (vm_circ + vm_axial + straightness)  |  LOO R2=0.596"
    else:
        model_file  = "svr_final_model.pkl"
        prep_file   = "svr_preprocessor.pkl"
        feat_file   = "svr_feature_names.pkl"
        model_label = "Primary multimodal (2-photon + PA mechanics)  |  LOO R2=0.834"

    print(f"Model: {model_label}")
    model         = joblib.load(model_file)
    preprocessor  = joblib.load(prep_file)
    feature_names = joblib.load(feat_file)
    print(f"Features: {feature_names}")

    print(f"\nLoading: {data_path}  (sheet: {sheet})")
    df = load_and_prep(data_path, sheet=sheet)

    missing = [f for f in feature_names if f not in df.columns]
    if missing:
        raise ValueError(f"Missing features in input data: {missing}")

    extra_cols = [c for c in ["Age_num"] if c in df.columns]
    sub = df[feature_names + extra_cols].copy()

    n_before = len(sub)
    sub = sub.dropna(subset=feature_names)
    if n_before - len(sub) > 0:
        print(f"Warning: dropped {n_before - len(sub)} rows with missing values")

    X = sub[feature_names].values
    print(f"Predicting on n={len(X)} samples")

    predictions = model.predict(preprocess(X, preprocessor))

    results = sub.copy()
    results["Predicted_Age_mo"] = predictions
    for c in ["Sample Name", "Sex", "Group"]:
        if c in df.columns:
            results[c] = df.loc[sub.index, c]

    if "Age_num" in sub.columns:
        true_age  = sub["Age_num"].values
        residuals = predictions - true_age
        mae = np.mean(np.abs(residuals))
        r2  = 1 - np.sum((true_age - predictions)**2) / np.sum((true_age - true_age.mean())**2)
        results["True_Age_mo"] = true_age
        results["Residual_mo"] = residuals
        print(f"\nR2  = {r2:.3f}")
        print(f"MAE = {mae:.2f} mo")

    # Average across pressure measurements (same Sample Name within same Group)
    group_cols = [c for c in ["Group", "Sample Name"] if c in results.columns]
    if group_cols:
        agg = {c: "first" for c in ["Sex", "True_Age_mo"] if c in results.columns}
        agg["Predicted_Age_mo"] = "mean"
        if "Residual_mo" in results.columns:
            agg["Residual_mo"] = "mean"
        results_avg = results.groupby(group_cols, sort=False).agg(agg).reset_index()
    else:
        results_avg = results.copy()

    display_cols = [c for c in ["Group", "Sample Name", "Sex", "True_Age_mo",
                                 "Predicted_Age_mo", "Residual_mo"]
                    if c in results_avg.columns] or ["Predicted_Age_mo"]

    print("\nPredictions (averaged across pressures per sample):")
    print(results_avg[display_cols].to_string(index=False, float_format="%.2f"))

    if output_path:
        results_avg[display_cols].to_excel(output_path, index=False)
        print(f"\nSaved to: {output_path}")

    return results_avg


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("data", help="Input Excel file (.xlsx)")
    parser.add_argument("--lung", action="store_true",
                        help="Use lung model (2-photon + Rn)")
    parser.add_argument("--2photon", dest="twophoton", action="store_true",
                        help="Use 2-photon-only model")
    parser.add_argument("--sheet", default=0,
                        help="Sheet name or index (default: first sheet)")
    parser.add_argument("--output", "-o", default=None,
                        help="Output Excel path (.xlsx)")
    args = parser.parse_args()

    sheet = int(args.sheet) if str(args.sheet).isdigit() else args.sheet

    if args.twophoton:
        model_type = "2photon"
    elif args.lung:
        model_type = "lung"
    else:
        model_type = "multimodal"

    predict(args.data, model_type=model_type, sheet=sheet, output_path=args.output)