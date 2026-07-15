"""
Random Forest Model — Predicting Mouse Age from Vascular Parameters
Data: Aging_PA_03-20-2026.xlsx
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import (
    classification_report, confusion_matrix, ConfusionMatrixDisplay,
    mean_absolute_error, r2_score
)

# LOAD & CLEAN
df = pd.read_excel("Aging_PA_03-20-2026.xlsx", header=0)

# Drop "Average" summary rows and rows missing Age
df = df[~df["Sample Name"].str.contains("Average", na=False)]
df = df.dropna(subset=["Age"])

# Drop COLLAGEN BUNDLE WIDTH (entirely empty)
df = df.drop(columns=["COLLAGEN BUNDLE WIDTH"])

# Parse age in months from strings like "3mo", "6mo"
df["Age_months"] = df["Age"].str.replace("mo", "", regex=False).astype(int)

# Encode Sex: Male=0, Female=1
df["Sex_encoded"] = (df["Sex"].str.strip().str.lower() == "female").astype(int)

# Define features
FEATURE_COLS = [
    "Sex_encoded",
    "ADVENTITIA THICKNESS",
    "MEDIA THICKNESS",
    "ORIENTATION-θ",
    "ORIENTATION-κ",
    "ADVENTITIAL CELL COUNT",
    "SMC CELL COUNT",
    "ENDOTHELIAL CELL COUNT",
    "COLLAGEN STRAIGHTNESS",
]

# Drop rows missing any feature
df_model = df.dropna(subset=FEATURE_COLS).copy()
print(f"Samples used for modelling: {len(df_model)}")
print(f"Age distribution:\n{df_model['Age_months'].value_counts().sort_index()}\n")

X = df_model[FEATURE_COLS].values
y_reg = df_model["Age_months"].values
y_cls = df_model["Age_months"].astype(str) + "mo"

# TRAIN / TEST SPLIT
X_train, X_test, y_train_r, y_test_r, y_train_c, y_test_c = train_test_split(
    X, y_reg, y_cls, test_size=0.2, random_state=42, stratify=y_cls
)

# CLASSIFIER (predicts age group)
clf = RandomForestClassifier(
    n_estimators=300,
    max_depth=None,
    min_samples_leaf=2,
    random_state=42,
    class_weight="balanced",
)
clf.fit(X_train, y_train_c)
y_pred_c = clf.predict(X_test)

# Cross-validation accuracy
cv_acc = cross_val_score(clf, X, y_cls, cv=StratifiedKFold(5, shuffle=True, random_state=42), scoring="accuracy")
print(f"Classifier CV Accuracy: {cv_acc.mean():.3f} ± {cv_acc.std():.3f}")
print(classification_report(y_test_c, y_pred_c))

# REGRESSOR (predicts exact age in months)
reg = RandomForestRegressor(
    n_estimators=300,
    max_depth=None,
    min_samples_leaf=2,
    random_state=42,
)
reg.fit(X_train, y_train_r)
y_pred_r = reg.predict(X_test)

mae  = mean_absolute_error(y_test_r, y_pred_r)
r2   = r2_score(y_test_r, y_pred_r)
cv_r2 = cross_val_score(reg, X, y_reg, cv=5, scoring="r2")
print(f"\nRegressor  MAE : {mae:.2f} months")
print(f"Regressor  R²  : {r2:.3f}")
print(f"Regressor CV R²: {cv_r2.mean():.3f} ± {cv_r2.std():.3f}\n")

# PLOTS
AGE_ORDER = ["3mo", "6mo", "12mo", "18mo", "24mo", "27mo"]
palette   = sns.color_palette("viridis", len(AGE_ORDER))

fig = plt.figure(figsize=(18, 14))
fig.patch.set_facecolor("#f9f9f9")
gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.35)

# Feature Importance
ax1 = fig.add_subplot(gs[0, :2])
importances = pd.Series(reg.feature_importances_, index=FEATURE_COLS).sort_values()
colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(importances)))
importances.plot(kind="barh", ax=ax1, color=colors, edgecolor="white", linewidth=0.5)
ax1.set_title("Feature Importance (Random Forest Regressor)", fontsize=13, fontweight="bold", pad=10)
ax1.set_xlabel("Mean Decrease in Impurity", fontsize=10)
ax1.tick_params(axis="y", labelsize=9)
for bar, val in zip(ax1.patches, importances.values):
    ax1.text(bar.get_width() + 0.002, bar.get_y() + bar.get_height() / 2,
             f"{val:.3f}", va="center", fontsize=8)
ax1.set_facecolor("#f9f9f9")

# CV Accuracy Distribution
ax2 = fig.add_subplot(gs[0, 2])
ax2.bar(range(1, 6), cv_acc, color=palette[2], edgecolor="white", width=0.6)
ax2.axhline(cv_acc.mean(), color="#e74c3c", linewidth=1.5, linestyle="--", label=f"Mean = {cv_acc.mean():.2f}")
ax2.set_title("5-Fold CV Accuracy\n(Classifier)", fontsize=13, fontweight="bold", pad=10)
ax2.set_xlabel("Fold", fontsize=10)
ax2.set_ylabel("Accuracy", fontsize=10)
ax2.set_ylim(0, 1.05)
ax2.legend(fontsize=9)
ax2.set_facecolor("#f9f9f9")

# Confusion Matrix
ax3 = fig.add_subplot(gs[1, 0])
present = [a for a in AGE_ORDER if a in np.unique(np.concatenate([y_test_c, y_pred_c]))]
cm = confusion_matrix(y_test_c, y_pred_c, labels=present)
disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=present)
disp.plot(ax=ax3, cmap="Blues", colorbar=False)
ax3.set_title("Confusion Matrix\n(Classifier — Test Set)", fontsize=12, fontweight="bold", pad=8)
ax3.tick_params(axis="both", labelsize=8)

# Predicted vs Actual (Regressor)
ax4 = fig.add_subplot(gs[1, 1])
for i, age in enumerate([3, 6, 12, 18, 24, 27]):
    mask = y_test_r == age
    if mask.sum():
        ax4.scatter(y_test_r[mask], y_pred_r[mask], color=palette[i],
                    label=f"{age}mo", s=60, edgecolors="white", linewidths=0.5, alpha=0.9)
lims = [0, 30]
ax4.plot(lims, lims, "k--", linewidth=1, alpha=0.5, label="Perfect prediction")
ax4.set_xlabel("Actual Age (months)", fontsize=10)
ax4.set_ylabel("Predicted Age (months)", fontsize=10)
ax4.set_title(f"Predicted vs Actual Age\n(Regressor — R²={r2:.2f}, MAE={mae:.1f} mo)", fontsize=12, fontweight="bold", pad=8)
ax4.legend(fontsize=8, ncol=2)
ax4.set_facecolor("#f9f9f9")

# Residuals by Age Group
ax5 = fig.add_subplot(gs[1, 2])
residuals = y_pred_r - y_test_r
res_df = pd.DataFrame({"Age": y_test_r, "Residual": residuals})
for i, age in enumerate([3, 6, 12, 18, 24, 27]):
    subset = res_df[res_df["Age"] == age]["Residual"]
    if len(subset):
        ax5.scatter([age] * len(subset), subset, color=palette[i], s=55,
                    edgecolors="white", linewidths=0.5, alpha=0.85)
ax5.axhline(0, color="#e74c3c", linewidth=1.5, linestyle="--")
ax5.set_xlabel("Actual Age (months)", fontsize=10)
ax5.set_ylabel("Residual (Predicted − Actual)", fontsize=10)
ax5.set_title("Residuals by Age Group\n(Regressor)", fontsize=12, fontweight="bold", pad=8)
ax5.set_xticks([3, 6, 12, 18, 24, 27])
ax5.set_facecolor("#f9f9f9")

fig.suptitle("Random Forest — Mouse Vascular Aging Model", fontsize=16, fontweight="bold", y=1.01)
plt.savefig("random_forest_results.png", dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.close()
print("Saved: random_forest_results.png")

# SAVE PREDICTIONS
test_idx = df_model.index[
    df_model.index.isin(
        df_model.iloc[
            np.where(np.isin(np.arange(len(df_model)),
                             np.where(np.isin(X, X_test).all(axis=1))[0]))[0]
        ].index
    )
]

results = pd.DataFrame({
    "Sample Name": df_model.loc[df_model.index[
        [i for i in range(len(X)) if list(X[i]) in [list(r) for r in X_test]]
    ], "Sample Name"].values if False else df_model["Sample Name"].iloc[-len(y_test_r):].values,
    "Actual_Age_months": y_test_r,
    "Predicted_Age_months_regressor": np.round(y_pred_r, 1),
    "Predicted_Age_group_classifier": y_pred_c,
})
results.to_csv("predictions_test_set.csv", index=False)
print("Saved: predictions_test_set.csv")
print("\nDone.")
