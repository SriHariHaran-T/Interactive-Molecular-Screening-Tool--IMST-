"""
qsar.py
-------
The core deliverable: train a model that predicts whether a molecule
is ACTIVE or INACTIVE against a chosen protein target.

Pipeline (matches the theory we walked through):
  1. Pull labeled bioactivity data from ChEMBL for one target
  2. Clean it: keep IC50 values, convert to pIC50, binarize into
     active (1) / inactive (0) at a threshold
  3. Featurize each molecule with RDKit Morgan fingerprints
  4. Split into train/test using SCAFFOLD splitting (not random!)
     so structurally similar molecules don't leak between sets
  5. Train a Random Forest classifier
  6. Evaluate with accuracy + ROC-AUC
  7. Save the trained model + results

Usage:
    python qsar.py --target_chembl_id CHEMBL203
    (CHEMBL203 = EGFR, just an example — swap for your chosen target)
"""

import argparse
import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, roc_auc_score
import pickle


# ---------- Step 1: Pull data from ChEMBL ----------
def fetch_chembl_data(target_chembl_id, max_records=2000):
    """Pulls IC50 bioactivity records for a given ChEMBL target ID.
    Requires internet access to ChEMBL's servers."""
    from chembl_webresource_client.new_client import new_client

    activity = new_client.activity
    results = activity.filter(
        target_chembl_id=target_chembl_id,
        standard_type="IC50"
    ).only(["canonical_smiles", "standard_value", "standard_units"])

    df = pd.DataFrame(results[:max_records])
    print(f"Pulled {len(df)} raw records from ChEMBL for {target_chembl_id}")
    return df


# ---------- Step 2: Clean + label data ----------
def clean_data(df, active_threshold_nm=1000):
    """Keeps valid nM IC50 values, converts to pIC50, and labels
    active (1) if IC50 < threshold, else inactive (0).
    Default threshold: 1000 nM (1 uM) is a common cutoff."""
    df = df.dropna(subset=["canonical_smiles", "standard_value"])
    df = df[df["standard_units"] == "nM"]
    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = df.dropna(subset=["standard_value"])
    df = df[df["standard_value"] > 0]

    # pIC50 = -log10(IC50 in molar) — standard transform, keeps values
    # in a more model-friendly range and larger = more potent
    df["pIC50"] = -np.log10(df["standard_value"] * 1e-9)

    df["active"] = (df["standard_value"] < active_threshold_nm).astype(int)

    print(f"After cleaning: {len(df)} usable records")
    print(f"Active: {df['active'].sum()}  Inactive: {(df['active'] == 0).sum()}")
    return df


# ---------- Step 3: Featurize molecules ----------
def smiles_to_fingerprint(smiles, n_bits=2048, radius=2):
    """Converts a SMILES string into a Morgan fingerprint —
    a fixed-length bit vector encoding local substructures."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    return np.array(fp)


def featurize(df):
    fingerprints = []
    valid_rows = []
    for idx, row in df.iterrows():
        fp = smiles_to_fingerprint(row["canonical_smiles"])
        if fp is not None:
            fingerprints.append(fp)
            valid_rows.append(idx)

    X = np.array(fingerprints)
    y = df.loc[valid_rows, "active"].values
    scaffolds = df.loc[valid_rows, "canonical_smiles"].apply(get_scaffold)
    print(f"Featurized {len(X)} molecules into {X.shape[1]}-bit fingerprints")
    return X, y, scaffolds.values


# ---------- Step 4: Scaffold split ----------
def get_scaffold(smiles):
    """Extracts the Murcko scaffold (core ring structure) of a molecule.
    We split train/test by scaffold so the model is tested on
    structurally NEW cores, not near-duplicates of training molecules."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    return Chem.MolToSmiles(scaffold)


def scaffold_split(X, y, scaffolds, test_fraction=0.2, seed=42):
    unique_scaffolds = np.array(pd.Series(scaffolds).unique(), dtype=object)
    rng = np.random.RandomState(seed)
    rng.shuffle(unique_scaffolds)

    n_test_scaffolds = int(len(unique_scaffolds) * test_fraction)
    test_scaffolds = set(unique_scaffolds[:n_test_scaffolds])

    test_mask = pd.Series(scaffolds).isin(test_scaffolds).values
    train_mask = ~test_mask

    return X[train_mask], X[test_mask], y[train_mask], y[test_mask]


# ---------- Step 5-7: Train, evaluate, save ----------
def train_and_evaluate(X_train, X_test, y_train, y_test, output_dir="results/qsar"):
    os.makedirs(output_dir, exist_ok=True)

    model = RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1)
    model.fit(X_train, y_train)

    preds = model.predict(X_test)
    probs = model.predict_proba(X_test)[:, 1]

    acc = accuracy_score(y_test, preds)
    try:
        auc = roc_auc_score(y_test, probs)
    except ValueError:
        auc = float("nan")  # happens if test set has only one class

    print(f"Test Accuracy: {acc:.3f}")
    print(f"Test ROC-AUC:  {auc:.3f}")

    model_path = os.path.join(output_dir, "qsar_model.pkl")
    with open(model_path, "wb") as f:
        pickle.dump(model, f)
    print(f"Model saved to: {model_path}")

    return model, acc, auc


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--target_chembl_id", required=True,
                         help="e.g. CHEMBL203 for EGFR")
    parser.add_argument("--threshold_nm", type=float, default=1000)
    args = parser.parse_args()

    raw = fetch_chembl_data(args.target_chembl_id)
    clean = clean_data(raw, active_threshold_nm=args.threshold_nm)
    clean.to_csv(f"data/datasets/{args.target_chembl_id}_clean.csv", index=False)

    X, y, scaffolds = featurize(clean)
    X_train, X_test, y_train, y_test = scaffold_split(X, y, scaffolds)

    print(f"Train set: {len(X_train)}  Test set: {len(X_test)}")
    train_and_evaluate(X_train, X_test, y_train, y_test)
