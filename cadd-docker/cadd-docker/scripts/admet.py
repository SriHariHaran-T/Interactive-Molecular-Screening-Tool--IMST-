"""
admet.py
--------
Computes free, instant ADMET-relevant properties for a molecule
using RDKit descriptors, and checks two rule-based heuristics:

  - Lipinski's Rule of Five  -> general oral drug-likeness
  - Veber's Rule             -> specifically: will it likely be
                                 absorbed through the human gut?
                                 (this is the exact question the
                                 professor asked about)

This is the "free tier" of ADMET prediction — pure rules, no
trained model needed yet. For toxicity endpoints specifically
(Tox21, ClinTox), you'd extend this with a DeepChem pretrained
model — noted at the bottom.

Usage (single molecule, typed directly):
    python admet.py "CC(=O)OC1=CC=CC=C1C(=O)O"   (aspirin)

Usage (many molecules at once, from a CSV file):
    python admet.py --csv data/datasets/drugs.csv
    (the CSV needs a column named "smiles" — a "name" column is
    optional but makes the results easier to read)
"""

import sys
import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski


def compute_admet_descriptors(smiles, name=None):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")

    descriptors = {
        "Name": name if name else "",
        "SMILES": smiles,
        "MolWeight": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),          # topological polar surface area
                                                  # (relates to absorption/BBB penetration)
        "HBondDonors": Lipinski.NumHDonors(mol),
        "HBondAcceptors": Lipinski.NumHAcceptors(mol),
        "RotatableBonds": Descriptors.NumRotatableBonds(mol),
    }

    # ---- Lipinski's Rule of Five ----
    # General "does this look like a normal oral drug" heuristic
    violations = 0
    if descriptors["MolWeight"] > 500:
        violations += 1
    if descriptors["LogP"] > 5:
        violations += 1
    if descriptors["HBondDonors"] > 5:
        violations += 1
    if descriptors["HBondAcceptors"] > 10:
        violations += 1

    descriptors["LipinskiViolations"] = violations
    descriptors["LikelyOralDrug"] = violations <= 1  # 1 violation is generally tolerated

    # ---- Veber's Rule ----
    # Specifically predicts oral bioavailability / gut absorption.
    # A molecule needs BOTH conditions true to be considered
    # likely well-absorbed through the intestine:
    #   - TPSA <= 140 (too polar = struggles to cross the gut wall)
    #   - Rotatable bonds <= 10 (too flexible = poor absorption)
    tpsa_ok = descriptors["TPSA"] <= 140
    flexibility_ok = descriptors["RotatableBonds"] <= 10

    descriptors["PassesVeberRule"] = tpsa_ok and flexibility_ok
    descriptors["LikelyGutAbsorbed"] = descriptors["PassesVeberRule"]

    return descriptors


def run_on_list(smiles_list, names=None, output_path="results/admet/admet_results.csv"):
    if names is None:
        names = [None] * len(smiles_list)

    rows = []
    for smiles, name in zip(smiles_list, names):
        try:
            rows.append(compute_admet_descriptors(smiles, name))
        except ValueError as e:
            print(f"Skipping invalid molecule: {e}")

    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    print(df.to_string(index=False))
    print(f"\nSaved to: {output_path}")
    return df


def run_on_csv(csv_path, output_path="results/admet/admet_results.csv"):
    input_df = pd.read_csv(csv_path)
    if "smiles" not in input_df.columns:
        raise ValueError('CSV must have a column named "smiles"')

    smiles_list = input_df["smiles"].tolist()
    names = input_df["name"].tolist() if "name" in input_df.columns else None
    return run_on_list(smiles_list, names=names, output_path=output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--csv", help="Path to a CSV file with a 'smiles' column")
    parser.add_argument("smiles_args", nargs="*", help="One or more SMILES strings")
    args = parser.parse_args()

    if args.csv:
        run_on_csv(args.csv)
    elif args.smiles_args:
        run_on_list(args.smiles_args)
    else:
        print('Usage: python admet.py "<SMILES>" ["<SMILES2>" ...]')
        print('   or: python admet.py --csv data/datasets/drugs.csv')
        sys.exit(1)

    # ---- Extending to real toxicity prediction (optional, later) ----
    # import deepchem as dc
    # tasks, datasets, transformers = dc.molnet.load_tox21(featurizer='ECFP')
    # model = dc.models.MultitaskClassifier(n_tasks=len(tasks), n_features=1024)
    # model.fit(datasets[0])  # train[0] is the training set
    # -> then model.predict() on new molecules for tox endpoints
