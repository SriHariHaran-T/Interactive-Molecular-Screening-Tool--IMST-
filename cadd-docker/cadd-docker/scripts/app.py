"""
app.py
------
A simple interactive app: pick a drug (or type your own SMILES),
and see:
  1. Whether it's likely to be an effective ORAL drug (absorbed
     by the human gut) — using the same Lipinski + Veber rules
     from admet.py
  2. Its interactive 3D structure

Run with:
    streamlit run app.py

This opens a browser tab automatically — it's a real local web app,
not a chart printed to the terminal.
"""

import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski


# ---------- Reuse the same ADMET logic from admet.py ----------
def compute_admet_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    descriptors = {
        "MolWeight": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "HBondDonors": Lipinski.NumHDonors(mol),
        "HBondAcceptors": Lipinski.NumHAcceptors(mol),
        "RotatableBonds": Descriptors.NumRotatableBonds(mol),
    }

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
    descriptors["LikelyOralDrug"] = violations <= 1

    tpsa_ok = descriptors["TPSA"] <= 140
    flexibility_ok = descriptors["RotatableBonds"] <= 10
    descriptors["LikelyGutAbsorbed"] = tpsa_ok and flexibility_ok

    return descriptors


def generate_3d_molblock(smiles):
    """Builds a 3D structure from SMILES and returns it as a
    'molblock' — a text format 3Dmol.js (the viewer library) can read."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return Chem.MolToMolBlock(mol)


def render_3d_viewer(molblock, height=400):
    """Embeds a 3Dmol.js viewer directly in the Streamlit app via HTML.
    3Dmol.js is a free JavaScript library for viewing molecules,
    loaded from a CDN — no installation needed."""
    html = f"""
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
    <div id="viewer" style="height: {height}px; width: 100%; position: relative;"></div>
    <script>
      let viewer = $3Dmol.createViewer("viewer", {{backgroundColor: "white"}});
      viewer.addModel(`{molblock}`, "mol");
      viewer.setStyle({{}}, {{stick: {{}}, sphere: {{scale: 0.3}}}});
      viewer.zoomTo();
      viewer.render();
    </script>
    """
    st.components.v1.html(html, height=height)


# ---------- The App ----------
st.set_page_config(page_title="Interactive Molecular Screening Tool (IMST)", layout="centered")
st.title("Interactive Molecular Screening Tool (IMST)")
st.write(
    "Pick a drug, or enter your own molecule, to check whether it's "
    "likely to be effectively absorbed as an oral drug and see its "
    "3D structure."
)

import os

# Load the starter drug list
current_dir = os.path.dirname(os.path.abspath(__file__))
csv_path = os.path.join(current_dir, "..", "data", "datasets", "drugs.csv")
drugs_df = pd.read_csv(csv_path)

mode = st.radio("Choose input method:", ["Pick from list", "Enter my own SMILES"])

if mode == "Pick from list":
    drug_name = st.selectbox("Select a drug:", drugs_df["name"].tolist())
    smiles = drugs_df.loc[drugs_df["name"] == drug_name, "smiles"].values[0]
else:
    smiles = st.text_input("Enter a SMILES string:", "CC(=O)OC1=CC=CC=C1C(=O)O")
    drug_name = "Custom molecule"

if smiles:
    result = compute_admet_descriptors(smiles)

    if result is None:
        st.error("Couldn't parse that SMILES string — check it's valid.")
    else:
        st.subheader(drug_name)

        col1, col2 = st.columns(2)
        with col1:
            if result["LikelyGutAbsorbed"]:
                st.success("✅ Likely absorbed by the gut")
            else:
                st.warning("⚠️ Likely POOR gut absorption")
        with col2:
            if result["LikelyOralDrug"]:
                st.success("✅ Fits general oral drug profile")
            else:
                st.warning("⚠️ Fails oral drug-likeness checks")

        st.write("**Computed properties:**")
        st.table(pd.DataFrame([result]).T.rename(columns={0: "Value"}))

        st.write("**3D Structure:**")
        try:
            molblock = generate_3d_molblock(smiles)
            render_3d_viewer(molblock)
        except Exception as e:
            st.error(f"Couldn't generate 3D structure: {e}")
