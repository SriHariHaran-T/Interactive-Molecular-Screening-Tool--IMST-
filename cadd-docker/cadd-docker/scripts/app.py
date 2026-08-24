"""
app.py
------
Interactive Molecular Screening Tool (IMST)

Three tabs:
  1. ADMET Profiling  - oral drug-likeness + 3D structure viewer
  2. Batch Docking    - dock library compounds against HIV-1 Protease (1HVR)
  3. Live Docking     - dock a user-supplied SMILES against HIV-1 Protease

Run with:
    python -m streamlit run app.py
"""

import os
import sys
import json
import time
import hashlib
import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Lipinski

# ---------------------------------------------------------------------------
# Project paths (resolved relative to this file; no hard-coded user paths)
# ---------------------------------------------------------------------------
SCRIPT_DIR     = os.path.dirname(os.path.abspath(__file__))
DATA_DIR       = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "data"))
RESULTS_DIR    = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "results"))
RECEPTOR_PDB   = os.path.join(DATA_DIR, "proteins", "1HVR_clean.pdb")
RECEPTOR_PDBQT = os.path.join(DATA_DIR, "proteins", "1HVR_clean.pdbqt")
CSV_PATH       = os.path.join(DATA_DIR, "datasets", "drugs.csv")
LIGANDS_DIR    = os.path.join(DATA_DIR, "ligands")
DOCKING_DIR    = os.path.join(RESULTS_DIR, "docking")

# HIV-1 Protease 1HVR docking box
# Center derived from the co-crystallised XK2 inhibitor in the PDB structure.
HIV_CENTER = (-9.19, 15.91, 27.95)
HIV_BOX    = (20, 20, 20)

if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from prepare_ligand import prepare_ligand
from docking        import run_docking
from contacts       import analyze_contacts


# ---------------------------------------------------------------------------
# Preset molecule library
# Each entry: (display_name, smiles, brief_note)
# SMILES are verified canonical strings from PubChem.
# ---------------------------------------------------------------------------
PRESET_MOLECULES = [
    ("Aspirin",     "CC(=O)OC1=CC=CC=C1C(=O)O",
     "Non-steroidal anti-inflammatory and analgesic agent."),
    ("Caffeine",    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
     "Central nervous system stimulant found in coffee and tea."),
    ("Ibuprofen",   "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
     "Non-steroidal anti-inflammatory drug (NSAID)."),
    ("Paracetamol", "CC(=O)Nc1ccc(O)cc1",
     "Analgesic and antipyretic; also known as acetaminophen."),
    ("Nicotine",    "CN1CCCC1c1cccnc1",
     "Alkaloid found in tobacco; acts on nicotinic acetylcholine receptors."),
]
PRESET_NAMES  = [m[0] for m in PRESET_MOLECULES]
PRESET_LOOKUP = {m[0]: (m[1], m[2]) for m in PRESET_MOLECULES}


# ---------------------------------------------------------------------------
# ADMET calculations  (logic unchanged from original)
# ---------------------------------------------------------------------------

def compute_admet_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    desc = {
        "Molecular weight (Da)": Descriptors.MolWt(mol),
        "LogP":                  Descriptors.MolLogP(mol),
        "TPSA (A^2)":            Descriptors.TPSA(mol),
        "H-bond donors":         Lipinski.NumHDonors(mol),
        "H-bond acceptors":      Lipinski.NumHAcceptors(mol),
        "Rotatable bonds":       Descriptors.NumRotatableBonds(mol),
    }
    violations = sum([
        Descriptors.MolWt(mol)           > 500,
        Descriptors.MolLogP(mol)         > 5,
        Lipinski.NumHDonors(mol)         > 5,
        Lipinski.NumHAcceptors(mol)      > 10,
    ])
    desc["Lipinski violations"]          = violations
    desc["Passes Lipinski heuristic"]    = violations <= 1
    desc["Passes Veber heuristic"]       = (
        Descriptors.TPSA(mol) <= 140 and
        Descriptors.NumRotatableBonds(mol) <= 10
    )
    return desc


def generate_3d_molblock(smiles):
    """Generate a 3D conformer from SMILES and return as a molblock string."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    return Chem.MolToMolBlock(mol)


# ---------------------------------------------------------------------------
# 3Dmol.js rendering helpers
# ---------------------------------------------------------------------------

def _read_file_text(path):
    with open(path, "r") as f:
        return f.read()


def render_3d_viewer(molblock, height=400):
    """Render a single molecule in 3Dmol.js (used in ADMET tab)."""
    safe = molblock.replace("`", "'")
    html = f"""
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
    <div id="viewer_admet" style="height:{height}px;width:100%;position:relative;
         border:1px solid #30363d;border-radius:6px;"></div>
    <script>
      let v = $3Dmol.createViewer("viewer_admet", {{backgroundColor:"#0d1117"}});
      v.addModel(`{safe}`, "mol");
      v.setStyle({{}}, {{stick:{{}}, sphere:{{scale:0.3}}}});
      v.zoomTo();
      v.render();
    </script>
    """
    st.components.v1.html(html, height=height)


def render_docking_viewer(protein_pdb_path, ligand_pdbqt_path, contacts, height=520):
    """
    Render the protein and actual Vina docked ligand pose in 3Dmol.js.
    Contacts are drawn as distance lines.
    """
    protein_text = _read_file_text(protein_pdb_path).replace("`", "'").replace("\\", "\\\\")

    # Extract only the first MODEL block from the docked PDBQT (best pose)
    ligand_lines = []
    with open(ligand_pdbqt_path) as f:
        for line in f:
            if line.startswith("ENDMDL"):
                break
            if line.startswith("ATOM") or line.startswith("HETATM"):
                ligand_lines.append(line[:66])
    ligand_text = "".join(ligand_lines).replace("`", "'").replace("\\", "\\\\")

    # Build contact coordinate data for JavaScript
    hbond_pairs, contact_pairs, seen = [], [], set()
    for c in contacts:
        key = (tuple(c["ligand_coord"]), tuple(c["protein_coord"]))
        if key in seen:
            continue
        seen.add(key)
        lc, pc = c["ligand_coord"], c["protein_coord"]
        entry = {"start": lc, "end": pc, "dist": c["distance"], "res": c["protein_res"]}
        if c["type"] == "hbond":
            hbond_pairs.append(entry)
        else:
            contact_pairs.append(entry)

    hbond_js   = json.dumps(hbond_pairs[:20])
    contact_js = json.dumps(contact_pairs[:30])
    viewer_id  = f"dockview_{int(time.time() * 1000) % 99999}"

    html = f"""
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
    <div id="{viewer_id}" style="height:{height}px;width:100%;position:relative;
         border:1px solid #30363d;border-radius:6px;"></div>
    <script>
    (function(){{
      let viewer = $3Dmol.createViewer("{viewer_id}", {{backgroundColor:"#0d1117"}});

      // Protein - spectrum cartoon
      viewer.addModel(`{protein_text}`, "pdb");
      viewer.setStyle({{model:0}}, {{cartoon:{{color:"spectrum", opacity:0.85}}}});

      // Docked ligand - actual Vina pose (green sticks)
      viewer.addModel(`{ligand_text}`, "pdb");
      viewer.setStyle({{model:1}}, {{
        stick:  {{radius:0.2, colorscheme:"greenCarbon"}},
        sphere: {{scale:0.25, colorscheme:"greenCarbon"}}
      }});

      // Putative H-bond lines (gold)
      let hbonds = {hbond_js};
      hbonds.forEach(function(h) {{
        viewer.addCylinder({{
          start:   {{x:h.start[0], y:h.start[1], z:h.start[2]}},
          end:     {{x:h.end[0],   y:h.end[1],   z:h.end[2]}},
          radius:  0.06, color:"#FFD700",
          fromCap: 1, toCap:1, opacity:0.85
        }});
        viewer.addLabel(h.dist.toFixed(2)+" A", {{
          position: {{
            x:(h.start[0]+h.end[0])/2,
            y:(h.start[1]+h.end[1])/2,
            z:(h.start[2]+h.end[2])/2
          }},
          fontSize:10, fontColor:"#FFD700",
          backgroundColor:"rgba(0,0,0,0.5)"
        }});
      }});

      // Hydrophobic contact lines (light blue, thin)
      let contacts = {contact_js};
      contacts.forEach(function(c) {{
        viewer.addCylinder({{
          start:  {{x:c.start[0], y:c.start[1], z:c.start[2]}},
          end:    {{x:c.end[0],   y:c.end[1],   z:c.end[2]}},
          radius: 0.03, color:"#00BFFF", opacity:0.45
        }});
      }});

      viewer.zoomTo({{model:1}});
      viewer.render();
    }})();
    </script>
    """
    st.components.v1.html(html, height=height + 20)


# ---------------------------------------------------------------------------
# Docking pipeline helpers
# ---------------------------------------------------------------------------

def _safe_slug(name):
    return "".join(c if c.isalnum() else "_" for c in name).strip("_")


def run_single_drug(name, smiles):
    """
    Prepare one drug from SMILES and dock against 1HVR.
    Returns (name, smiles, score, docked_path, error_msg).
    One compound failure does not stop the batch.
    """
    slug = _safe_slug(name)
    try:
        ligand_pdbqt = prepare_ligand(smiles, slug)
    except Exception as e:
        return name, smiles, None, None, f"Ligand preparation failed: {e}"

    docked_path, score = run_docking(
        receptor=RECEPTOR_PDBQT, ligand=ligand_pdbqt,
        center=HIV_CENTER, box_size=HIV_BOX, exhaustiveness=8
    )
    if docked_path is None:
        return name, smiles, None, None, "AutoDock Vina produced no result"
    return name, smiles, score, docked_path, None


# ---------------------------------------------------------------------------
# Page configuration
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="IMST - Interactive Molecular Screening Tool",
    page_icon="🧬",
    layout="wide",
)

st.markdown("""
<style>
  .stApp { background:#0d1117; color:#e6edf3; }
  .stTabs [data-baseweb="tab-list"] { gap:6px; }
  .stTabs [data-baseweb="tab"] {
      border-radius:5px 5px 0 0; padding:7px 16px;
      background:#161b22; color:#8b949e; font-size:0.9rem;
  }
  .stTabs [aria-selected="true"] {
      background:#1f6feb !important; color:#ffffff !important;
  }
  div[data-testid="stMetricValue"] { font-size:1.4rem; }
  .guide-term  { font-weight:600; color:#79c0ff; }
  .disclaimer  { font-size:0.82rem; color:#8b949e; border-left:3px solid #30363d;
                 padding-left:10px; margin-top:8px; }
</style>
""", unsafe_allow_html=True)


st.title("Interactive Molecular Screening Tool (IMST)")
st.caption(
    "Computational ADMET profiling and molecular docking against HIV-1 Protease (1HVR). "
    "All results are computational predictions only."
)

# ---------------------------------------------------------------------------
# User Guide (collapsible, near top)
# ---------------------------------------------------------------------------
with st.expander("How to Use this Tool", expanded=False):
    st.markdown("""
**Overview**

This tool lets you explore how a small molecule (ligand) might interact with a
protein target (HIV-1 Protease) using computational methods.
It is intended as an educational and research demonstration.
All outputs are computational predictions, not experimental measurements.

---

**Step-by-step workflow**

1. **Select or enter a molecule.**
   Use the preset library or type a SMILES string.

2. **Review ADMET properties** (ADMET Profiling tab).
   The tool computes structural descriptors and checks two rule-based heuristics:
   Lipinski's Rule of Five (predicts oral drug-likeness) and Veber's rule
   (predicts intestinal absorption). These are screening filters, not safety tests.

3. **Run docking** (Live Custom Docking or Batch Docking tab).
   The ligand is prepared as a 3D structure, converted to the PDBQT format
   required by AutoDock Vina, and docked against the HIV-1 Protease receptor.

4. **Examine the predicted pose** in the 3D viewer.
   The viewer shows the receptor (protein) and the predicted docked pose of the
   ligand. Rotate and zoom with the mouse.

5. **Review predicted contacts.**
   The tool identifies close-distance pairs between ligand and protein atoms.
   Putative hydrogen bonds and hydrophobic contacts are highlighted.

---

**Key terms**

<span class="guide-term">SMILES</span> - Simplified Molecular Input Line Entry System.
A text notation that encodes a molecular structure, e.g. *CC(=O)OC1=CC=CC=C1C(=O)O* for aspirin.

<span class="guide-term">Ligand</span> - The small molecule being tested. In drug discovery, ligands
are candidate compounds that may bind to a protein target.

<span class="guide-term">Protein / target</span> - The biological macromolecule the ligand is docked
against. Here the target is HIV-1 Protease, a key enzyme in HIV replication, supplied as a
3D structural model (PDB ID: 1HVR).

<span class="guide-term">Binding pocket</span> - The cavity on the protein surface where a ligand
is predicted (or observed experimentally) to bind.

<span class="guide-term">PDBQT</span> - A molecular file format used by AutoDock Vina.
It extends the standard PDB format with partial atomic charges and atom types.

<span class="guide-term">AutoDock Vina</span> - An open-source molecular docking program that
searches for favorable binding poses of a ligand within a defined region of a protein structure.

<span class="guide-term">Docking affinity (kcal/mol)</span> - A predicted binding energy score.
More negative values indicate a more favorable predicted pose in this computational model.
This is NOT equivalent to experimentally measured binding affinity.

<span class="guide-term">Molecular contacts</span> - Atom pairs from the ligand and protein that
are within a defined distance threshold. These suggest possible physical interactions but
are not experimentally confirmed.

<span class="guide-term">Putative hydrogen bond</span> - A predicted hydrogen-bond interaction
inferred from distance (less than 3.5 A) and atom type (N, O, or S pairs).
These are computational approximations.

---

<div class="disclaimer">
All ADMET results, docking scores, binding poses, and contact predictions produced by this
tool are computational estimates based on structural heuristics and force-field scoring.
They do not constitute experimental evidence, safety data, toxicity assessments, clinical
evaluation, or regulatory approval of any kind. Do not make clinical, pharmaceutical, or
safety decisions based on these results.
</div>
""", unsafe_allow_html=True)


drugs_df = pd.read_csv(CSV_PATH)

tab1, tab2, tab3 = st.tabs([
    "ADMET Profiling",
    "Batch Target Docking",
    "Live Custom Docking",
])


# ===========================================================================
# TAB 1 - ADMET Profiling
# ===========================================================================
with tab1:
    st.subheader("ADMET Profiling")
    st.write(
        "Compute structural descriptors and apply rule-based drug-likeness heuristics "
        "for any small molecule. Enter a SMILES string, select a preset compound, "
        "or choose from the drug library."
    )

    st.markdown(
        '<div class="disclaimer">Results are rule-based computational heuristics '
        "(Lipinski Rule of Five, Veber's rule). They do not establish safety, "
        "efficacy, toxicity, clinical suitability, or regulatory approval.</div>",
        unsafe_allow_html=True
    )
    st.markdown("")

    # Input method
    input_mode = st.radio(
        "Input method:",
        ["Preset compounds", "Drug library", "Enter SMILES"],
        horizontal=True,
        key="admet_mode",
    )

    smiles       = ""
    drug_name    = "Custom molecule"
    preset_note  = ""

    if input_mode == "Preset compounds":
        selected_preset = st.selectbox(
            "Select a compound:", PRESET_NAMES, key="admet_preset"
        )
        p_smiles, p_note = PRESET_LOOKUP[selected_preset]
        st.markdown(f"**Name:** {selected_preset}")
        st.code(p_smiles, language=None)
        st.caption(p_note)
        st.caption(
            "Note: listing a compound here does not imply it is safe, "
            "suitable as a drug, or recommended for any use."
        )
        smiles    = p_smiles
        drug_name = selected_preset

    elif input_mode == "Drug library":
        selected_drug = st.selectbox(
            "Select a drug:", drugs_df["name"].tolist(), key="admet_drug"
        )
        smiles    = drugs_df.loc[drugs_df["name"] == selected_drug, "smiles"].values[0]
        drug_name = selected_drug

    else:  # Enter SMILES
        smiles    = st.text_input(
            "SMILES string:",
            value="CC(=O)OC1=CC=CC=C1C(=O)O",
            key="admet_smiles",
        )
        drug_name = "Custom molecule"

    if smiles:
        result = compute_admet_descriptors(smiles)
        if result is None:
            st.error("Cannot parse this SMILES string. Please check that it is valid.")
        else:
            st.markdown("---")
            st.subheader(drug_name)

            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Mol. Weight",    f"{result['Molecular weight (Da)']:.1f} Da")
            c2.metric("LogP",           f"{result['LogP']:.2f}")
            c3.metric("TPSA",           f"{result['TPSA (A^2)']:.1f} A\u00b2")
            c4.metric("H-bond donors",  result["H-bond donors"])

            col_l, col_v = st.columns(2)
            with col_l:
                if result["Passes Lipinski heuristic"]:
                    st.success("Passes Lipinski drug-likeness heuristic (Rule of Five)")
                else:
                    st.warning(
                        f"Fails Lipinski heuristic "
                        f"({result['Lipinski violations']} violation(s))"
                    )
            with col_v:
                if result["Passes Veber heuristic"]:
                    st.success("Passes Veber oral absorption heuristic")
                else:
                    st.warning("Fails Veber oral absorption heuristic")

            with st.expander("Full descriptor table"):
                st.table(
                    pd.DataFrame.from_dict(result, orient="index", columns=["Value"])
                )

            st.markdown("**3D Structure (energy-minimised conformer)**")
            try:
                molblock = generate_3d_molblock(smiles)
                render_3d_viewer(molblock)
            except Exception as e:
                st.error(f"Could not generate 3D structure: {e}")


# ===========================================================================
# TAB 2 - Batch Target Docking
# ===========================================================================
with tab2:
    st.subheader("Batch Docking - HIV-1 Protease (1HVR)")
    st.write(
        "Dock all compounds in the drug library against HIV-1 Protease "
        "and rank them by predicted binding affinity. "
        "Each compound is independently prepared and docked using AutoDock Vina. "
        "Estimated runtime: 20-90 seconds total."
    )

    st.markdown(
        '<div class="disclaimer">Predicted affinities are AutoDock Vina scores. '
        "A more negative score indicates a more favorable predicted pose "
        "within this model. These scores are not equivalent to experimental "
        "binding affinities and do not establish drug efficacy.</div>",
        unsafe_allow_html=True
    )
    st.markdown("")

    col_params, col_run = st.columns([3, 1])
    with col_params:
        st.markdown(
            f"**Target:** HIV-1 Protease (PDB: 1HVR)  \n"
            f"**Box center:** {HIV_CENTER}  \n"
            f"**Box size:** {HIV_BOX} A  \n"
            f"**Exhaustiveness:** 8"
        )
    with col_run:
        run_batch = st.button(
            "Run Batch Docking", type="primary",
            key="run_batch", use_container_width=True
        )

    if run_batch:
        if not os.path.exists(RECEPTOR_PDBQT):
            st.error(
                f"Receptor PDBQT not found at: `{RECEPTOR_PDBQT}`\n"
                "Run `prepare_1hvr.py` first to generate the receptor file."
            )
        else:
            results_rows = []
            prog        = st.progress(0, text="Starting...")
            status_box  = st.empty()
            n           = len(drugs_df)

            for i, row in drugs_df.iterrows():
                name, smi = row["name"], row["smiles"]
                status_box.write(f"Docking {name} ({i+1}/{n})...")
                name_r, smi_r, score, docked, err = run_single_drug(name, smi)
                results_rows.append({
                    "Compound":                     name_r,
                    "SMILES":                       smi_r[:40] + "..." if len(smi_r) > 40 else smi_r,
                    "Predicted affinity (kcal/mol)": score if score is not None else "n/a",
                    "Status":                       "OK" if err is None else f"Failed: {err}",
                })
                label = f"{name_r}: {score:.3f} kcal/mol" if score is not None else name_r
                prog.progress((i + 1) / n, text=label)

            status_box.empty()
            prog.empty()

            df_res  = pd.DataFrame(results_rows)
            ok_mask = df_res["Predicted affinity (kcal/mol)"] != "n/a"
            df_ok   = df_res[ok_mask].copy()
            df_fail = df_res[~ok_mask].copy()
            df_ok["Predicted affinity (kcal/mol)"] = df_ok["Predicted affinity (kcal/mol)"].astype(float)
            df_ok = df_ok.sort_values("Predicted affinity (kcal/mol)")

            st.success(f"Batch docking complete. {len(df_ok)} of {n} compounds succeeded.")
            st.subheader("Ranked by Predicted Binding Affinity")
            st.dataframe(
                df_ok.reset_index(drop=True),
                use_container_width=True,
                column_config={
                    "Predicted affinity (kcal/mol)": st.column_config.NumberColumn(format="%.3f")
                },
            )
            if len(df_fail):
                with st.expander(f"{len(df_fail)} compound(s) failed"):
                    st.dataframe(df_fail.reset_index(drop=True), use_container_width=True)


# ===========================================================================
# TAB 3 - Live Custom Docking
# ===========================================================================
with tab3:
    st.subheader("Live Custom Docking - HIV-1 Protease (1HVR)")
    st.write(
        "Enter a SMILES string to perform molecular docking against HIV-1 Protease. "
        "The tool prepares a 3D ligand structure, runs AutoDock Vina, and displays "
        "the predicted docked pose with contact analysis."
    )

    st.markdown(
        '<div class="disclaimer">The predicted affinity and docked pose are '
        "AutoDock Vina computational outputs. They are not experimental measurements "
        "and do not imply drug activity, safety, or suitability.</div>",
        unsafe_allow_html=True
    )
    st.markdown("")

    st.markdown(
        f"**Target:** HIV-1 Protease (PDB: 1HVR)  \n"
        f"**Box center:** {HIV_CENTER} | **Box size:** {HIV_BOX} A | **Exhaustiveness:** 8"
    )
    st.markdown("")

    # Molecule input: preset or custom SMILES
    input_src = st.radio(
        "Molecule source:",
        ["Preset compound", "Enter SMILES manually"],
        horizontal=True,
        key="live_src",
    )

    if input_src == "Preset compound":
        preset_sel  = st.selectbox("Select compound:", PRESET_NAMES, key="live_preset")
        live_smiles_default, live_note = PRESET_LOOKUP[preset_sel]
        st.markdown(f"**Name:** {preset_sel}")
        st.code(live_smiles_default, language=None)
        st.caption(live_note)
        smiles_input = live_smiles_default
    else:
        smiles_input = st.text_input(
            "SMILES string:",
            value="CC(=O)OC1=CC=CC=C1C(=O)O",
            placeholder="e.g. CC(=O)OC1=CC=CC=C1C(=O)O",
            key="live_smiles",
        )

    dock_btn = st.button("Dock and Analyze", type="primary", key="dock_live")

    if dock_btn and smiles_input.strip():
        smiles_clean = smiles_input.strip()

        if Chem.MolFromSmiles(smiles_clean) is None:
            st.error("Invalid SMILES string. Please correct it and try again.")
            st.stop()

        if not os.path.exists(RECEPTOR_PDBQT):
            st.error(
                f"Receptor file not found: `{RECEPTOR_PDBQT}`\n"
                "Run `prepare_1hvr.py` to generate the receptor PDBQT."
            )
            st.stop()

        slug = "live_" + hashlib.md5(smiles_clean.encode()).hexdigest()[:8]

        # Step 1: Ligand preparation
        with st.spinner("Step 1 of 3 - Preparing ligand (RDKit + Open Babel)..."):
            try:
                ligand_pdbqt = prepare_ligand(smiles_clean, slug, output_dir=LIGANDS_DIR)
                st.success(f"Ligand PDBQT ready: {os.path.basename(ligand_pdbqt)}")
            except Exception as e:
                st.error(f"Ligand preparation failed:\n```\n{e}\n```")
                st.stop()

        # Step 2: Docking
        with st.spinner("Step 2 of 3 - Running AutoDock Vina (exhaustiveness=8, up to ~60 s)..."):
            docked_path, best_score = run_docking(
                receptor=RECEPTOR_PDBQT, ligand=ligand_pdbqt,
                center=HIV_CENTER, box_size=HIV_BOX, exhaustiveness=8,
            )

        if docked_path is None or best_score is None:
            st.error(
                "AutoDock Vina did not produce a result. "
                "Verify that `vina` is accessible on PATH and that both PDBQT files are valid.\n\n"
                f"Ligand: `{ligand_pdbqt}`\n\nReceptor: `{RECEPTOR_PDBQT}`"
            )
            st.stop()

        # Step 3: Contact analysis
        with st.spinner("Step 3 of 3 - Analysing predicted contacts..."):
            try:
                contacts = analyze_contacts(RECEPTOR_PDB, docked_path)
                seen_disp, unique_contacts = set(), []
                for c in contacts:
                    key = (tuple(c["ligand_coord"]), c["protein_res"], c["protein_atom"])
                    if key not in seen_disp:
                        seen_disp.add(key)
                        unique_contacts.append(c)
            except Exception as e:
                st.warning(f"Contact analysis failed (3D viewer will still display): {e}")
                unique_contacts = []

        # Results
        st.markdown("---")
        col_score, col_detail = st.columns([1, 3])
        with col_score:
            st.metric(
                label="Predicted docking affinity",
                value=f"{best_score:.3f} kcal/mol",
                delta="Lower = more favorable predicted pose",
                delta_color="off",
            )
        with col_detail:
            st.markdown(
                f"**Target:** HIV-1 Protease / 1HVR  \n"
                f"**SMILES:** `{smiles_clean}`  \n"
                f"**Output pose file:** `{os.path.basename(docked_path)}`"
            )

        # 3D Docking Visualization
        st.markdown("### 3D Docking Visualization")
        st.caption(
            "Green sticks: docked ligand pose (AutoDock Vina output). "
            "Spectrum cartoon: HIV-1 Protease. "
            "Gold lines: putative hydrogen bonds. "
            "Light blue lines: hydrophobic contacts. "
            "Drag to rotate; scroll to zoom."
        )
        render_docking_viewer(RECEPTOR_PDB, docked_path, unique_contacts, height=560)

        # Contact tables
        if unique_contacts:
            hbonds = [c for c in unique_contacts if c["type"] == "hbond"]
            hydros = [c for c in unique_contacts if c["type"] == "hydrophobic"]

            col_h, col_c = st.columns(2)
            with col_h:
                st.markdown(f"**Putative hydrogen bonds ({len(hbonds)})**")
                if hbonds:
                    st.dataframe(
                        pd.DataFrame([{
                            "Protein residue": c["protein_res"],
                            "Protein atom":    c["protein_atom"],
                            "Ligand atom":     c["ligand_atom"],
                            "Distance (A)":    c["distance"],
                        } for c in hbonds]).sort_values("Distance (A)"),
                        use_container_width=True, hide_index=True,
                    )
                else:
                    st.info("No putative H-bonds found (criterion: N/O/S pair, distance < 3.5 A).")

            with col_c:
                st.markdown(f"**Hydrophobic contacts ({len(hydros)})**")
                if hydros:
                    st.dataframe(
                        pd.DataFrame([{
                            "Protein residue": c["protein_res"],
                            "Protein atom":    c["protein_atom"],
                            "Ligand atom":     c["ligand_atom"],
                            "Distance (A)":    c["distance"],
                        } for c in hydros[:30]]).sort_values("Distance (A)"),
                        use_container_width=True, hide_index=True,
                    )

            st.markdown(
                '<div class="disclaimer">'
                "Contacts are identified by distance thresholds only: putative H-bonds are "
                "N/O/S atom pairs within 3.5 A; hydrophobic contacts are all atom pairs within "
                "4.0 A. These are computational approximations and are not experimentally "
                "confirmed interactions.</div>",
                unsafe_allow_html=True
            )
        else:
            st.info("No contacts were returned by the contact analyser.")
