# Interactive Molecular Screening Tool (IMST)

Interactive Molecular Screening Tool (IMST) is a Streamlit-based computational platform for preliminary molecular screening, ADMET analysis, molecular docking, protein-ligand contact analysis, and interactive 3D visualization.

The application combines cheminformatics and molecular docking tools into a single web interface, allowing users to examine small molecules and perform docking against a selected protein target.

> **Important:** IMST provides computational predictions and rule-based heuristics. Its results do not establish clinical safety, efficacy, toxicity, therapeutic suitability, or regulatory approval.

## Features

### ADMET Profiling

Computes commonly used molecular descriptors using RDKit:

- Molecular Weight
- LogP
- Topological Polar Surface Area (TPSA)
- Hydrogen Bond Donors
- Hydrogen Bond Acceptors
- Rotatable Bonds

The application also evaluates common Lipinski and Veber structural heuristics.

These rules are presented as preliminary drug-likeness and absorption heuristics rather than definitive assessments of drug safety or efficacy.

### Molecular Structure Generation

Users can enter a SMILES string to generate and visualize a 3D molecular structure.

The application uses RDKit for molecular structure processing and 3D conformer generation.

### Molecular Docking

IMST supports molecular docking against the HIV-1 Protease target represented by PDB structure **1HVR**.

The docking pipeline includes:

1. SMILES input
2. Ligand preparation using RDKit
3. Molecular format conversion using Open Babel
4. PDBQT preparation
5. Molecular docking using AutoDock Vina
6. Docking affinity extraction
7. Docked-pose visualization

The current live docking configuration uses:

- **Target:** HIV-1 Protease (1HVR)
- **Box center:** `(-9.19, 15.91, 27.95)`
- **Box size:** `20 × 20 × 20 Å`
- **Vina exhaustiveness:** `8`

### Protein Preparation

The provided `1hvr.cif` structure is processed to generate a clean protein structure for docking.

The workflow prepares:

```text
1hvr.cif
    ↓
1HVR_clean.pdb
    ↓
1HVR_clean.pdbqt
````

### Ligand Preparation

Ligands can be supplied as SMILES strings.

The preparation workflow is:

```text
SMILES
   ↓
RDKit
   ↓
3D molecular structure
   ↓
SDF
   ↓
Open Babel
   ↓
PDBQT
```

### Protein-Ligand Contact Analysis

After docking, the resulting pose can be analyzed for geometric interactions between the ligand and protein.

The contact-analysis component uses RDKit and Biopython to identify:

* Close protein-ligand contacts
* Putative hydrogen-bond interactions
* Geometric interaction information

These interactions are computationally inferred from molecular geometry and should not be interpreted as experimentally confirmed interactions.

### Interactive 3D Visualization

The docking result can be visualized interactively using `3Dmol.js`.

The visualization can display:

* HIV-1 Protease
* The actual AutoDock Vina docked ligand pose
* Putative hydrogen bonds
* Hydrophobic contacts
* Protein-ligand geometry

## Technology Stack

### Frontend / Application

* Streamlit
* Python
* 3Dmol.js

### Cheminformatics

* RDKit
* Open Babel

### Molecular Docking

* AutoDock Vina

### Structural Analysis

* Biopython
* RDKit

### Data Processing

* Pandas

## Prerequisites

The application requires Python and the following Python packages:

* `streamlit`
* `pandas`
* `rdkit`
* `biopython`
* `scikit-learn`
* `numpy`
* `chembl_webresource_client`
* `deepchem`

The application also requires the following external command-line tools:

* Open Babel
* AutoDock Vina

Both executables must be available through the system `PATH`. For Streamlit Cloud deployment, these are specified in `packages.txt`.

Verify the installations using:

```bash
obabel -V
vina --version
```

## Installation

Install the Python dependencies using:

```bash
pip install -r requirements.txt
```

Make sure Open Babel and AutoDock Vina are installed and accessible from the command line.

## Running the Application

Navigate to the `scripts` directory:

```bash
cd scripts
```

Then start the Streamlit application:

```bash
python -m streamlit run app.py
```

Alternatively:

```bash
streamlit run app.py
```

The application will normally be available at:

```text
http://localhost:8501
```

## Application Workflow

The overall computational workflow is:

```text
                    User Input
                        │
                        ▼
                  SMILES String
                        │
                        ▼
                      RDKit
                        │
                        ▼
                Ligand Preparation
                        │
                        ▼
                   Open Babel
                        │
                        ▼
                     PDBQT
                        │
                        ▼
                 AutoDock Vina
                        │
                  ┌─────┴─────┐
                  ▼           ▼
             Docking      Docked Pose
              Score           │
                  │           ▼
                  │     Contact Analysis
                  │           │
                  └─────┬─────┘
                        ▼
                 3D Visualization
                        │
                        ▼
                  Streamlit UI
```

## Project Structure

```text
IMST/
├── scripts/
│   ├── app.py
│   ├── docking.py
│   ├── prepare_ligand.py
│   ├── prepare_1hvr.py
│   └── contacts.py
│
├── data/
│   ├── datasets/
│   │   └── drugs.csv
│   └── proteins/
│       ├── 1HVR_clean.pdb
│       └── 1HVR_clean.pdbqt
│
├── 1hvr.cif
├── requirements.txt
└── README.md
```

## Current Target

The current docking demonstration uses:

**HIV-1 Protease (1HVR)**

The docking search region is centered around the binding site identified from the co-crystallized ligand.

This configuration is intended for computational demonstration and preliminary screening.

## Limitations

IMST is a computational screening and visualization platform. Its results have important limitations:

* Lipinski and Veber rules are structural heuristics.
* Docking scores are computational estimates rather than experimental binding affinities.
* Putative molecular contacts are inferred from geometry.
* Docking does not demonstrate biological activity.
* A favorable docking score does not establish that a compound is an effective or safe drug.
* Experimental validation is required before making biological or clinical conclusions.

## Web Application

The deployed Streamlit application is available at:

[https://imstproject.streamlit.app/](https://imstproject.streamlit.app/)

