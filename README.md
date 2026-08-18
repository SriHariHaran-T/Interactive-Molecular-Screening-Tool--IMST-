# Interactive Molecular Screening Tool (IMST)

This is a local interactive web application that allows you to pick a drug (or enter your own SMILES string) and see:
1. Whether it's likely to be an effective ORAL drug (absorbed by the human gut) based on standard Lipinski and Veber rules.
2. Its interactive 3D structure.

## Features
- **ADMET Descriptors**: Computes Molecular Weight, LogP, TPSA, H-Bond Donors, H-Bond Acceptors, and Rotatable Bonds using RDKit.
- **Drug-likeness Rules**: Evaluates standard oral drug criteria to check for likelihood of gut absorption and general oral drug profile.
- **3D Molecular Viewer**: Generates a 3D structure from a SMILES string and renders it interactively using `3Dmol.js`.

## Prerequisites
To run this application, you need Python installed along with the following packages:
- `streamlit`
- `pandas`
- `rdkit`

You can install the dependencies via:
```bash
pip install -r requirements.txt
```

## Running the App

Navigate to the `scripts` directory and start the Streamlit app:

```bash
cd scripts
python -m streamlit run app.py
```

Or just:
```bash
cd scripts
streamlit run app.py
```

This will automatically open a browser tab (typically at `http://localhost:8501`) where you can interact with the app.

## Project Structure
- `scripts/app.py`: The main Streamlit web application.
- `data/datasets/drugs.csv`: A sample dataset of drugs to choose from.
- `requirements.txt`: Python package dependencies.

## Access Project Using URL 

https://imstproject.streamlit.app/