# CADD Docker Container — Project Scaffold

## What's in here
- `Dockerfile` + `requirements.txt` — everything needed to build the container
- `scripts/` — the 5 pipeline scripts, each runnable standalone
- `data/` — put your protein PDBs, ligand files, and datasets here
- `results/` — outputs land here (docking poses, QSAR models, ADMET tables)

## Recommended order to actually run things

1. **`admet.py`** — no external data needed, works instantly.
   ```
   python scripts/admet.py "CC(=O)OC1=CC=CC=C1C(=O)O"
   ```

2. **`prepare_ligand.py`** — turn any SMILES into a 3D structure.
   ```
   python scripts/prepare_ligand.py "CC(=O)OC1=CC=CC=C1C(=O)O" aspirin
   ```

3. **`qsar.py`** — the core deliverable. Needs internet access to ChEMBL.
   Pick a target ID from https://www.ebi.ac.uk/chembl/ (search a protein,
   the ID looks like `CHEMBL203`).
   ```
   python scripts/qsar.py --target_chembl_id CHEMBL203
   ```

4. **`prepare_protein.py`** + **`docking.py`** — only once you have a
   specific PDB structure and binding site coordinates. More advanced,
   treat as a stretch goal.

## Building the Docker image
```
docker build -t cadd-container .
docker run -it -v $(pwd)/data:/app/data -v $(pwd)/results:/app/results cadd-container
```
(the `-v` flags mount your local data/results folders into the container,
so results persist after the container exits)

## Status
- `admet.py`, `prepare_ligand.py`, `qsar.py` — tested and working
- `prepare_protein.py`, `docking.py` — written, needs a real PDB + Vina
  installed to test end-to-end (works inside the Docker container)
