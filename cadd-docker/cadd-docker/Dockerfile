# =========================================================
# CADD Docker Container
# Python 3.x + RDKit, Biopython, Open Babel, AutoDock Vina,
# DeepChem, scikit-learn, pandas, NumPy
# =========================================================
FROM python:3.10-slim

# System-level dependencies
# - build-essential: needed to compile some Python packages
# - openbabel: file format conversion (SMILES <-> PDB <-> PDBQT etc.)
# - autodock-vina: the docking engine itself
# - wget: for pulling structures/data inside the container if needed
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    openbabel \
    autodock-vina \
    wget \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Python dependencies (see requirements.txt for versions)
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the project in
COPY scripts/ ./scripts/
COPY data/ ./data/
COPY results/ ./results/

# Default: drop into a shell so you can run scripts interactively
CMD ["/bin/bash"]
