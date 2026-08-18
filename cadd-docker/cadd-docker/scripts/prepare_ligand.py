"""
prepare_ligand.py
------------------
Takes a molecule as a SMILES string (a text representation of
structure) and turns it into a 3D structure ready for docking.

Steps:
  1. Parse the SMILES string into an RDKit molecule object
  2. Add explicit hydrogens (SMILES usually omits them)
  3. Generate a 3D conformation (SMILES is 2D/flat by default)
  4. Optimize the geometry with a force field (cleans up bond
     lengths/angles so it's physically reasonable)
  5. Save as an SDF file, which Open Babel then converts to .pdbqt

Usage:
    python prepare_ligand.py "CC(=O)OC1=CC=CC=C1C(=O)O" aspirin
    (that SMILES string is aspirin)
"""

import sys
from rdkit import Chem
from rdkit.Chem import AllChem


def prepare_ligand(smiles, name, output_dir="data/ligands"):
    # Step 1: parse SMILES into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")

    # Step 2: add explicit hydrogens
    mol = Chem.AddHs(mol)

    # Step 3: generate a 3D conformation
    # ETKDG is RDKit's standard, reliable 3D embedding algorithm
    AllChem.EmbedMolecule(mol, randomSeed=42)

    # Step 4: optimize geometry with the MMFF force field
    AllChem.MMFFOptimizeMolecule(mol)

    # Step 5: save as SDF (a standard 3D molecule file format)
    output_path = f"{output_dir}/{name}.sdf"
    writer = Chem.SDWriter(output_path)
    writer.write(mol)
    writer.close()

    print(f"3D ligand structure saved to: {output_path}")
    print("Next step (run in terminal, needs Open Babel):")
    print(f"  obabel {output_path} -O {output_path.replace('.sdf', '.pdbqt')}")
    print("  (converts to the .pdbqt format Vina requires)")

    return output_path


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print('Usage: python prepare_ligand.py "<SMILES>" <molecule_name>')
        sys.exit(1)

    smiles_arg = sys.argv[1]
    name_arg = sys.argv[2]
    prepare_ligand(smiles_arg, name_arg)
