"""
prepare_protein.py
-------------------
Takes a raw PDB file (as downloaded from rcsb.org) and cleans it up
so it's ready for docking:
  1. Removes water molecules (HOH) — they're not part of the actual
     protein and just get in the way during docking.
  2. Removes any co-crystallized ligand (we dock our OWN ligand later,
     not the one that happened to be in the crystal structure).
  3. Saves a clean protein-only PDB file.

Converting that clean PDB -> .pdbqt (the format Vina needs, which
includes partial charges and torsion info) is best done with
Open Babel from the command line, shown at the bottom of this file.

Usage:
    python prepare_protein.py data/proteins/1HVR.pdb
"""

import sys
from Bio.PDB import PDBParser, PDBIO, Select


class ProteinOnlySelect(Select):
    """Keeps only standard amino acid residues — drops water (HOH)
    and any hetero-atoms (co-crystallized ligands, ions, etc.)."""

    def accept_residue(self, residue):
        # residue.id[0] == " " means it's a standard amino acid residue
        # anything else (e.g. "W" for water, "H_LIG" for a hetero group)
        # gets excluded
        return residue.id[0] == " "


def clean_protein(input_pdb_path, output_pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb_path)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb_path, ProteinOnlySelect())

    print(f"Cleaned protein saved to: {output_pdb_path}")
    print("Next step (run in terminal, needs Open Babel):")
    print(f"  obabel {output_pdb_path} -O {output_pdb_path.replace('.pdb', '.pdbqt')} -xr")
    print("  (-xr adds hydrogens and prepares it as a rigid receptor for Vina)")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python prepare_protein.py <path_to_raw_pdb_file>")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = input_path.replace(".pdb", "_clean.pdb")
    clean_protein(input_path, output_path)
