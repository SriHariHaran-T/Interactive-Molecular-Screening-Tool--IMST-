import os
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

def prepare_1hvr(cif_path, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    parser = MMCIFParser()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', PDBConstructionWarning)
        structure = parser.get_structure('1hvr', cif_path)

    # In 1HVR, the inhibitor is XK2
    ligand_resname = 'XK2'

    # Calculate bounding box from XK2
    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname == ligand_resname:
                    for atom in residue:
                        coords.append(atom.get_coord())

    if coords:
        coords = np.array(coords)
        center = coords.mean(axis=0)
        print(f"Calculated Center of {ligand_resname}: {center}")
    else:
        print(f"Ligand {ligand_resname} not found!")
        return None

    # Save protein only (no water, no ligand)
    class ProteinSelect(Select):
        def accept_residue(self, residue):
            if residue.id[0] != ' ':
                return 0
            return 1

    io = PDBIO()
    io.set_structure(structure)
    pdb_out = os.path.join(output_dir, "1HVR_clean.pdb")
    io.save(pdb_out, ProteinSelect())
    print(f"Saved protein to {pdb_out}")
    return pdb_out, center

if __name__ == "__main__":
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(current_dir, "..", "..", ".."))
    cif_path = os.path.join(project_root, "1hvr.cif")
    out_dir = os.path.join(current_dir, "..", "data", "proteins")
    
    prepare_1hvr(cif_path, out_dir)
