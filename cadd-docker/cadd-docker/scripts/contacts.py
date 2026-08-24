import json
import math
from Bio.PDB import PDBParser
from rdkit import Chem

def analyze_contacts(protein_pdb, ligand_pdbqt):
    # Parse protein with Biopython
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", protein_pdb)
    
    # Extract protein atoms
    protein_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                # Skip water or heteroatoms
                if residue.id[0] == ' ':
                    for atom in residue:
                        protein_atoms.append({
                            'coord': atom.get_coord(),
                            'element': atom.element,
                            'name': atom.get_name(),
                            'resname': residue.get_resname(),
                            'resid': residue.id[1],
                            'chain': chain.id
                        })
                        
    # Parse ligand with RDKit from PDBQT (treat as PDB for coords)
    # RDKit can read the first model of a PDBQT as PDB if we ignore some lines or just let it try
    # To be safe, we extract the coordinates manually since PDBQT is PDB-like
    ligand_atoms = []
    with open(ligand_pdbqt, 'r') as f:
        for line in f:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                element = line[77:].strip()
                name = line[12:16].strip()
                ligand_atoms.append({
                    'coord': (x, y, z),
                    'element': element,
                    'name': name
                })
            elif line.startswith("ENDMDL"):
                break # Only process the first (best) pose

    contacts = []
    
    for l_atom in ligand_atoms:
        for p_atom in protein_atoms:
            # Calculate distance
            dx = l_atom['coord'][0] - p_atom['coord'][0]
            dy = l_atom['coord'][1] - p_atom['coord'][1]
            dz = l_atom['coord'][2] - p_atom['coord'][2]
            dist = math.sqrt(dx*dx + dy*dy + dz*dz)
            
            if dist < 4.0:
                contact_type = "hydrophobic"
                
                # Simple H-bond heuristic (<3.5A and involving N/O/S)
                hbond_elements = ['N', 'O', 'S']
                l_elem = l_atom['element'].upper()
                p_elem = p_atom['element'].upper()
                
                if dist < 3.5 and l_elem in hbond_elements and p_elem in hbond_elements:
                    contact_type = "hbond"
                    
                contacts.append({
                    'type': contact_type,
                    'distance': round(dist, 2),
                    'ligand_atom': l_atom['name'],
                    'protein_res': f"{p_atom['resname']}{p_atom['resid']}",
                    'protein_atom': p_atom['name'],
                    'protein_coord': [float(p_atom['coord'][0]), float(p_atom['coord'][1]), float(p_atom['coord'][2])],
                    'ligand_coord': [float(l_atom['coord'][0]), float(l_atom['coord'][1]), float(l_atom['coord'][2])]
                })
                
    return contacts

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python contacts.py <protein.pdb> <ligand.pdbqt>")
        sys.exit(1)
        
    protein = sys.argv[1]
    ligand = sys.argv[2]
    result = analyze_contacts(protein, ligand)
    print(json.dumps(result, indent=2))
