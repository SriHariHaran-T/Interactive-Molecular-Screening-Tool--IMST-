"""
docking.py
----------
Runs AutoDock Vina to dock a prepared ligand (.pdbqt) against a
prepared protein/receptor (.pdbqt), inside a defined 3D search box.

You need, before running this:
  - A receptor .pdbqt file      (from prepare_protein.py + obabel)
  - A ligand .pdbqt file        (from prepare_ligand.py + obabel)
  - The center + size of a search box around the binding site
    (you find this by looking at where a co-crystallized ligand
    sits in the original PDB file, or with a pocket-finder tool
    like fpocket)

This script just wraps the Vina command line — Vina itself does
the actual pose search + scoring, this is glue code.

Usage:
    python docking.py --receptor data/proteins/1HVR_clean.pdbqt \
                       --ligand data/ligands/aspirin.pdbqt \
                       --center_x 10.0 --center_y 5.0 --center_z 15.0 \
                       --size_x 20 --size_y 20 --size_z 20
"""

import argparse
import subprocess
import os


def run_docking(receptor, ligand, center, box_size,
                 output_dir="results/docking", exhaustiveness=8):
    os.makedirs(output_dir, exist_ok=True)

    ligand_name = os.path.splitext(os.path.basename(ligand))[0]
    out_path = os.path.join(output_dir, f"{ligand_name}_docked.pdbqt")
    log_path = os.path.join(output_dir, f"{ligand_name}_log.txt")

    cmd = [
        "vina",
        "--receptor", receptor,
        "--ligand", ligand,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(box_size[0]),
        "--size_y", str(box_size[1]),
        "--size_z", str(box_size[2]),
        "--exhaustiveness", str(exhaustiveness),
        "--out", out_path,
    ]

    print("Running:", " ".join(cmd))
    with open(log_path, "w") as log_file:
        result = subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT)

    if result.returncode != 0:
        print(f"Vina failed — check the log at {log_path}")
        return None

    print(f"Docking complete. Poses saved to: {out_path}")
    print(f"Full log (includes scores) saved to: {log_path}")
    print("Lower (more negative) score = predicted stronger binding.")
    return out_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--receptor", required=True)
    parser.add_argument("--ligand", required=True)
    parser.add_argument("--center_x", type=float, required=True)
    parser.add_argument("--center_y", type=float, required=True)
    parser.add_argument("--center_z", type=float, required=True)
    parser.add_argument("--size_x", type=float, default=20)
    parser.add_argument("--size_y", type=float, default=20)
    parser.add_argument("--size_z", type=float, default=20)
    args = parser.parse_args()

    run_docking(
        receptor=args.receptor,
        ligand=args.ligand,
        center=(args.center_x, args.center_y, args.center_z),
        box_size=(args.size_x, args.size_y, args.size_z),
    )
