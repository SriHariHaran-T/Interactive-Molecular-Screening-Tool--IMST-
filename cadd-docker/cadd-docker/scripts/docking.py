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
import shutil


import re
import traceback

def run_docking(receptor, ligand, center, box_size,
                 output_dir=None, exhaustiveness=8):
    if output_dir is None:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        output_dir = os.path.abspath(os.path.join(current_dir, "..", "results", "docking"))
    os.makedirs(output_dir, exist_ok=True)

    # 1. Verify Vina in PATH
    vina_path = shutil.which("vina")
    print(f"shutil.which('vina'): {vina_path}")
    if vina_path is None:
        path_env = os.environ.get("PATH", "")
        raise RuntimeError(
            f"AutoDock Vina executable ('vina') not found on PATH.\n"
            f"shutil.which('vina') returned None.\n"
            f"PATH environment variable: {path_env}"
        )

    # 2. Verify receptor and ligand PDBQT files
    if not os.path.exists(receptor):
        raise FileNotFoundError(f"Receptor PDBQT file not found: {receptor}")
    if os.path.getsize(receptor) == 0:
        raise ValueError(f"Receptor PDBQT file is empty (0 bytes): {receptor}")

    if not os.path.exists(ligand):
        raise FileNotFoundError(f"Ligand PDBQT file not found: {ligand}")
    if os.path.getsize(ligand) == 0:
        raise ValueError(f"Ligand PDBQT file is empty (0 bytes): {ligand}")

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
    
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True
    )
    
    print(f"Vina exit code: {result.returncode}")
    print(f"Vina STDOUT:\n{result.stdout}")
    print(f"Vina STDERR:\n{result.stderr}")
    
    with open(log_path, "w") as log_file:
        log_file.write(result.stdout or "")
        if result.stderr:
            log_file.write("\nSTDERR:\n" + result.stderr)
            
    if result.returncode != 0:
        raise RuntimeError(
            f"AutoDock Vina failed with return code {result.returncode}\n"
            f"Command: {' '.join(cmd)}\n\n"
            f"STDOUT:\n{result.stdout if result.stdout else '(empty)'}\n\n"
            f"STDERR:\n{result.stderr if result.stderr else '(empty)'}"
        )

    # Parse score from stdout
    # Vina output table looks like:
    # mode |   affinity | dist from best mode
    #      | (kcal/mol) | rmsd l.b.| rmsd u.b.
    # -----+------------+----------+----------
    #    1         -7.4      0.000      0.000
    best_score = None
    for line in result.stdout.split('\n'):
        line = line.strip()
        if line.startswith('1') and len(line.split()) >= 2:
            try:
                best_score = float(line.split()[1])
                break
            except ValueError:
                pass
                
    if best_score is None:
        raise RuntimeError(
            f"AutoDock Vina ran (exit code {result.returncode}) but no affinity score could be parsed from output.\n"
            f"Command: {' '.join(cmd)}\n\n"
            f"STDOUT:\n{result.stdout if result.stdout else '(empty)'}\n\n"
            f"STDERR:\n{result.stderr if result.stderr else '(empty)'}"
        )
                
    print(f"Docking complete. Poses saved to: {out_path}")
    print(f"Best score: {best_score}")
    return out_path, best_score


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
