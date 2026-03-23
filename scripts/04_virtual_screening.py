#!/usr/bin/env python3
"""Virtual screening: dock compound library against BfmS top pocket using AutoDock Vina."""

import os
import sys
import subprocess
import tempfile
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

BASE_DIR = os.path.join(os.path.dirname(__file__), '..')
DATA_DIR = os.path.join(BASE_DIR, 'data')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
POSES_DIR = os.path.join(RESULTS_DIR, 'poses')
VINA_BIN = os.path.join(BASE_DIR, 'vina')

# Pocket 1 from P2Rank (top scored pocket on chain C)
POCKET_CENTER = (-73.92, 33.01, -21.74)
BOX_SIZE = (25.0, 25.0, 25.0)


def smiles_to_pdbqt(smiles):
    """Convert SMILES to PDBQT string using Meeko."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result != 0:
        result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if result != 0:
            return None

    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            pass

    try:
        prep = MoleculePreparation()
        mol_setups = prep.prepare(mol)
        if not mol_setups:
            return None
        pdbqt_str, is_ok, err = PDBQTWriterLegacy.write_string(mol_setups[0])
        if is_ok:
            return pdbqt_str
    except Exception:
        pass

    return None


def dock_single(args):
    """Dock a single compound. Returns (name, smiles, affinity, pose_path) or None."""
    smiles, name, receptor_path, vina_bin, poses_dir, center, box = args

    try:
        pdbqt_str = smiles_to_pdbqt(smiles)
        if pdbqt_str is None:
            return None

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdbqt', delete=False) as f:
            f.write(pdbqt_str)
            ligand_path = f.name

        pose_path = os.path.join(poses_dir, f"{name}_pose.pdbqt")

        cmd = [
            vina_bin,
            '--receptor', receptor_path,
            '--ligand', ligand_path,
            '--center_x', str(center[0]),
            '--center_y', str(center[1]),
            '--center_z', str(center[2]),
            '--size_x', str(box[0]),
            '--size_y', str(box[1]),
            '--size_z', str(box[2]),
            '--exhaustiveness', '8',
            '--num_modes', '3',
            '--out', pose_path,
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

        try:
            os.unlink(ligand_path)
        except OSError:
            pass

        if result.returncode != 0:
            return None

        # Parse best affinity from Vina output
        # Look for line like "   1       -5.007          0          0"
        in_results = False
        for line in result.stdout.split('\n'):
            stripped = line.strip()
            if stripped.startswith('-----+'):
                in_results = True
                continue
            if in_results and stripped:
                parts = stripped.split()
                if len(parts) >= 2:
                    try:
                        mode = int(parts[0])
                        affinity = float(parts[1])
                        if mode == 1:
                            return (name, smiles, affinity, pose_path)
                    except (ValueError, IndexError):
                        continue

        return None

    except Exception:
        return None


def run_virtual_screening(compound_file, receptor_path, n_workers=4):
    """Run virtual screening on all compounds."""
    os.makedirs(POSES_DIR, exist_ok=True)

    compounds = []
    with open(compound_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                compounds.append((parts[0], parts[1]))
            elif parts:
                compounds.append((parts[0], f"CPD_{len(compounds):04d}"))

    print(f"Loaded {len(compounds)} compounds for screening")
    print(f"Pocket center: {POCKET_CENTER}")
    print(f"Box size: {BOX_SIZE}")
    print(f"Receptor: {receptor_path}")
    print(f"Workers: {n_workers}")
    print()

    # Test first compound to verify pipeline
    print("Testing first compound...")
    test_args = (compounds[0][0], compounds[0][1], receptor_path, VINA_BIN,
                 POSES_DIR, POCKET_CENTER, BOX_SIZE)
    test_result = dock_single(test_args)
    if test_result:
        print(f"  Test dock OK: {test_result[1]} -> {test_result[2]:.1f} kcal/mol")
    else:
        print("  Test dock FAILED - checking PDBQT generation...")
        pdbqt = smiles_to_pdbqt(compounds[0][0])
        print(f"  PDBQT generated: {pdbqt is not None}")
        if pdbqt:
            print(f"  PDBQT preview:\n{pdbqt[:300]}")

    args_list = [
        (smi, name, receptor_path, VINA_BIN, POSES_DIR, POCKET_CENTER, BOX_SIZE)
        for smi, name in compounds
    ]

    results = []
    completed = 0
    failed = 0

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(dock_single, args): args[1] for args in args_list}
        for future in as_completed(futures):
            name = futures[future]
            completed += 1
            try:
                result = future.result()
                if result is not None:
                    results.append(result)
                    if len(results) <= 5 or len(results) % 50 == 0:
                        print(f"  Hit #{len(results)}: {result[0]} = {result[2]:.1f} kcal/mol")
                else:
                    failed += 1
            except Exception:
                failed += 1

            if completed % 100 == 0:
                print(f"  Progress: {completed}/{len(compounds)} | hits: {len(results)} | failed: {failed}")

    results.sort(key=lambda x: x[2])

    print(f"\nScreening complete!")
    print(f"  Total: {len(compounds)} | Docked: {len(results)} | Failed: {failed}")
    if results:
        print(f"  Best: {results[0][2]:.1f} kcal/mol ({results[0][0]})")
        print(f"  Worst: {results[-1][2]:.1f} kcal/mol")

    # Save all results
    results_csv = os.path.join(RESULTS_DIR, 'docking_results.csv')
    with open(results_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['rank', 'compound_id', 'smiles', 'affinity_kcal_mol', 'pose_file'])
        for i, (name, smiles, affinity, pose_path) in enumerate(results, 1):
            writer.writerow([i, name, smiles, f"{affinity:.1f}", os.path.basename(pose_path)])
    print(f"  All results: {results_csv}")

    # Save top 50
    top50_csv = os.path.join(RESULTS_DIR, 'top50_hits.csv')
    top50 = results[:50]
    with open(top50_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['rank', 'compound_id', 'smiles', 'affinity_kcal_mol', 'pose_file'])
        for i, (name, smiles, affinity, pose_path) in enumerate(top50, 1):
            writer.writerow([i, name, smiles, f"{affinity:.1f}", os.path.basename(pose_path)])
    print(f"  Top 50: {top50_csv}")

    return results


if __name__ == '__main__':
    receptor_path = os.path.join(DATA_DIR, '3KLN_receptor.pdbqt')
    compound_file = os.path.join(DATA_DIR, 'compound_library.smi')
    os.makedirs(RESULTS_DIR, exist_ok=True)
    run_virtual_screening(compound_file, receptor_path, n_workers=6)
