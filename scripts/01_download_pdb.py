#!/usr/bin/env python3
"""Download and prepare PDB structure 3KLN (BfmS histidine kinase)."""

import os
import requests
from Bio.PDB import PDBParser, PDBIO, Select

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')
PDB_ID = '3KLN'

class ProteinSelect(Select):
    """Select only protein atoms (remove water, ligands, etc.)."""
    def accept_residue(self, residue):
        return residue.id[0] == ' '

def download_pdb(pdb_id, output_dir):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    output_path = os.path.join(output_dir, f"{pdb_id}.pdb")

    print(f"Downloading {pdb_id} from RCSB...")
    response = requests.get(url)
    response.raise_for_status()

    with open(output_path, 'w') as f:
        f.write(response.text)

    print(f"Saved raw PDB to {output_path}")
    return output_path

def clean_pdb(input_path, output_dir):
    """Remove waters and heteroatoms, keep only protein chain A."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_path)

    # Print structure info
    for model in structure:
        for chain in model:
            residues = [r for r in chain if r.id[0] == ' ']
            print(f"  Chain {chain.id}: {len(residues)} residues")

    # Save clean protein
    io = PDBIO()
    io.set_structure(structure)
    clean_path = os.path.join(output_dir, "3KLN_clean.pdb")
    io.save(clean_path, ProteinSelect())
    print(f"Saved clean protein to {clean_path}")
    return clean_path

if __name__ == '__main__':
    os.makedirs(DATA_DIR, exist_ok=True)
    raw_path = download_pdb(PDB_ID, DATA_DIR)
    clean_path = clean_pdb(raw_path, DATA_DIR)
    print("Done.")
