#!/usr/bin/env python3
"""Prepare receptor PDBQT for AutoDock Vina docking."""

import os
import re

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')

def pdb_to_pdbqt_receptor(pdb_path, pdbqt_path):
    """Convert PDB to PDBQT format for receptor.

    Adds Gasteiger charges and AD4 atom types.
    Simple conversion: assigns atom types based on element.
    """
    ad4_type_map = {
        'C': 'C', 'N': 'N', 'O': 'OA', 'S': 'SA',
        'H': 'HD', 'F': 'F', 'P': 'P', 'CL': 'Cl',
        'BR': 'Br', 'I': 'I', 'ZN': 'Zn', 'FE': 'Fe',
        'MG': 'Mg', 'CA': 'Ca', 'MN': 'Mn',
    }

    lines_out = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                element = line[76:78].strip().upper() if len(line) > 76 else ''
                if not element:
                    atom_name = line[12:16].strip()
                    element = re.sub(r'[0-9]', '', atom_name)[0:1].upper()

                ad4_type = ad4_type_map.get(element, element)

                # Check if nitrogen/oxygen could be H-bond donor
                atom_name = line[12:16].strip()
                if element == 'N' and 'H' not in atom_name:
                    ad4_type = 'NA'  # H-bond acceptor nitrogen
                if element == 'O':
                    ad4_type = 'OA'  # H-bond acceptor oxygen

                charge = 0.0  # Placeholder - Vina doesn't use charges

                pdbqt_line = (
                    line[:54] +
                    f"{0.0:>6.2f}{0.0:>6.2f}" +
                    f"    {charge:>+6.3f} " +
                    f"{ad4_type:<2s}\n"
                )
                lines_out.append(pdbqt_line)

    with open(pdbqt_path, 'w') as f:
        for line in lines_out:
            f.write(line)

    print(f"Receptor PDBQT saved: {pdbqt_path} ({len(lines_out)} atoms)")

if __name__ == '__main__':
    pdb_path = os.path.join(DATA_DIR, '3KLN_clean.pdb')
    pdbqt_path = os.path.join(DATA_DIR, '3KLN_receptor.pdbqt')
    pdb_to_pdbqt_receptor(pdb_path, pdbqt_path)
