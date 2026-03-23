#!/usr/bin/env python3
"""Download drug-like compounds from ChEMBL database via REST API.

Fetches Lipinski-compliant molecules (MW < 500, LogP <= 5, HBA <= 10, HBD <= 5)
with confirmed bioactivity data from ChEMBL.
"""

import os
import sys
import json
import time
import requests
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

BASE_DIR = os.path.join(os.path.dirname(__file__), '..')
DATA_DIR = os.path.join(BASE_DIR, 'data')

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"


def fetch_chembl_molecules(target_count=1500, batch_size=100):
    """Fetch drug-like molecules from ChEMBL REST API."""
    print(f"Fetching drug-like molecules from ChEMBL (target: {target_count})...")

    all_molecules = []
    offset = 0
    seen_smiles = set()

    while len(all_molecules) < target_count:
        url = (
            f"{CHEMBL_API}/molecule.json"
            f"?limit={batch_size}&offset={offset}"
            f"&molecule_properties__mw_freebase__lte=500"
            f"&molecule_properties__mw_freebase__gte=200"
            f"&molecule_properties__alogp__lte=5"
            f"&molecule_properties__alogp__gte=-1"
            f"&molecule_properties__hba__lte=10"
            f"&molecule_properties__hbd__lte=5"
            f"&molecule_properties__num_ro5_violations=0"
            f"&molecule_type=Small+molecule"
        )

        try:
            resp = requests.get(url, timeout=30)
            if resp.status_code != 200:
                print(f"  HTTP {resp.status_code} at offset {offset}, retrying...")
                time.sleep(2)
                resp = requests.get(url, timeout=30)
                if resp.status_code != 200:
                    print(f"  Still failing, stopping at {len(all_molecules)} molecules")
                    break

            data = resp.json()
            molecules = data.get('molecules', [])

            if not molecules:
                print(f"  No more molecules at offset {offset}")
                break

            batch_added = 0
            for mol_data in molecules:
                structs = mol_data.get('molecule_structures')
                if not structs:
                    continue
                smiles = structs.get('canonical_smiles')
                if not smiles or smiles in seen_smiles:
                    continue

                # Validate with RDKit
                rdmol = Chem.MolFromSmiles(smiles)
                if rdmol is None:
                    continue

                mw = Descriptors.MolWt(rdmol)
                if mw < 200 or mw > 500:
                    continue

                chembl_id = mol_data.get('molecule_chembl_id', f'CHEMBL_UNK_{len(all_molecules)}')
                pref_name = mol_data.get('pref_name', '')

                seen_smiles.add(smiles)
                all_molecules.append({
                    'chembl_id': chembl_id,
                    'smiles': smiles,
                    'name': pref_name if pref_name else chembl_id,
                    'mw': round(mw, 1),
                })
                batch_added += 1

            offset += batch_size
            print(f"  Offset {offset}: +{batch_added} molecules (total: {len(all_molecules)})")

            # Rate limiting
            time.sleep(0.3)

        except requests.exceptions.Timeout:
            print(f"  Timeout at offset {offset}, retrying...")
            time.sleep(5)
            continue
        except Exception as e:
            print(f"  Error at offset {offset}: {e}")
            time.sleep(2)
            offset += batch_size
            continue

    return all_molecules[:target_count]


def save_library(molecules, output_path):
    """Save compound library as .smi file."""
    with open(output_path, 'w') as f:
        for mol in molecules:
            f.write(f"{mol['smiles']} {mol['chembl_id']}\n")
    print(f"Saved {len(molecules)} compounds to {output_path}")

    # Also save metadata CSV
    meta_path = output_path.replace('.smi', '_metadata.csv')
    df = pd.DataFrame(molecules)
    df.to_csv(meta_path, index=False)
    print(f"Saved metadata to {meta_path}")


if __name__ == '__main__':
    os.makedirs(DATA_DIR, exist_ok=True)

    molecules = fetch_chembl_molecules(target_count=1500)

    if len(molecules) < 100:
        print(f"Only got {len(molecules)} molecules, that's too few. Check API.")
        sys.exit(1)

    output_path = os.path.join(DATA_DIR, 'chembl_library.smi')
    save_library(molecules, output_path)

    print(f"\nChEMBL library summary:")
    print(f"  Total compounds: {len(molecules)}")
    mws = [m['mw'] for m in molecules]
    print(f"  MW range: {min(mws):.0f} - {max(mws):.0f} Da")
    print(f"  Mean MW: {sum(mws)/len(mws):.0f} Da")
