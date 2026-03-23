#!/usr/bin/env python3
"""Download ZINC drug-like compound subsets for virtual screening.

Downloads SMILES from ZINC20 drug-like subsets via the ZINC API.
Uses multiple tranches to get a diverse drug-like compound set.
"""

import os
import requests
import time

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')

# ZINC20 drug-like tranches - diverse MW/logP bins
# Format: ZINC-downloader-2D for SMILES
ZINC_TRANCHES = [
    # MW 200-250, logP -1 to 5 (various)
    "https://zinc20.docking.org/tranches/download?tranch_type=2d&mwt=200&logp=0&col=smiles&fmt=txt",
    "https://zinc20.docking.org/tranches/download?tranch_type=2d&mwt=200&logp=1&col=smiles&fmt=txt",
    "https://zinc20.docking.org/tranches/download?tranch_type=2d&mwt=250&logp=0&col=smiles&fmt=txt",
    "https://zinc20.docking.org/tranches/download?tranch_type=2d&mwt=250&logp=1&col=smiles&fmt=txt",
    "https://zinc20.docking.org/tranches/download?tranch_type=2d&mwt=300&logp=1&col=smiles&fmt=txt",
    "https://zinc20.docking.org/tranches/download?tranch_type=2d&mwt=300&logp=2&col=smiles&fmt=txt",
    "https://zinc20.docking.org/tranches/download?tranch_type=2d&mwt=350&logp=1&col=smiles&fmt=txt",
    "https://zinc20.docking.org/tranches/download?tranch_type=2d&mwt=350&logp=2&col=smiles&fmt=txt",
    "https://zinc20.docking.org/tranches/download?tranch_type=2d&mwt=350&logp=3&col=smiles&fmt=txt",
    "https://zinc20.docking.org/tranches/download?tranch_type=2d&mwt=400&logp=2&col=smiles&fmt=txt",
]


def download_zinc_smiles_from_api(target_count=1000):
    """Download drug-like SMILES from ZINC using the substances API."""
    print(f"Downloading drug-like compounds from ZINC20 API...")

    all_smiles = []
    page = 1
    per_page = 100

    while len(all_smiles) < target_count:
        url = (
            f"https://zinc20.docking.org/substances.txt"
            f"?page={page}&count={per_page}"
            f"&mwt_gt=200&mwt_lt=500"
            f"&logp_gt=-1&logp_lt=5"
            f"&purchasability=in-stock"
        )
        try:
            resp = requests.get(url, timeout=30)
            if resp.status_code != 200:
                print(f"  Page {page}: HTTP {resp.status_code}, stopping")
                break
            lines = [l.strip() for l in resp.text.strip().split('\n') if l.strip() and not l.startswith('smiles')]
            if not lines:
                break
            for line in lines:
                parts = line.split()
                if parts:
                    all_smiles.append((parts[0], parts[1] if len(parts) > 1 else f"ZINC_{len(all_smiles)}"))
            print(f"  Page {page}: got {len(lines)} compounds (total: {len(all_smiles)})")
            page += 1
            time.sleep(0.5)
        except Exception as e:
            print(f"  Page {page} error: {e}")
            break

    return all_smiles[:target_count]


def generate_diverse_druglike_library(count=1000):
    """Generate a diverse drug-like compound library using known pharmacophore scaffolds.

    Uses RDKit to enumerate variants of known histidine kinase inhibitor scaffolds
    and drug-like fragments from literature.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    import random

    print("Generating diverse drug-like compound library...")

    # Known HK inhibitor scaffolds and drug-like fragments
    # These are based on published histidine kinase inhibitor motifs
    core_scaffolds = [
        # Imidazole-based (HK mimics)
        "c1cnc[nH]1",
        "c1ccc2[nH]cnc2c1",  # benzimidazole
        "c1ccc2nc[nH]c2c1",
        # Pyridine/pyrimidine cores
        "c1ccncc1",
        "c1ccncn1",
        "c1ncncn1",  # triazine
        # Quinoline/quinazoline
        "c1ccc2ncccc2c1",
        "c1ccc2ncncc2c1",
        # Thiophene/furan
        "c1ccsc1",
        "c1ccoc1",
        # Indole
        "c1ccc2[nH]ccc2c1",
        # Oxazole/thiazole
        "c1cocn1",
        "c1cscn1",
        # Pyrazole
        "c1cn[nH]c1",
        # Biphenyl
        "c1ccc(-c2ccccc2)cc1",
        # Naphthalene
        "c1ccc2ccccc2c1",
        # Phenothiazine
        "c1ccc2c(c1)Sc1ccccc1N2",
        # Benzothiazole
        "c1ccc2scnc2c1",
        # Acridine
        "c1ccc2nc3ccccc3cc2c1",
    ]

    # Functional group decorations
    r_groups = [
        "O", "N", "F", "Cl", "Br",
        "C(=O)O", "C(=O)N", "S(=O)(=O)N",
        "OC", "NC", "NCC", "OCCO",
        "C(F)(F)F", "C#N",
        "C(=O)OCC", "NS(=O)(=O)c1ccccc1",
        "Oc1ccccc1", "Nc1ccccc1",
        "C(=O)Nc1ccccc1",
        "OCc1ccccc1",
        "CC(=O)N",
        "C(=O)NCCO",
        "S(=O)(=O)NC",
    ]

    # Linker groups
    linkers = [
        "", "C", "CC", "CCC", "C(=O)", "C(=O)N",
        "NC(=O)", "OC", "CO", "CS", "SC",
        "NCC", "CCN", "COC", "CCNCC",
    ]

    compounds = set()
    attempts = 0
    max_attempts = count * 20

    random.seed(42)

    while len(compounds) < count and attempts < max_attempts:
        attempts += 1
        scaffold = random.choice(core_scaffolds)
        mol = Chem.MolFromSmiles(scaffold)
        if mol is None:
            continue

        # Randomly add 1-3 substituents
        n_subs = random.randint(1, 3)
        smiles = scaffold
        for _ in range(n_subs):
            rg = random.choice(r_groups)
            linker = random.choice(linkers)

            # Build decorated molecule
            decorated = f"{scaffold}{linker}{rg}"
            mol = Chem.MolFromSmiles(decorated)

            if mol is not None:
                smiles = decorated
                scaffold = smiles

        # Validate drug-likeness (Lipinski)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hba = Descriptors.NumHAcceptors(mol)
        hbd = Descriptors.NumHDonors(mol)

        if 150 < mw < 600 and -2 < logp < 6 and hba <= 10 and hbd <= 5:
            canon = Chem.MolToSmiles(mol)
            if canon not in compounds:
                compounds.add(canon)

    print(f"Generated {len(compounds)} unique drug-like compounds (from {attempts} attempts)")
    return list(compounds)


if __name__ == '__main__':
    os.makedirs(DATA_DIR, exist_ok=True)

    # First try ZINC API
    zinc_compounds = download_zinc_smiles_from_api(target_count=500)

    # Generate additional compounds to reach target
    generated = generate_diverse_druglike_library(count=1000)

    # Combine all compounds
    output_path = os.path.join(DATA_DIR, 'compound_library.smi')
    all_smiles = set()

    with open(output_path, 'w') as f:
        # Add ZINC compounds
        for smi, name in zinc_compounds:
            if smi not in all_smiles:
                all_smiles.add(smi)
                f.write(f"{smi} {name}\n")

        # Add generated compounds
        for i, smi in enumerate(generated):
            if smi not in all_smiles:
                all_smiles.add(smi)
                f.write(f"{smi} GEN_{i:04d}\n")

    print(f"\nTotal unique compounds in library: {len(all_smiles)}")
    print(f"Saved to {output_path}")
