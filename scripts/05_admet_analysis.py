#!/usr/bin/env python3
"""Run ADMET-AI predictions on top docking hits and filter for drug-likeness."""

import os
import sys
import csv
import pandas as pd
import numpy as np

BASE_DIR = os.path.join(os.path.dirname(__file__), '..')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')


def compute_rdkit_descriptors(smiles):
    """Compute key molecular descriptors using RDKit."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    return {
        'MW': round(Descriptors.MolWt(mol), 2),
        'LogP': round(Descriptors.MolLogP(mol), 2),
        'HBA': Descriptors.NumHAcceptors(mol),
        'HBD': Descriptors.NumHDonors(mol),
        'TPSA': round(Descriptors.TPSA(mol), 2),
        'RotBonds': Descriptors.NumRotatableBonds(mol),
        'AromaticRings': Descriptors.NumAromaticRings(mol),
        'HeavyAtoms': Descriptors.HeavyAtomCount(mol),
        'QED': round(Descriptors.qed(mol), 3),
    }


def check_lipinski(desc):
    """Check Lipinski Rule of 5."""
    violations = 0
    if desc['MW'] > 500: violations += 1
    if desc['LogP'] > 5: violations += 1
    if desc['HBA'] > 10: violations += 1
    if desc['HBD'] > 5: violations += 1
    return violations


def run_admet_predictions(smiles_list):
    """Run ADMET-AI predictions on a list of SMILES."""
    try:
        from admet_ai import ADMETModel
        print("Running ADMET-AI predictions...")
        model = ADMETModel()
        results = model.predict(smiles=smiles_list)
        print(f"  ADMET-AI completed for {len(smiles_list)} compounds")
        return results
    except Exception as e:
        print(f"  ADMET-AI error: {e}")
        print("  Falling back to RDKit-based ADMET estimation...")
        return None


def estimate_admet_rdkit(smiles):
    """Estimate ADMET properties using RDKit descriptors as proxy."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}

    tpsa = Descriptors.TPSA(mol)
    logp = Descriptors.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    hbd = Descriptors.NumHDonors(mol)

    # Estimate absorption (Caco-2 permeability proxy)
    # High TPSA = low permeability
    caco2_est = 'High' if tpsa < 90 else ('Medium' if tpsa < 140 else 'Low')

    # Estimate oral bioavailability (Veber rules)
    oral_ba = 'Good' if (tpsa <= 140 and Descriptors.NumRotatableBonds(mol) <= 10) else 'Poor'

    # Estimate BBB penetration
    bbb = 'Yes' if (tpsa < 70 and logp > 0 and mw < 400) else 'No'

    # Estimate hERG liability (logP-based rough estimate)
    herg_risk = 'High' if logp > 3.7 else ('Medium' if logp > 2 else 'Low')

    # CYP inhibition risk
    cyp_risk = 'High' if logp > 4 else ('Medium' if logp > 2.5 else 'Low')

    # Hepatotoxicity estimate
    hepatotox = 'Risk' if mw > 600 or logp > 5 else 'Low'

    return {
        'Caco2_Permeability': caco2_est,
        'Oral_Bioavailability': oral_ba,
        'BBB_Penetration': bbb,
        'hERG_Risk': herg_risk,
        'CYP_Inhibition_Risk': cyp_risk,
        'Hepatotoxicity': hepatotox,
    }


def run_analysis():
    """Run complete ADMET analysis on top 50 docking hits."""
    top50_path = os.path.join(RESULTS_DIR, 'top50_hits.csv')

    if not os.path.exists(top50_path):
        print(f"Error: {top50_path} not found. Run docking first.")
        sys.exit(1)

    df = pd.read_csv(top50_path)
    print(f"Loaded {len(df)} top hits for ADMET analysis")

    # Compute RDKit descriptors
    print("\nComputing molecular descriptors...")
    descriptors = []
    for _, row in df.iterrows():
        desc = compute_rdkit_descriptors(row['smiles'])
        if desc:
            desc['compound_id'] = row['compound_id']
            desc['smiles'] = row['smiles']
            desc['affinity_kcal_mol'] = row['affinity_kcal_mol']
            desc['Lipinski_Violations'] = check_lipinski(desc)
            descriptors.append(desc)

    desc_df = pd.DataFrame(descriptors)
    print(f"  Computed descriptors for {len(desc_df)} compounds")

    # Run ADMET-AI
    smiles_list = desc_df['smiles'].tolist()
    admet_results = run_admet_predictions(smiles_list)

    if admet_results is not None and isinstance(admet_results, pd.DataFrame):
        # Merge ADMET-AI results
        # Select key ADMET columns
        key_cols = [c for c in admet_results.columns if any(
            kw in c.lower() for kw in ['caco', 'bioavail', 'bbb', 'herg', 'cyp', 'hepat', 'toxic', 'solub', 'clearance']
        )]
        if key_cols:
            admet_subset = admet_results[key_cols].copy()
            admet_subset.index = range(len(admet_subset))
            desc_df = pd.concat([desc_df, admet_subset], axis=1)
            print(f"  Merged {len(key_cols)} ADMET-AI features")
    else:
        # Use RDKit-based estimates
        print("  Using RDKit-based ADMET estimates...")
        admet_estimates = []
        for smi in smiles_list:
            admet_estimates.append(estimate_admet_rdkit(smi))
        admet_est_df = pd.DataFrame(admet_estimates)
        desc_df = pd.concat([desc_df.reset_index(drop=True), admet_est_df], axis=1)

    # Score and rank
    # Composite score: weighted combination of affinity, QED, Lipinski compliance
    desc_df['affinity_score'] = (desc_df['affinity_kcal_mol'] - desc_df['affinity_kcal_mol'].max()) / \
                                 (desc_df['affinity_kcal_mol'].min() - desc_df['affinity_kcal_mol'].max())
    desc_df['lipinski_score'] = 1 - (desc_df['Lipinski_Violations'] / 4)
    desc_df['composite_score'] = (
        0.5 * desc_df['affinity_score'] +
        0.3 * desc_df['QED'] +
        0.2 * desc_df['lipinski_score']
    )
    desc_df = desc_df.sort_values('composite_score', ascending=False)
    desc_df['final_rank'] = range(1, len(desc_df) + 1)

    # Filter: keep only drug-like compounds
    druglike = desc_df[desc_df['Lipinski_Violations'] <= 1].copy()
    print(f"\nDrug-like compounds (Lipinski violations <= 1): {len(druglike)}/{len(desc_df)}")

    # Save full results
    full_path = os.path.join(RESULTS_DIR, 'admet_full_results.csv')
    desc_df.to_csv(full_path, index=False)
    print(f"Full results saved to {full_path}")

    # Save filtered drug-like candidates
    candidates_path = os.path.join(RESULTS_DIR, 'final_candidates.csv')
    druglike.to_csv(candidates_path, index=False)
    print(f"Drug-like candidates saved to {candidates_path}")

    # Print top 10 summary
    print("\n" + "=" * 80)
    print("TOP 10 DRUG CANDIDATES FOR BfmS INHIBITION")
    print("=" * 80)
    for _, row in druglike.head(10).iterrows():
        print(f"\n  Rank {int(row['final_rank'])}: {row['compound_id']}")
        print(f"    SMILES: {row['smiles']}")
        print(f"    Binding Affinity: {row['affinity_kcal_mol']:.1f} kcal/mol")
        print(f"    MW: {row['MW']:.1f} | LogP: {row['LogP']:.2f} | QED: {row['QED']:.3f}")
        print(f"    HBA: {row['HBA']} | HBD: {row['HBD']} | TPSA: {row['TPSA']:.1f}")
        print(f"    Lipinski Violations: {row['Lipinski_Violations']}")
        print(f"    Composite Score: {row['composite_score']:.3f}")
        if 'Oral_Bioavailability' in row:
            print(f"    Oral Bioavailability: {row['Oral_Bioavailability']} | hERG Risk: {row['hERG_Risk']}")

    return desc_df


if __name__ == '__main__':
    run_analysis()
