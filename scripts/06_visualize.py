#!/usr/bin/env python3
"""Generate publication-ready visualizations of docking results."""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors

BASE_DIR = os.path.join(os.path.dirname(__file__), '..')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
FIGURES_DIR = os.path.join(RESULTS_DIR, 'figures')


def plot_binding_affinity_distribution(df, output_path):
    """Plot distribution of binding affinities across all screened compounds."""
    fig, ax = plt.subplots(figsize=(10, 6))
    affinities = df['affinity_kcal_mol'].values

    ax.hist(affinities, bins=30, color='#2196F3', edgecolor='white', alpha=0.85)
    ax.axvline(x=np.percentile(affinities, 5), color='red', linestyle='--',
               linewidth=2, label=f'Top 5% threshold ({np.percentile(affinities, 5):.1f} kcal/mol)')
    ax.set_xlabel('Binding Affinity (kcal/mol)', fontsize=13)
    ax.set_ylabel('Count', fontsize=13)
    ax.set_title('Virtual Screening: Binding Affinity Distribution\nBfmS Histidine Kinase (PDB: 3KLN)', fontsize=14)
    ax.legend(fontsize=11)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def plot_top_candidates_dashboard(candidates_df, output_path):
    """Create a multi-panel dashboard of top candidates."""
    top = candidates_df.head(10)

    fig = plt.figure(figsize=(16, 14))
    gs = GridSpec(3, 2, figure=fig, hspace=0.35, wspace=0.3)

    # Panel 1: Binding affinity bar chart
    ax1 = fig.add_subplot(gs[0, 0])
    colors = plt.cm.RdYlGn_r(np.linspace(0.1, 0.9, len(top)))
    bars = ax1.barh(range(len(top)), top['affinity_kcal_mol'].values, color=colors, edgecolor='white')
    ax1.set_yticks(range(len(top)))
    ax1.set_yticklabels(top['compound_id'].values, fontsize=9)
    ax1.set_xlabel('Binding Affinity (kcal/mol)', fontsize=11)
    ax1.set_title('Top 10 Candidates: Binding Affinity', fontsize=12, fontweight='bold')
    ax1.invert_yaxis()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Panel 2: QED drug-likeness scores
    ax2 = fig.add_subplot(gs[0, 1])
    qed_colors = plt.cm.viridis(top['QED'].values / top['QED'].max())
    ax2.barh(range(len(top)), top['QED'].values, color=qed_colors, edgecolor='white')
    ax2.set_yticks(range(len(top)))
    ax2.set_yticklabels(top['compound_id'].values, fontsize=9)
    ax2.set_xlabel('QED Score', fontsize=11)
    ax2.set_title('Drug-likeness (QED)', fontsize=12, fontweight='bold')
    ax2.axvline(x=0.5, color='red', linestyle='--', alpha=0.5, label='QED=0.5')
    ax2.legend(fontsize=9)
    ax2.invert_yaxis()
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # Panel 3: MW vs LogP scatter
    ax3 = fig.add_subplot(gs[1, 0])
    scatter = ax3.scatter(
        candidates_df['LogP'], candidates_df['MW'],
        c=candidates_df['affinity_kcal_mol'], cmap='RdYlGn',
        s=60, alpha=0.7, edgecolors='gray', linewidth=0.5
    )
    # Highlight top 3
    top3 = candidates_df.head(3)
    ax3.scatter(top3['LogP'], top3['MW'], s=150, facecolors='none',
                edgecolors='red', linewidth=2.5, label='Top 3')
    ax3.set_xlabel('LogP', fontsize=11)
    ax3.set_ylabel('Molecular Weight (Da)', fontsize=11)
    ax3.set_title('Chemical Space: MW vs LogP', fontsize=12, fontweight='bold')
    # Lipinski boundaries
    ax3.axhline(y=500, color='gray', linestyle=':', alpha=0.5)
    ax3.axvline(x=5, color='gray', linestyle=':', alpha=0.5)
    ax3.legend(fontsize=9)
    plt.colorbar(scatter, ax=ax3, label='Affinity (kcal/mol)')
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    # Panel 4: TPSA vs Affinity
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.scatter(candidates_df['TPSA'], candidates_df['affinity_kcal_mol'],
                c='#2196F3', s=50, alpha=0.6, edgecolors='gray', linewidth=0.5)
    ax4.scatter(top3['TPSA'], top3['affinity_kcal_mol'],
                s=150, facecolors='none', edgecolors='red', linewidth=2.5, label='Top 3')
    ax4.axhline(y=-7.0, color='green', linestyle='--', alpha=0.5, label='Strong binder threshold')
    ax4.set_xlabel('TPSA (A²)', fontsize=11)
    ax4.set_ylabel('Binding Affinity (kcal/mol)', fontsize=11)
    ax4.set_title('TPSA vs Binding Affinity', fontsize=12, fontweight='bold')
    ax4.legend(fontsize=9)
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)

    # Panel 5: Composite score
    ax5 = fig.add_subplot(gs[2, 0])
    colors5 = plt.cm.plasma(top['composite_score'].values / top['composite_score'].max())
    ax5.barh(range(len(top)), top['composite_score'].values, color=colors5, edgecolor='white')
    ax5.set_yticks(range(len(top)))
    ax5.set_yticklabels(top['compound_id'].values, fontsize=9)
    ax5.set_xlabel('Composite Score', fontsize=11)
    ax5.set_title('Overall Ranking Score\n(50% Affinity + 30% QED + 20% Lipinski)', fontsize=12, fontweight='bold')
    ax5.invert_yaxis()
    ax5.spines['top'].set_visible(False)
    ax5.spines['right'].set_visible(False)

    # Panel 6: ADMET summary heatmap
    ax6 = fig.add_subplot(gs[2, 1])
    admet_cols = ['Oral_Bioavailability', 'hERG_Risk', 'CYP_Inhibition_Risk', 'Hepatotoxicity']
    available_cols = [c for c in admet_cols if c in candidates_df.columns]
    if available_cols:
        # Map to numeric for heatmap
        mapping = {'Good': 1, 'High': 0, 'Medium': 0.5, 'Low': 1, 'Poor': 0, 'Yes': 1, 'No': 0, 'Risk': 0}
        heatmap_data = top[available_cols].copy()
        for col in available_cols:
            heatmap_data[col] = heatmap_data[col].map(mapping).fillna(0.5)

        im = ax6.imshow(heatmap_data.values, cmap='RdYlGn', aspect='auto', vmin=0, vmax=1)
        ax6.set_xticks(range(len(available_cols)))
        ax6.set_xticklabels([c.replace('_', '\n') for c in available_cols], fontsize=9)
        ax6.set_yticks(range(len(top)))
        ax6.set_yticklabels(top['compound_id'].values, fontsize=9)
        ax6.set_title('ADMET Profile\n(Green=Favorable, Red=Unfavorable)', fontsize=12, fontweight='bold')

        # Add text annotations
        for i in range(len(top)):
            for j in range(len(available_cols)):
                val = top.iloc[i][available_cols[j]]
                ax6.text(j, i, str(val), ha='center', va='center', fontsize=8)

    fig.suptitle('BfmS Histidine Kinase Drug Discovery Campaign\nVirtual Screening Results Summary',
                 fontsize=16, fontweight='bold', y=1.02)

    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def plot_top3_structures(candidates_df, output_path):
    """Draw 2D structures of top 3 candidates."""
    top3 = candidates_df.head(3)

    mols = []
    legends = []
    for _, row in top3.iterrows():
        mol = Chem.MolFromSmiles(row['smiles'])
        if mol:
            mols.append(mol)
            legends.append(
                f"{row['compound_id']}\n"
                f"Affinity: {row['affinity_kcal_mol']:.1f} kcal/mol\n"
                f"QED: {row['QED']:.3f} | MW: {row['MW']:.0f}"
            )

    if mols:
        img = Draw.MolsToGridImage(
            mols, molsPerRow=3, subImgSize=(400, 400),
            legends=legends
        )
        img.save(output_path)
        print(f"  Saved: {output_path}")


def plot_pocket_context(output_path):
    """Create a schematic of the BfmS pocket with key residues."""
    fig, ax = plt.subplots(figsize=(10, 8))

    # Draw protein outline (schematic)
    theta = np.linspace(0, 2 * np.pi, 100)
    protein_x = 5 + 4 * np.cos(theta)
    protein_y = 5 + 3.5 * np.sin(theta)
    ax.fill(protein_x, protein_y, color='#E3F2FD', edgecolor='#1565C0', linewidth=2)
    ax.text(5, 8, 'BfmS Histidine Kinase\n(PDB: 3KLN)', ha='center', fontsize=14, fontweight='bold')

    # Draw pocket cavity
    pocket_theta = np.linspace(-0.5, 0.5, 50)
    pocket_x = 8.5 + 1.5 * np.cos(pocket_theta)
    pocket_y = 5 + 2.0 * np.sin(pocket_theta)
    ax.fill(pocket_x, pocket_y, color='#FFF3E0', edgecolor='#E65100', linewidth=2)
    ax.text(9, 5, 'Active\nSite\nPocket', ha='center', fontsize=10, fontweight='bold', color='#E65100')

    # Key residues from P2Rank pocket 1
    residues = [
        ('C101', 8.2, 6.2), ('C153', 8.8, 6.5), ('C154', 9.3, 5.8),
        ('C155', 9.5, 5.2), ('C212', 9.2, 4.3), ('C215', 8.7, 3.8),
        ('C216', 8.3, 4.2), ('C66', 7.8, 4.8), ('C69', 7.5, 5.5),
        ('C72', 7.8, 6.0), ('C98', 8.0, 5.5),
    ]
    for name, x, y in residues:
        ax.plot(x, y, 'o', color='#C62828', markersize=8)
        ax.text(x + 0.15, y + 0.15, name, fontsize=7, color='#C62828')

    # Ligand representation
    lig_x = [9.0, 9.3, 9.1, 8.8, 9.0]
    lig_y = [5.0, 5.3, 5.6, 5.3, 5.0]
    ax.fill(lig_x, lig_y, color='#4CAF50', alpha=0.6, edgecolor='#2E7D32', linewidth=2)
    ax.text(9.0, 5.3, 'Ligand', ha='center', fontsize=8, color='white', fontweight='bold')

    # Pocket info
    ax.text(0.5, 1.5,
            'P2Rank Pocket 1\n'
            'Center: (-73.9, 33.0, -21.7) A\n'
            'Score: 4.24\n'
            'Surface atoms: 22\n'
            'Key residues: C66, C69, C72, C98,\n'
            'C101, C153-C155, C212, C215-C216',
            fontsize=10, verticalalignment='bottom',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', edgecolor='gray'))

    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('BfmS Binding Pocket Architecture\nAcinetobacter baumannii (WHO Critical Priority)',
                 fontsize=14, fontweight='bold', pad=20)

    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def generate_all_figures():
    """Generate all publication figures."""
    os.makedirs(FIGURES_DIR, exist_ok=True)

    # Load results
    full_results_path = os.path.join(RESULTS_DIR, 'docking_results.csv')
    candidates_path = os.path.join(RESULTS_DIR, 'final_candidates.csv')

    if not os.path.exists(candidates_path):
        # Use admet_full_results if final_candidates not yet made
        candidates_path = os.path.join(RESULTS_DIR, 'admet_full_results.csv')

    if os.path.exists(full_results_path):
        all_df = pd.read_csv(full_results_path)
        print(f"Loaded {len(all_df)} docking results")
    else:
        print(f"Error: {full_results_path} not found")
        sys.exit(1)

    if os.path.exists(candidates_path):
        cand_df = pd.read_csv(candidates_path)
        print(f"Loaded {len(cand_df)} candidates with ADMET data")
    else:
        print(f"Warning: {candidates_path} not found, using docking results only")
        cand_df = all_df

    print("\nGenerating figures...")

    # Figure 1: Binding affinity distribution
    plot_binding_affinity_distribution(all_df, os.path.join(FIGURES_DIR, 'fig1_affinity_distribution.png'))

    # Figure 2: Top candidates dashboard
    plot_top_candidates_dashboard(cand_df, os.path.join(FIGURES_DIR, 'fig2_candidates_dashboard.png'))

    # Figure 3: Top 3 2D structures
    plot_top3_structures(cand_df, os.path.join(FIGURES_DIR, 'fig3_top3_structures.png'))

    # Figure 4: Pocket context
    plot_pocket_context(os.path.join(FIGURES_DIR, 'fig4_pocket_architecture.png'))

    print("\nAll figures generated successfully!")


if __name__ == '__main__':
    generate_all_figures()
