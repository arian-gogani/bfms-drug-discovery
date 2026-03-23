# BfmS Drug Discovery Pipeline - Checkpoint

**Last updated:** 2026-03-23
**Git branch:** main
**Latest commit:** f1ec2e6 (pushed to GitHub)
**Repository:** https://github.com/arian-gogani/bfms-drug-discovery

---

## Project Goal

Design novel small molecule inhibitor candidates for the **BfmS histidine kinase** in *Acinetobacter baumannii* (WHO #1 critical priority pathogen) using AI-driven structure-based virtual screening, then produce results publishable on bioRxiv.

---

## Pipeline Status: ALL STEPS COMPLETE

| Step | Status | Key Result |
|------|--------|------------|
| 1. Install dependencies | DONE | Python 3.14 venv with rdkit, meeko, admet-ai, openmm, biopython, etc. |
| 2. Download PDB 3KLN | DONE | 4 chains (A:217, B:215, C:213, D:173 residues) |
| 3. P2Rank pocket detection | DONE | 5 pockets found; Pocket 1 score 4.24 selected |
| 4. Prepare receptor PDBQT | DONE | 6,640 atoms in receptor file |
| 5. Build compound library (ChEMBL) | DONE | 1,500 Lipinski-compliant compounds from ChEMBL API |
| 6. Virtual screening (Vina) | DONE | 1,459/1,500 docked (97.3% success) |
| 7. ADMET-AI analysis | DONE | 50/50 top hits pass Lipinski, 30 ADMET endpoints |
| 8. Generate figures | DONE | 4 publication-ready PNG figures at 300 DPI |
| 9. Write paper | DONE | 11-page paper (paper.md + paper.pdf), 24 references |
| 10. Push to GitHub | DONE | 3 commits on main |

---

## Key Results

### Screening Summary
- **Library:** 1,500 ChEMBL drug-like compounds (MW 200-500, zero Ro5 violations)
- **Docked:** 1,459 (97.3%), 41 failed (charge state / 3D generation issues)
- **Affinity range:** -9.3 to -4.0 kcal/mol
- **Mean affinity:** -6.48 kcal/mol (SD: 0.74)
- **Top 5% threshold:** -7.8 kcal/mol
- **Compounds < -7.0 kcal/mol:** 305 (20.9%)
- **Compounds < -8.0 kcal/mol:** 26 (1.8%)

### Top 5 Candidates (ranked by composite score)

| Rank | ChEMBL ID | Affinity | MW | QED | Composite |
|------|-----------|----------|-----|------|-----------|
| 1 | CHEMBL7029 | -9.3 | 487.4 | 0.402 | 0.821 |
| 2 | CHEMBL414184 | -9.1 | 473.4 | 0.427 | 0.757 |
| 3 | CHEMBL6748 | -9.1 | 473.4 | 0.424 | 0.756 |
| 4 | CHEMBL7360 | -8.9 | 422.9 | 0.492 | 0.705 |
| 5 | CHEMBL415478 | -8.8 | 475.4 | 0.407 | 0.644 |

### Key Finding
Top 7 candidates converge on an **oxadiazolinedione warhead** motif linked to trifluoromethylphenyl-oxazole ether groups. CHEMBL7062 (rank 7, -8.2 kcal/mol, QED 0.898) is an alternative **diaminopyrimidine** scaffold.

### Binding Pocket (P2Rank Pocket 1)
- **Center:** (-73.92, 33.01, -21.74) Angstrom
- **Chain:** C
- **Score:** 4.24
- **Surface atoms:** 22
- **Key residues:** C66, C69, C72, C98, C101, C153, C154, C155, C212, C215, C216

---

## File Locations

### Data (`data/`)
| File | Description |
|------|-------------|
| `3KLN.pdb` | Raw PDB structure from RCSB |
| `3KLN_clean.pdb` | Cleaned protein (water/heteroatoms removed) |
| `3KLN_receptor.pdbqt` | Receptor in PDBQT format for Vina |
| `chembl_library.smi` | 1,500 ChEMBL SMILES + IDs (active library) |
| `chembl_library_metadata.csv` | ChEMBL compound metadata (ID, SMILES, name, MW) |
| `compound_library.smi` | Old synthetic library (1,000 compounds, superseded) |
| `vina_config.txt` | Vina docking configuration |
| `p2rank_output/` | Full P2Rank output (gitignored) |

### Results (`results/`)
| File | Description |
|------|-------------|
| `docking_results.csv` | All 1,459 docked compounds ranked by affinity |
| `top50_hits.csv` | Top 50 by raw binding affinity |
| `admet_full_results.csv` | Top 50 with RDKit descriptors + 30 ADMET-AI features |
| `final_candidates.csv` | 50 drug-like candidates (Lipinski <= 1 violation), ranked by composite score |
| `poses/` | 1,459 PDBQT binding pose files (gitignored, ~14 MB) |
| `figures/fig1_affinity_distribution.png` | Affinity distribution histogram |
| `figures/fig2_candidates_dashboard.png` | 6-panel candidate dashboard |
| `figures/fig3_top3_structures.png` | 2D structures of top 3 |
| `figures/fig4_pocket_architecture.png` | Pocket schematic with residues |

### Scripts (`scripts/`)
| File | Description |
|------|-------------|
| `01_download_pdb.py` | Download + clean PDB 3KLN |
| `02_prepare_receptor.py` | Convert PDB to PDBQT for Vina |
| `03_download_zinc.py` | Generate synthetic compound library (legacy) |
| `04_virtual_screening.py` | AutoDock Vina parallel docking (accepts library path as CLI arg) |
| `05_admet_analysis.py` | ADMET-AI + RDKit descriptor computation + composite scoring |
| `06_visualize.py` | Generate all 4 publication figures |
| `07_download_chembl.py` | Download drug-like compounds from ChEMBL API |

### Paper
| File | Description |
|------|-------------|
| `paper.md` | Full bioRxiv paper in Markdown (26 KB) |
| `paper.pdf` | Formatted PDF, 11 pages (66 KB) |

---

## Environment Setup (to resume)

```bash
# The venv already exists at .venv/
source .venv/bin/activate

# Key installed packages:
# rdkit, meeko, gemmi, admet-ai (torch, chemprop), openmm,
# biopython, numpy, pandas, matplotlib, requests, scipy, weasyprint

# External tools (not in venv):
# ./vina                     - AutoDock Vina v1.2.5 (macOS aarch64 binary)
# ./p2rank_2.4.2/            - P2Rank v2.4.2 (requires Java 17)
# Java 17 at /opt/homebrew/opt/openjdk@17
# pandoc (for markdown->html)
# poppler (for pdfinfo/pdftoppm)
# pango (for weasyprint PDF generation)
```

---

## Git History

```
f1ec2e6 Rescreen with ChEMBL library: 1459 docked, best hit -9.3 kcal/mol
677fdff Add bioRxiv-ready research paper (paper.md + paper.pdf)
016c3c0 Complete BfmS drug discovery pipeline: 999 compounds docked, top hit -9.0 kcal/mol
754b734 init
```

---

## What Needs to Happen Next

### Immediate computational extensions
1. **Expanded ChEMBL screen** - Screen more of ChEMBL (currently 1,500 of ~2.4M). Increase to 10-50K compounds. Script `07_download_chembl.py` supports this by changing `target_count`.
2. **Molecular dynamics** - Run MD simulations on top 3-5 complexes using OpenMM (already installed) to validate binding stability and identify key interactions.
3. **Focused library around oxadiazolinedione** - Use the SAR insight to enumerate analogs of the top scaffold and screen them.
4. **Multi-pocket screening** - Dock against all 5 P2Rank pockets, not just Pocket 1.
5. **Protein-protein interaction analysis** - Investigate BfmS-BfmR interface as alternative target.

### Paper improvements
6. **Add supplementary tables** - Full ADMET profiles for all 50 candidates.
7. **Molecular interaction diagrams** - Generate 2D interaction maps (hydrogen bonds, hydrophobic contacts) for top 3 poses using PLIP or ProLIF.
8. **Comparison table** - Compare ChEMBL hits vs synthetic library hits side-by-side.

### Experimental validation (wet lab)
9. **Procure top 10 compounds** - ChEMBL compounds are commercially available; check Sigma, MolPort, Enamine.
10. **BfmS autophosphorylation assay** - Test inhibition of BfmS kinase activity.
11. **Biofilm formation assay** - Crystal violet or confocal microscopy on A. baumannii.
12. **MIC/MBC determination** - Minimum inhibitory concentration against CRAB strains.

---

## How to Rerun the Pipeline

```bash
# Full pipeline from scratch:
source .venv/bin/activate
python scripts/01_download_pdb.py
python scripts/02_prepare_receptor.py
python scripts/07_download_chembl.py          # or scripts/03_download_zinc.py for synthetic
python scripts/04_virtual_screening.py data/chembl_library.smi
python scripts/05_admet_analysis.py
python scripts/06_visualize.py

# Regenerate PDF:
pandoc paper.md -o paper.html -s --metadata title="BfmS Drug Discovery"
python -c "from weasyprint import HTML; HTML('paper.html').write_pdf('paper.pdf')"
rm paper.html
```

---

## Composite Score Formula

```
Composite Score = 0.5 * Affinity_norm + 0.3 * QED + 0.2 * Lipinski_score
```
- `Affinity_norm`: min-max normalized (0=weakest, 1=strongest in top 50)
- `QED`: quantitative estimate of drug-likeness (0-1)
- `Lipinski_score`: 1 - (violations / 4)
- Filter: compounds with >1 Lipinski violation excluded
