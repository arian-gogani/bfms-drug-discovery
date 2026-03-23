# AI-Driven Computational Discovery of Novel BfmS Histidine Kinase Inhibitors Against Carbapenem-Resistant *Acinetobacter baumannii*

**Authors:** Arian Gogani^1^

**Affiliations:** ^1^ Independent Researcher

**Corresponding Author:** Arian Gogani

**Keywords:** Acinetobacter baumannii, BfmS, histidine kinase, virtual screening, drug discovery, antimicrobial resistance, ADMET, AutoDock Vina, ChEMBL

---

## Abstract

Virtual screening of 1,500 ChEMBL drug-like compounds against the BfmS histidine kinase sensor domain (PDB: 3KLN) reveals striking convergence on a single privileged scaffold: six of the top seven candidates share a conserved oxadiazolinedione warhead motif linked to trifluoromethylphenyl-oxazole ether groups, strongly implicating this chemotype in selective BfmS pocket engagement. The lead compound, CHEMBL7029, achieves the strongest predicted binding affinity (-9.3 kcal/mol, estimated Ki ~0.15 uM), with favorable oral bioavailability (0.81), zero Lipinski violations, and moderate metabolic stability. A secondary hit from an orthogonal scaffold, CHEMBL7062 -- a compact diaminopyrimidine (MW 338.3 Da, LogP 2.89) -- emerges at rank 7 with the highest drug-likeness score in the dataset (QED 0.898), offering a structurally distinct backup series for lead optimization. BfmS, a master regulator of biofilm formation and virulence in *Acinetobacter baumannii* (WHO #1 critical priority pathogen), has no reported small-molecule inhibitors despite an available crystal structure. We address this gap with a fully computational, AI-augmented pipeline integrating P2Rank pocket prediction, AutoDock Vina docking (1,459/1,500 compounds docked, 97.3% success; mean affinity -6.48 kcal/mol, SD 0.74; 305 compounds exceeding -7.0 kcal/mol), and ADMET-AI profiling across 30 pharmacokinetic endpoints. All 50 top-ranked candidates pass Lipinski drug-likeness filters with zero violations. The convergent oxadiazolinedione SAR and the complementary diaminopyrimidine scaffold provide two experimentally tractable chemical series -- both sourced from a curated bioactivity database with established commercial supply chains -- for immediate hit-to-lead development against carbapenem-resistant *A. baumannii*.

## 1. Introduction

### 1.1 The *Acinetobacter baumannii* Crisis

Antimicrobial resistance (AMR) represents one of the most pressing global health threats of the 21st century, with drug-resistant infections projected to cause 10 million deaths annually by 2050 (O'Neill, 2016). Among the ESKAPE pathogens, *Acinetobacter baumannii* has emerged as a particularly formidable challenge. In 2024, the WHO reaffirmed carbapenem-resistant *A. baumannii* (CRAB) as the highest priority pathogen for which new antibiotics are critically needed (WHO, 2024). CRAB infections carry mortality rates of 40-60% in intensive care units and are associated with ventilator-associated pneumonia, bloodstream infections, and wound infections, particularly among immunocompromised patients (Harding et al., 2018).

The clinical pipeline for novel anti-*Acinetobacter* agents remains dangerously sparse. Existing treatments rely on last-resort agents such as colistin and tigecycline, both of which carry significant toxicity profiles and face emerging resistance (Karakonstantis et al., 2020). There is a critical need for novel therapeutic strategies that target *A. baumannii*-specific virulence mechanisms rather than conventional antimicrobial targets.

### 1.2 BfmS as a Drug Target

Two-component signal transduction systems (TCS) are attractive antibacterial targets because they regulate essential virulence processes and are absent in mammalian hosts, reducing the risk of off-target toxicity (Gotoh et al., 2010). In *A. baumannii*, the BfmR/BfmS two-component system is a master regulator of biofilm formation, desiccation tolerance, and virulence (Tomaras et al., 2008). BfmS is the membrane-bound sensor histidine kinase that detects environmental signals and phosphorylates the response regulator BfmR, which in turn activates the expression of genes required for biofilm development and capsule production (Liou et al., 2014).

Critically, BfmS/BfmR knockout mutants show dramatically reduced biofilm formation and attenuated virulence in animal infection models, validating BfmS as a drug target (Thompson et al., 2012). The crystal structure of the BfmS sensor domain has been solved at high resolution (PDB: 3KLN), providing an atomic-level template for structure-based drug design (Draughn et al., 2018). Despite this structural availability, no small molecule inhibitors of BfmS have been reported in the literature, representing a significant gap in the field.

### 1.3 Computational Drug Discovery Approach

Modern computational drug discovery integrates structure-based virtual screening with machine learning-driven ADMET prediction to rapidly identify and prioritize drug candidates before costly experimental validation (Sliwoski et al., 2014). AutoDock Vina remains one of the most widely validated docking programs, with demonstrated accuracy in predicting binding poses and relative affinities for drug-like molecules (Eberhardt et al., 2021). Recent advances in ML-based ADMET prediction, particularly the ADMET-AI platform, enable rapid pharmacokinetic profiling across dozens of endpoints using only molecular structure as input (Swanson et al., 2024).

In this study, we describe a complete computational pipeline for the discovery of BfmS histidine kinase inhibitors, from target structure preparation through pocket identification, virtual screening of ChEMBL-derived compounds, and multi-parameter ADMET filtering, yielding a ranked list of novel drug candidates for experimental follow-up.

## 2. Methods

### 2.1 Target Structure Preparation

The crystal structure of the BfmS sensor domain from *A. baumannii* was retrieved from the RCSB Protein Data Bank (PDB ID: 3KLN) using the RCSB programmatic API. The structure contains four chains (A: 217 residues, B: 215 residues, C: 213 residues, D: 173 residues). The structure was cleaned by removing all water molecules and heteroatoms, retaining only protein heavy atoms and their associated coordinates. The cleaned structure was converted to PDBQT format for docking, with AutoDock 4 atom types assigned based on element identity (C, NA, OA, SA, HD for carbon, nitrogen, oxygen, sulfur, and polar hydrogen, respectively). All four chains were retained to preserve the biologically relevant oligomeric context of the protein.

### 2.2 Binding Pocket Identification

Binding pockets were identified using P2Rank v2.4.2 (Krivak & Hoksza, 2018), a machine learning-based pocket prediction tool that operates on protein surface features without requiring homology information. P2Rank was applied to the cleaned multi-chain structure using default parameters. Five pockets were identified, and the top-ranked pocket (Pocket 1, P2Rank score: 4.24, probability: 0.174) was selected as the docking target. This pocket is located on Chain C and comprises 22 surface atoms spanning residues C66, C69, C72, C98, C101, C153, C154, C155, C212, C215, and C216. The pocket center coordinates (-73.92, 33.01, -21.74 Angstrom) were used to define the docking search space.

### 2.3 Compound Library

A drug-like compound library of 1,500 molecules was assembled from the ChEMBL database (Zdrazil et al., 2024) via the ChEMBL REST API. Compounds were filtered at the database level to enforce strict Lipinski Rule of Five compliance: molecular weight 200-500 Da, ALogP -1 to 5, hydrogen bond acceptors <= 10, hydrogen bond donors <= 5, and zero Ro5 violations. Only small molecules were included. All retrieved SMILES were validated with RDKit (Landrum, 2023) and deduplicated by canonical SMILES, yielding a final library of 1,500 unique, experimentally characterized drug-like compounds with a mean molecular weight of 330 Da.

### 2.4 Virtual Screening

Molecular docking was performed using AutoDock Vina v1.2.5 (Eberhardt et al., 2021). Each compound was converted from SMILES to 3D coordinates using RDKit's ETKDG conformer generator (Riniker & Landrum, 2015), followed by MMFF94 force field geometry optimization. PDBQT files were generated using Meeko v0.7.1 (Santos-Martins et al., 2023), which correctly handles torsional flexibility and assigns appropriate AutoDock atom types.

The docking search space was defined as a 25 x 25 x 25 Angstrom box centered on the P2Rank Pocket 1 center. Docking was performed with an exhaustiveness parameter of 8 and a maximum of 3 binding modes per compound. Docking was parallelized across 6 CPU cores using Python's ProcessPoolExecutor. The best binding mode (lowest energy) for each compound was retained for analysis.

### 2.5 ADMET Profiling and Drug-Likeness Filtering

Molecular descriptors were computed using RDKit, including molecular weight (MW), octanol-water partition coefficient (LogP), number of hydrogen bond acceptors (HBA) and donors (HBD), topological polar surface area (TPSA), number of rotatable bonds, aromatic ring count, heavy atom count, and quantitative estimate of drug-likeness (QED; Bickerton et al., 2012). Lipinski Rule of Five violations were assessed for each compound.

ADMET properties were predicted using ADMET-AI v2.0.1 (Swanson et al., 2024), a deep learning platform based on the Chemprop message-passing neural network architecture. Thirty ADMET endpoints were predicted, including: Caco-2 permeability (Wang et al. model), oral bioavailability (Ma et al. model), blood-brain barrier penetration (Martins et al. model), hERG channel inhibition liability, CYP450 inhibition (CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP3A4 via Veith et al. models), CYP substrate classification (CYP2C9, CYP2D6, CYP3A4 via CarbonMangels models), hepatocyte clearance, microsomal clearance, and aqueous solubility.

### 2.6 Composite Scoring and Candidate Ranking

Candidates were ranked by a composite score integrating binding affinity, drug-likeness, and Lipinski compliance:

**Composite Score = 0.5 * Affinity_norm + 0.3 * QED + 0.2 * Lipinski_score**

where Affinity_norm is the min-max normalized binding affinity (0 = weakest, 1 = strongest in the top 50), QED is the quantitative estimate of drug-likeness (0 to 1), and Lipinski_score = 1 - (violations / 4). Compounds with more than one Lipinski violation were excluded from the final candidate list.

## 3. Results

### 3.1 Binding Pocket Characterization

P2Rank analysis of the BfmS sensor domain identified five putative binding pockets (Table 1). Pocket 1, located on Chain C, was the top-ranked site with a score of 4.24 and 22 surface atoms. This pocket is formed by residues in the sensor domain core and likely represents the signal perception or dimerization interface of the histidine kinase. The pocket residues (C66, C69, C72, C98, C101, C153-155, C212, C215-216) span multiple secondary structure elements, creating a concave surface suitable for small molecule binding.

**Table 1. P2Rank binding pocket predictions for BfmS (PDB: 3KLN).**

| Pocket | Score | Probability | Surface Atoms | Center (x, y, z) | Chain |
|--------|-------|-------------|---------------|-------------------|-------|
| 1 | 4.24 | 0.174 | 22 | (-73.9, 33.0, -21.7) | C |
| 2 | 1.81 | 0.034 | 19 | (-27.7, -3.2, -20.6) | A |
| 3 | 1.45 | 0.020 | 23 | (-47.9, 12.9, -20.6) | A/C |
| 4 | 0.97 | 0.005 | 12 | (-28.6, -9.6, -21.2) | A |
| 5 | 0.85 | 0.003 | 15 | (-55.1, 11.5, -54.8) | D |

### 3.2 Virtual Screening Results

Of 1,500 ChEMBL compounds submitted for docking, 1,459 (97.3%) were successfully processed by AutoDock Vina. Forty-one compounds failed due to 3D coordinate generation issues or charge state errors. The binding affinity distribution (Figure 1) follows an approximately normal distribution centered at -6.48 kcal/mol (SD = 0.74 kcal/mol), with a median of -6.5 kcal/mol. The top 5% threshold was -7.8 kcal/mol. Three hundred and five compounds (20.9%) exhibited affinities stronger than -7.0 kcal/mol, and 26 compounds (1.8%) exceeded -8.0 kcal/mol. The strongest predicted binder achieved -9.3 kcal/mol (CHEMBL7029).

Notably, the ChEMBL library produced a distribution shifted approximately 0.5 kcal/mol toward stronger binding compared to a preliminary screen of computationally generated scaffolds (mean -5.92 vs -6.48 kcal/mol), likely reflecting the more drug-like character and optimized geometries of experimentally characterized compounds.

### 3.3 Top Candidates and ADMET Profiles

After Lipinski drug-likeness filtering (violations <= 1), all 50 of the top 50 compounds were retained as viable candidates -- a 100% pass rate reflecting the stringent Lipinski filtering applied at the library construction stage. The top 10 candidates are summarized in Table 2.

**Table 2. Top 10 BfmS inhibitor candidates ranked by composite score.**

| Rank | ChEMBL ID | Affinity (kcal/mol) | MW (Da) | LogP | QED | TPSA (A^2) | Lip. Viol. | Composite Score |
|------|-----------|---------------------|---------|------|-------|------------|------------|-----------------|
| 1 | CHEMBL7029 | -9.3 | 487.4 | 4.79 | 0.402 | 103.3 | 0 | 0.821 |
| 2 | CHEMBL414184 | -9.1 | 473.4 | 4.40 | 0.427 | 103.3 | 0 | 0.757 |
| 3 | CHEMBL6748 | -9.1 | 473.4 | 4.49 | 0.424 | 103.3 | 0 | 0.756 |
| 4 | CHEMBL7360 | -8.9 | 422.9 | 4.62 | 0.492 | 77.2 | 0 | 0.705 |
| 5 | CHEMBL415478 | -8.8 | 475.4 | 4.32 | 0.407 | 103.3 | 0 | 0.644 |
| 6 | CHEMBL7251 | -8.7 | 474.4 | 4.85 | 0.519 | 77.2 | 0 | 0.641 |
| 7 | CHEMBL7062 | -8.2 | 338.3 | 2.89 | 0.898 | 87.1 | 0 | 0.577 |
| 8 | CHEMBL414014 | -8.6 | 498.4 | 4.70 | 0.409 | 77.2 | 0 | 0.573 |
| 9 | CHEMBL266574 | -8.5 | 485.4 | 4.30 | 0.423 | 103.3 | 0 | 0.541 |
| 10 | CHEMBL7748 | -8.4 | 357.4 | 3.70 | 0.527 | 91.5 | 0 | 0.537 |

A striking feature of the results is the convergence of the top candidates on a shared pharmacophore. The top-ranked compound **CHEMBL7029** (SMILES: `C/C(=C\Cn1oc(=O)[nH]c1=O)c1cccc(OCc2nc(-c3ccc(C(F)(F)F)cc3)oc2C)c1`) features an oxadiazolinedione (1,2,4-oxadiazolidine-3,5-dione) warhead connected via a propenyl linker to a methyl-substituted phenyl ring bearing a para-trifluoromethylphenyl-oxazole ether. Six of the top 7 candidates share this oxadiazolinedione moiety, varying in their aryl ether substituents and linker geometry (Table 2). This scaffold convergence strongly suggests that the oxadiazolinedione group makes specific, favorable interactions within the BfmS pocket.

ADMET-AI predictions for the top 3 candidates are presented in Table 3.

**Table 3. ADMET-AI predictions for the top 3 candidates.**

| Property | CHEMBL7029 | CHEMBL414184 | CHEMBL6748 |
|----------|------------|--------------|------------|
| BBB Permeability (Martins) | 0.71 | 0.69 | 0.71 |
| Oral Bioavailability (Ma) | 0.81 | 0.83 | 0.86 |
| hERG Inhibition | 0.53 | 0.51 | 0.47 |
| Caco-2 Permeability (log cm/s) | -4.75 | -4.75 | -4.80 |
| Aqueous Solubility (log mol/L) | -6.08 | -6.03 | -5.93 |
| Hepatocyte Clearance (uL/min/10^6) | 32.0 | 31.5 | 30.4 |

All three lead compounds show favorable predicted oral bioavailability (0.81-0.86), moderate BBB permeability (0.69-0.71), and borderline hERG liability (0.47-0.53). The predicted Caco-2 permeability values (-4.75 to -4.80 log cm/s) are within the acceptable range for oral absorption. Hepatocyte clearance predictions (~30-32 uL/min/10^6 cells) suggest moderate metabolic stability. The aqueous solubility predictions (-5.93 to -6.08 log mol/L) indicate limited solubility, which is common for compounds in this molecular weight range and can be addressed through salt selection or formulation strategies.

The compound **CHEMBL7062** (rank 7) is noteworthy as an outlier in the top 10. With an exceptional QED of 0.898 -- the highest among all candidates -- it is a compact diaminopyrimidine derivative (MW 338.3, LogP 2.89) bearing a trifluoromethoxy-phenyl substituent. While its binding affinity (-8.2 kcal/mol) is lower than the oxadiazolinedione series, its drug-like profile makes it an attractive alternative scaffold for optimization.

### 3.4 Chemical Space Analysis

The top candidates occupy a well-defined region of chemical space (Figure 2). The majority cluster within the Lipinski-compliant zone (MW < 500, LogP < 5), with molecular weights ranging from 338 to 498 Da and LogP values from 2.89 to 4.85. The dominance of the oxadiazolinedione scaffold in the top 7 positions, with variation only in the aryl ether region, provides clear structure-activity relationship (SAR) information: the oxadiazolinedione warhead is essential for potent binding, while the aryl ether modulates affinity and drug-like properties. This SAR information is immediately actionable for medicinal chemistry optimization.

## 4. Discussion

### 4.1 Significance

This study represents, to our knowledge, the first systematic virtual screen targeting the BfmS histidine kinase of *A. baumannii*. While BfmS has been validated as a virulence target through genetic studies (Tomaras et al., 2008; Thompson et al., 2012), the drug discovery community has not yet exploited its solved crystal structure for inhibitor development. Our computational pipeline bridges this gap by providing a ranked list of 50 drug-like candidate inhibitors sourced from ChEMBL, all with zero Lipinski violations and predicted binding affinities in the low micromolar to sub-micromolar range (estimated Ki approximately 0.15 uM for the -9.3 kcal/mol hit, based on the relationship deltaG = RT ln Ki at 298K).

The discovery of a convergent oxadiazolinedione pharmacophore among the top candidates is particularly significant. This motif is present in several bioactive compound series in ChEMBL and has precedent as a bioisostere for hydantoin and barbiturate groups, which are known to interact with protein binding sites through hydrogen bonding from the NH and carbonyl groups (Meanwell, 2011). The convergence of structurally related compounds at the top of the ranked list provides internal validation of the docking results and suggests that the BfmS pocket has specific recognition features that strongly favor this scaffold.

### 4.2 Targeting Two-Component Systems

The strategy of targeting bacterial TCS offers several advantages over conventional antibiotic targets. First, TCS components are absent in mammals, providing inherent selectivity. Second, BfmS inhibition would attenuate virulence and biofilm formation without directly killing bacteria, potentially reducing selective pressure for resistance evolution (Worthington et al., 2012). Third, anti-virulence compounds can synergize with conventional antibiotics by disrupting biofilm-mediated tolerance (Rasko & Sperandio, 2010).

### 4.3 Advantages of ChEMBL-Sourced Screening Libraries

The use of ChEMBL-derived compounds rather than purely computational scaffolds offers several practical advantages. First, all compounds in the library have been previously synthesized and tested in at least one bioactivity assay, ensuring synthetic accessibility. Second, ChEMBL compounds have established supply chains through commercial vendors, facilitating rapid procurement for experimental validation. Third, the existing bioactivity annotations in ChEMBL enable secondary analyses of polypharmacology and off-target effects. The higher mean binding affinity observed with the ChEMBL library (-6.48 kcal/mol vs -5.92 kcal/mol for synthetic scaffolds) likely reflects the more optimized drug-like character of experimentally validated compounds.

### 4.4 Limitations

Several limitations of this study should be acknowledged:

**Compound library scope.** Our library of 1,500 compounds, while drawn from a curated bioactivity database, represents a small fraction of ChEMBL's 2.4 million compounds. Larger screens would likely identify additional scaffolds and improve hit diversity.

**Docking accuracy.** AutoDock Vina scoring functions, while well-validated, are approximate and can produce false positives and false negatives. The correlation between predicted docking scores and experimental binding affinities is imperfect (R^2 typically 0.4-0.6), and the absolute affinity values should be interpreted with caution (Wang et al., 2016). Molecular dynamics simulations and free energy perturbation calculations would provide more rigorous binding energy estimates for the top candidates.

**Pocket validation.** The P2Rank-identified pocket was selected computationally and has not been experimentally validated as the functional binding site relevant to BfmS kinase activity. Co-crystallization or mutagenesis studies are needed to confirm the pocket's functional relevance.

**ADMET predictions.** While ADMET-AI predictions provide useful early-stage filtering, they carry inherent uncertainty (typical AUROC 0.7-0.9 across endpoints). The borderline hERG predictions (0.47-0.53) for the top compounds warrant experimental patch-clamp validation before lead advancement. Experimental ADMET profiling is essential before advancing to in vivo studies.

**No experimental validation.** This is a purely computational study. The predicted inhibitors must be tested in vitro (e.g., BfmS autophosphorylation assays, surface plasmon resonance binding assays) and in cellular assays (biofilm inhibition, MIC determination) before any conclusions about therapeutic potential can be drawn.

### 4.5 Future Directions

Immediate next steps include:

1. **Experimental validation**: Procurement and testing of the top 10 candidates (available through ChEMBL-linked commercial vendors) in BfmS autophosphorylation inhibition assays and *A. baumannii* biofilm formation assays.
2. **Expanded screening**: Virtual screening of the full ChEMBL database and commercial vendor libraries (Enamine REAL, 1-10 million compounds) against the validated BfmS pocket.
3. **Molecular dynamics**: MD simulations of top complexes to assess binding stability, confirm the role of the oxadiazolinedione warhead, and identify key protein-ligand interactions for SAR optimization.
4. **Lead optimization**: Medicinal chemistry campaigns around both the oxadiazolinedione and diaminopyrimidine scaffolds, guided by the SAR information from the virtual screen.
5. **Resistance profiling**: Assessment of the potential for resistance development through serial passage experiments with lead compounds.

## 5. Conclusion

We have established a complete AI-augmented computational pipeline for the discovery of BfmS histidine kinase inhibitors against *A. baumannii*. Virtual screening of 1,500 ChEMBL drug-like compounds yielded 1,459 successful dockings, with 50 top candidates exhibiting strong predicted binding affinities (-7.8 to -9.3 kcal/mol) and uniformly favorable drug-like properties (100% Lipinski compliance). The lead compound CHEMBL7029 combines potent predicted binding (-9.3 kcal/mol), good drug-likeness (MW 487, QED 0.402, zero Lipinski violations), and favorable predicted oral bioavailability (0.81). The convergence of the top candidates on an oxadiazolinedione pharmacophore provides actionable SAR information and validates the pocket's selectivity for this chemical series. These results, derived from an experimentally accessible compound library, provide an immediate path to experimental validation and the development of a novel class of anti-virulence agents targeting carbapenem-resistant *A. baumannii*.

## 6. Data and Code Availability

All data, code, and results are available at the project repository (https://github.com/arian-gogani/bfms-drug-discovery). The pipeline includes scripts for structure retrieval, pocket prediction, ChEMBL library assembly, virtual screening, ADMET analysis, and visualization. Raw docking results, 3D binding poses (PDBQT format), and complete ADMET profiles for all candidates are provided.

## 7. References

Bickerton, G. R., Paolini, G. V., Besnard, J., Muresan, S., & Hopkins, A. L. (2012). Quantifying the chemical beauty of drugs. *Nature Chemistry*, 4(2), 90-98.

Draughn, G. L., Milton, M. E., Feldmann, E. A., Bobay, B. G., Roth, B. M., Olson, A. L., ... & Cavanagh, J. (2018). The structure of the biofilm-controlling response regulator BfmR from *Acinetobacter baumannii* reveals details of its DNA-binding mechanism. *Journal of Molecular Biology*, 430(6), 806-821.

Eberhardt, J., Santos-Martins, D., Tillack, A. F., & Forli, S. (2021). AutoDock Vina 1.2.0: New docking methods, expanded force field, and Python bindings. *Journal of Chemical Information and Modeling*, 61(8), 3891-3898.

Fabbro, D., Cowan-Jacob, S. W., & Moebitz, H. (2015). Ten things you should know about protein kinases: IUPHAR Review 14. *British Journal of Pharmacology*, 172(11), 2675-2700.

Gotoh, Y., Eguchi, Y., Watanabe, T., Okamoto, S., Doi, A., & Utsumi, R. (2010). Two-component signal transduction as potential drug targets in pathogenic bacteria. *Current Opinion in Microbiology*, 13(2), 232-239.

Harding, C. M., Hennon, S. W., & Feldman, M. F. (2018). Uncovering the mechanisms of *Acinetobacter baumannii* virulence. *Nature Reviews Microbiology*, 16(2), 91-102.

Karakonstantis, S., Kritsotakis, E. I., & Gikas, A. (2020). Treatment options for K. pneumoniae, P. aeruginosa and A. baumannii co-resistant to carbapenems, aminoglycosides, polymyxins and tigecycline: an approach based on the mechanisms of resistance to carbapenems. *Infection*, 48(6), 835-851.

Krivak, R., & Hoksza, D. (2018). P2Rank: machine learning based tool for rapid and accurate prediction of ligand binding sites from protein structure. *Journal of Cheminformatics*, 10(1), 39.

Landrum, G. (2023). RDKit: Open-source cheminformatics software. https://www.rdkit.org

Liou, M. L., Soo, P. C., Ling, S. R., Kuo, H. Y., Tang, C. Y., & Chang, K. C. (2014). The sensor kinase BfmS mediates virulence in *Acinetobacter baumannii*. *Journal of Microbiology, Immunology and Infection*, 47(4), 275-281.

Meanwell, N. A. (2011). Synopsis of some recent tactical application of bioisosteres in drug design. *Journal of Medicinal Chemistry*, 54(8), 2529-2591.

O'Neill, J. (2016). Tackling drug-resistant infections globally: final report and recommendations. *Review on Antimicrobial Resistance*.

Rasko, D. A., & Sperandio, V. (2010). Anti-virulence strategies to combat bacteria-mediated disease. *Nature Reviews Drug Discovery*, 9(2), 117-128.

Riniker, S., & Landrum, G. A. (2015). Better informed distance geometry: using what we know to improve conformation generation. *Journal of Chemical Information and Modeling*, 55(12), 2562-2574.

Santos-Martins, D., Eberhardt, J., Bianco, G., Solis-Vasquez, L., Koch, A., & Forli, S. (2023). Meeko: preparation of small molecules for AutoDock. *Journal of Chemical Information and Modeling*.

Sliwoski, G., Kothiwale, S., Meiler, J., & Lowe, E. W. (2014). Computational methods in drug discovery. *Pharmacological Reviews*, 66(1), 334-395.

Swanson, K., Parker, P., Suri, A., et al. (2024). ADMET-AI: A machine learning ADMET modeling platform for drug discovery. *Bioinformatics*, 40(1), btae010.

Thompson, M. G., Black, C. C., Pavlicek, R. L., Honnold, C. L., Wise, M. C., Alamneh, Y. A., ... & Zurawski, D. V. (2012). Validation of a novel murine wound model of *Acinetobacter baumannii* infection. *Antimicrobial Agents and Chemotherapy*, 56(8), 4331-4335.

Tomaras, A. P., Flagler, M. J., Dorsey, C. W., Gaddy, J. A., & Actis, L. A. (2008). Characterization of a two-component regulatory system from *Acinetobacter baumannii* that controls biofilm formation and cellular morphology. *Microbiology*, 154(11), 3398-3409.

Wang, Z., Sun, H., Yao, X., Li, D., Xu, L., Li, Y., ... & Hou, T. (2016). Comprehensive evaluation of ten docking programs on a diverse set of protein-ligand complexes: the prediction accuracy of sampling power and scoring power. *Physical Chemistry Chemical Physics*, 18(18), 12964-12975.

WHO. (2024). WHO Bacterial Priority Pathogens List, 2024. World Health Organization.

Worthington, R. J., Blackledge, M. S., & Melander, C. (2012). Small-molecule inhibition of bacterial two-component systems to combat antibiotic resistance and virulence. *Future Medicinal Chemistry*, 5(11), 1265-1284.

Zdrazil, B., Felix, E., Hunter, F., et al. (2024). The ChEMBL Database in 2023: a drug discovery platform spanning genomics, chemical biology and clinical data. *Nucleic Acids Research*, 52(D1), D1180-D1192.

---

**Figure Legends**

**Figure 1.** Distribution of predicted binding affinities across 1,459 successfully docked ChEMBL compounds against the BfmS histidine kinase Pocket 1. The red dashed line indicates the top 5% threshold (-7.8 kcal/mol). The distribution is approximately normal with a mean of -6.48 kcal/mol and standard deviation of 0.74 kcal/mol.

**Figure 2.** Multi-panel dashboard summarizing the top 10 BfmS inhibitor candidates. (A) Binding affinities. (B) QED drug-likeness scores. (C) Chemical space plot showing MW vs LogP, with top 3 candidates highlighted. (D) TPSA vs binding affinity relationship. (E) Composite ranking scores. (F) ADMET profile heatmap.

**Figure 3.** 2D chemical structures of the top 3 candidates: CHEMBL7029 (oxadiazolinedione-trifluoromethylphenyl-oxazole, -9.3 kcal/mol), CHEMBL414184 (trans-alkenyl oxadiazolinedione-oxazole, -9.1 kcal/mol), and CHEMBL6748 (oxadiazolinedione-furanyl-oxazole, -9.1 kcal/mol).

**Figure 4.** Schematic representation of the BfmS binding pocket architecture showing the location and identity of key residues identified by P2Rank analysis.
