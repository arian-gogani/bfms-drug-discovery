#!/usr/bin/env python3
"""
Molecular Dynamics Simulation of CHEMBL7029 bound to BfmS (Chain C).

100 ns NPT production MD using OpenMM with Amber14 (protein) + MMFF94 (ligand).
"""

import os
import sys
import re
import time
import numpy as np
from pathlib import Path

SIM_TIME_NS = 100
TIMESTEP_FS = 1.0  # 1 fs, no HMR or constraints needed for ligand
TEMPERATURE_K = 300.0
PRESSURE_ATM = 1.0
PADDING_NM = 1.2
IONIC_STRENGTH_M = 0.15
SAVE_INTERVAL_PS = 100
LOG_INTERVAL_PS = 10
CHECKPOINT_INTERVAL_PS = 1000

PROJECT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT / "data"
RESULTS_DIR = PROJECT / "results"
MD_DIR = RESULTS_DIR / "md"
MD_DIR.mkdir(parents=True, exist_ok=True)

SMILES = "C/C(=C\\Cn1oc(=O)[nH]c1=O)c1cccc(OCc2nc(-c3ccc(C(F)(F)F)cc3)oc2C)c1"
POSE_FILE = RESULTS_DIR / "poses" / "CHEMBL7029_pose.pdbqt"
PROTEIN_FILE = DATA_DIR / "3KLN_clean.pdb"

RESUME = "--resume" in sys.argv


def parse_pdbqt_coords(pdbqt_path, model=1):
    coords = []
    current_model = 0
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                current_model = int(line.split()[1])
            if current_model != model and current_model > 0:
                continue
            if current_model == model and line.startswith("ENDMDL"):
                break
            if line.startswith("ATOM") or line.startswith("HETATM"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    return np.array(coords)


def parse_smiles_idx_mapping(pdbqt_path, model=1):
    mapping = {}
    current_model = 0
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                current_model = int(line.split()[1])
            if current_model != model and current_model > 0:
                continue
            if "SMILES IDX" in line:
                parts = line.split("SMILES IDX")[1].split()
                for i in range(0, len(parts), 2):
                    smiles_idx = int(parts[i]) - 1
                    pdbqt_num = int(parts[i+1]) - 1
                    mapping[smiles_idx] = pdbqt_num
    return mapping


def create_ligand(rdkit_mol=None):
    """Create RDKit molecule with docked coordinates."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    if rdkit_mol is None:
        mol = Chem.MolFromSmiles(SMILES)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)
    else:
        mol = rdkit_mol

    pdbqt_coords = parse_pdbqt_coords(POSE_FILE)
    smiles_to_pdbqt = parse_smiles_idx_mapping(POSE_FILE)

    conf = mol.GetConformer()
    heavy_idx = 0
    mapped = set()
    for atom_idx in range(mol.GetNumAtoms()):
        if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 1:
            continue
        if heavy_idx in smiles_to_pdbqt:
            pdbqt_idx = smiles_to_pdbqt[heavy_idx]
            if pdbqt_idx < len(pdbqt_coords):
                conf.SetAtomPosition(atom_idx, pdbqt_coords[pdbqt_idx].tolist())
                mapped.add(atom_idx)
        heavy_idx += 1

    ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
    if ff:
        for idx in mapped:
            ff.AddFixedPoint(idx)
        ff.Minimize(maxIts=200)

    print(f"  Ligand: {mol.GetNumAtoms()} atoms, {len(mapped)} heavy atoms mapped from docked pose")
    return mol


def extract_chain_c():
    chain_c_path = MD_DIR / "chain_c.pdb"
    lines = []
    with open(PROTEIN_FILE) as f:
        for line in f:
            if line.startswith("ATOM") and line[21] == "C":
                lines.append(line)
    with open(chain_c_path, "w") as f:
        for line in lines:
            f.write(line)
        if lines:
            last = lines[-1]
            serial = int(last[6:11]) + 1
            f.write(f"TER   {serial:5d}      {last[17:20]} {last[21]}{last[22:26]}\n")
        f.write("END\n")
    residues = set(int(l[22:26]) for l in lines)
    print(f"  Chain C: {len(residues)} residues")
    return chain_c_path


def get_ligand_coords_nm(rdkit_mol):
    """Get ligand coordinates in nm."""
    conf = rdkit_mol.GetConformer()
    coords = []
    for i in range(rdkit_mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        coords.append([pos.x / 10.0, pos.y / 10.0, pos.z / 10.0])
    return np.array(coords)


def add_ligand_forces(system, rdkit_mol, mmff_props, lig_indices):
    """Add MMFF94 bonded/nonbonded parameters for ligand atoms to the system."""
    import openmm as mm
    import openmm.unit as unit

    mol = rdkit_mol
    lj = {
        "H": (0.1069, 0.0657), "C": (0.1908, 0.3598),
        "N": (0.1824, 0.7113), "O": (0.1661, 0.8786),
        "F": (0.1750, 0.2552), "S": (0.2000, 1.0460),
    }

    nb = bond_f = angle_f = torsion_f = None
    for force in system.getForces():
        if isinstance(force, mm.NonbondedForce): nb = force
        elif isinstance(force, mm.HarmonicBondForce): bond_f = force
        elif isinstance(force, mm.HarmonicAngleForce): angle_f = force
        elif isinstance(force, mm.PeriodicTorsionForce): torsion_f = force

    # Set nonbonded params
    for ri, ti in enumerate(lig_indices):
        atom = mol.GetAtomWithIdx(ri)
        charge = mmff_props.GetMMFFPartialCharge(ri)
        sigma, eps = lj.get(atom.GetSymbol(), (0.1908, 0.3598))
        nb.setParticleParameters(ti, charge, sigma, eps)

    # MMFF94 unit conversions (Halgren 1996):
    # Bonds: E = 143.9325 * kb/2 * Δr² kcal/mol; kb in md/Å, Δr in Å
    #   OpenMM k (kJ/mol/nm²) = 143.9325 * kb * 4.184 * 100
    BOND_K_CONV = 143.9325 * 4.184 * 100
    # Angles: E = 0.043844 * ka/2 * Δθ² kcal/mol; ka in md·Å/rad², Δθ in DEGREES
    #   Convert deg->rad: factor of (180/π)²
    #   OpenMM k (kJ/mol/rad²) = 0.043844 * (180/π)² * ka * 4.184
    ANGLE_K_CONV = 0.043844 * (180.0 / np.pi) ** 2 * 4.184
    # Torsions: V in kcal/mol, E = V/2*(1+cos φ) or V/2*(1-cos 2φ) etc.
    #   OpenMM k (kJ/mol) = |V|/2 * 4.184
    TORSION_K_CONV = 4.184 / 2.0

    # Bonds (skip H-X bonds — those will be constrained instead)
    excluded_pairs = set()
    h_bond_pairs = set()
    for bond in mol.GetBonds():
        ri, rj = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        ti, tj = lig_indices[ri], lig_indices[rj]
        a1 = mol.GetAtomWithIdx(ri)
        a2 = mol.GetAtomWithIdx(rj)
        is_hbond = (a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1)
        params = mmff_props.GetMMFFBondStretchParams(mol, ri, rj)
        if params:
            _, kb, r0_a = params
            bond_f.addBond(ti, tj, r0_a / 10.0, kb * BOND_K_CONV)
        pair = (min(ti, tj), max(ti, tj))
        if pair not in excluded_pairs:
            nb.addException(ti, tj, 0, 1, 0)
            excluded_pairs.add(pair)

    # Angles
    for i in range(mol.GetNumAtoms()):
        nbrs = [n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors()]
        for ni in range(len(nbrs)):
            for nj in range(ni + 1, len(nbrs)):
                a1, a2, a3 = nbrs[ni], i, nbrs[nj]
                params = mmff_props.GetMMFFAngleBendParams(mol, a1, a2, a3)
                if params:
                    _, ka, theta0 = params
                    angle_f.addAngle(lig_indices[a1], lig_indices[a2], lig_indices[a3],
                                     theta0 * np.pi / 180.0, ka * ANGLE_K_CONV)
                # 1-3 exclusion
                pair = (min(lig_indices[nbrs[ni]], lig_indices[nbrs[nj]]),
                        max(lig_indices[nbrs[ni]], lig_indices[nbrs[nj]]))
                if pair not in excluded_pairs:
                    nb.addException(pair[0], pair[1], 0, 1, 0)
                    excluded_pairs.add(pair)

    # Torsions + 1-4 interactions
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        for ni in mol.GetAtomWithIdx(i).GetNeighbors():
            if ni.GetIdx() == j: continue
            for nj in mol.GetAtomWithIdx(j).GetNeighbors():
                if nj.GetIdx() == i: continue
                a1, a2, a3, a4 = ni.GetIdx(), i, j, nj.GetIdx()
                params = mmff_props.GetMMFFTorsionParams(mol, a1, a2, a3, a4)
                if params:
                    _, v1, v2, v3 = params
                    # MMFF: V1/2*(1+cosφ), V2/2*(1-cos2φ), V3/2*(1+cos3φ)
                    # V2 has opposite phase convention (1-cos vs 1+cos)
                    for per, v, base_phase in [(1, v1, 0.0), (2, v2, np.pi), (3, v3, 0.0)]:
                        if abs(v) > 1e-6:
                            phase = base_phase if v >= 0 else (np.pi - base_phase)
                            torsion_f.addTorsion(
                                lig_indices[a1], lig_indices[a2],
                                lig_indices[a3], lig_indices[a4],
                                per, phase, abs(v) * TORSION_K_CONV)
                # 1-4 scaled
                t1, t4 = lig_indices[a1], lig_indices[a4]
                pair = (min(t1, t4), max(t1, t4))
                if pair not in excluded_pairs:
                    q1 = mmff_props.GetMMFFPartialCharge(a1)
                    q4 = mmff_props.GetMMFFPartialCharge(a4)
                    e1 = mol.GetAtomWithIdx(a1).GetSymbol()
                    e4 = mol.GetAtomWithIdx(a4).GetSymbol()
                    s1, ep1 = lj.get(e1, (0.1908, 0.3598))
                    s4, ep4 = lj.get(e4, (0.1908, 0.3598))
                    nb.addException(t1, t4, q1 * q4 * 0.8333, (s1 + s4) / 2, (ep1 * ep4) ** 0.5 * 0.5)
                    excluded_pairs.add(pair)

    print(f"  Ligand forces: {len(excluded_pairs)} exclusions")


def setup_system():
    """Build the complete protein-ligand-solvent system."""
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    from rdkit.Chem import AllChem
    from pdbfixer import PDBFixer

    print("\n[1] Preparing ligand...")
    rdkit_mol = create_ligand()
    mmff_props = AllChem.MMFFGetMoleculeProperties(rdkit_mol)
    lig_coords_nm = get_ligand_coords_nm(rdkit_mol)

    print("\n[2] Extracting and fixing Chain C...")
    chain_c_path = extract_chain_c()
    fixer = PDBFixer(filename=str(chain_c_path))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixed_path = MD_DIR / "chain_c_fixed.pdb"
    with open(fixed_path, "w") as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)

    print("\n[3] Building protein + hydrogens...")
    forcefield = app.ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    protein_pdb = app.PDBFile(str(fixed_path))
    modeller = app.Modeller(protein_pdb.topology, protein_pdb.positions)
    modeller.addHydrogens(forcefield, pH=7.0)
    print(f"  Protein with H: {modeller.topology.getNumAtoms()} atoms")

    print("\n[4] Solvating...")
    modeller.addSolvent(
        forcefield, model="tip3p",
        padding=PADDING_NM * unit.nanometers,
        ionicStrength=IONIC_STRENGTH_M * unit.molar,
        positiveIon="Na+", negativeIon="Cl-"
    )
    print(f"  Solvated: {modeller.topology.getNumAtoms()} atoms")

    print("\n[5] Removing waters clashing with ligand...")
    pos_nm = np.array(modeller.positions.value_in_unit(unit.nanometers))
    clash_residues = []
    for residue in modeller.topology.residues():
        if residue.name == "HOH":
            for atom in residue.atoms():
                dists = np.linalg.norm(lig_coords_nm - pos_nm[atom.index], axis=1)
                if np.min(dists) < 0.25:  # 2.5 A
                    clash_residues.append(residue)
                    break
    if clash_residues:
        modeller.delete(clash_residues)
        print(f"  Removed {len(clash_residues)} waters, now {modeller.topology.getNumAtoms()} atoms")

    print("\n[6] Creating OpenMM system + injecting ligand...")
    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds,
        hydrogenMass=None
    )

    n_protein_water = system.getNumParticles()

    # Add ligand particles
    lig_indices = []
    for i in range(rdkit_mol.GetNumAtoms()):
        atom = rdkit_mol.GetAtomWithIdx(i)
        mass = atom.GetMass()
        system.addParticle(mass * unit.amu)
        lig_indices.append(n_protein_water + i)
        # Add dummy nonbonded params (will be overridden)
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                force.addParticle(0, 0.1, 0)
                break

    # Set proper ligand forces
    add_ligand_forces(system, rdkit_mol, mmff_props, lig_indices)

    # Build combined positions
    prot_pos = modeller.positions.value_in_unit(unit.nanometers)
    all_pos = [mm.Vec3(*p) for p in prot_pos]
    for c in lig_coords_nm:
        all_pos.append(mm.Vec3(*c))
    positions = all_pos * unit.nanometers

    # Add ligand to topology
    chain = modeller.topology.addChain("B")
    res = modeller.topology.addResidue("UNL", chain)
    elem_counts = {}
    lig_topo_atoms = []
    for i in range(rdkit_mol.GetNumAtoms()):
        atom = rdkit_mol.GetAtomWithIdx(i)
        elem = atom.GetSymbol()
        elem_counts[elem] = elem_counts.get(elem, 0) + 1
        name = f"{elem}{elem_counts[elem]}"
        a = modeller.topology.addAtom(name, app.element.Element.getBySymbol(elem), res)
        lig_topo_atoms.append(a)
    for bond in rdkit_mol.GetBonds():
        modeller.topology.addBond(
            lig_topo_atoms[bond.GetBeginAtomIdx()],
            lig_topo_atoms[bond.GetEndAtomIdx()])

    # Barostat
    system.addForce(mm.MonteCarloBarostat(
        PRESSURE_ATM * unit.atmospheres, TEMPERATURE_K * unit.kelvin, 25))

    topology = modeller.topology
    n_total = system.getNumParticles()
    print(f"  Final system: {n_total} particles")
    print(f"  Ligand indices: {lig_indices[0]}-{lig_indices[-1]}")

    # Save topology PDB
    with open(MD_DIR / "solvated_system.pdb", "w") as f:
        app.PDBFile.writeFile(topology, positions, f)

    # Save system XML for resume
    with open(MD_DIR / "system.xml", "w") as f:
        f.write(mm.XmlSerializer.serialize(system))

    return system, topology, positions


def get_platform():
    import openmm as mm
    for pname in ["OpenCL", "CPU", "Reference"]:
        try:
            platform = mm.Platform.getPlatformByName(pname)
            props = {"Precision": "mixed"} if pname == "OpenCL" else {}
            test_sys = mm.System()
            test_sys.addParticle(1.0)
            test_int = mm.VerletIntegrator(0.001)
            ctx = mm.Context(test_sys, test_int, platform, props)
            del ctx, test_int, test_sys
            return platform, props, pname
        except Exception:
            continue
    return mm.Platform.getPlatformByName("Reference"), {}, "Reference"


def main():
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit

    print("=" * 70)
    print("BfmS-CHEMBL7029 Molecular Dynamics Simulation (100 ns)")
    print("=" * 70)

    platform, props, pname = get_platform()
    print(f"Platform: {pname}")

    if RESUME and (MD_DIR / "checkpoint.chk").exists():
        print("\nResuming from checkpoint...")
        with open(MD_DIR / "system.xml") as f:
            system = mm.XmlSerializer.deserialize(f.read())
        topo_pdb = app.PDBFile(str(MD_DIR / "production_topology.pdb"))
        topology = topo_pdb.topology
        integrator = mm.LangevinMiddleIntegrator(
            TEMPERATURE_K * unit.kelvin, 1 / unit.picoseconds, TIMESTEP_FS * unit.femtoseconds)
        sim = app.Simulation(topology, system, integrator, platform, props)
        sim.loadCheckpoint(str(MD_DIR / "checkpoint.chk"))
    else:
        system, topology, positions = setup_system()
        integrator = mm.LangevinMiddleIntegrator(
            TEMPERATURE_K * unit.kelvin, 1 / unit.picoseconds, TIMESTEP_FS * unit.femtoseconds)
        sim = app.Simulation(topology, system, integrator, platform, props)
        sim.context.setPositions(positions)

        # Minimize
        print("\nMinimizing energy...")
        state = sim.context.getState(getEnergy=True)
        print(f"  Before: {state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole):.0f} kJ/mol")
        sim.minimizeEnergy(maxIterations=5000)
        state = sim.context.getState(getEnergy=True)
        print(f"  After:  {state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole):.0f} kJ/mol")

        # NVT equilibration
        print("\nNVT equilibration (100 ps)...")
        sim.context.setVelocitiesToTemperature(TEMPERATURE_K * unit.kelvin)
        sim.step(int(100e3 / TIMESTEP_FS))  # 100 ps
        state = sim.context.getState(getEnergy=True)
        print(f"  Energy: {state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole):.0f} kJ/mol")

        # NPT equilibration
        print("NPT equilibration (500 ps)...")
        sim.step(int(500e3 / TIMESTEP_FS))  # 500 ps
        state = sim.context.getState(getEnergy=True)
        print(f"  Energy: {state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole):.0f} kJ/mol")

    # Save reference topology
    state = sim.context.getState(getPositions=True)
    with open(MD_DIR / "production_topology.pdb", "w") as f:
        app.PDBFile.writeFile(topology, state.getPositions(), f)

    # Production reporters
    timestep = TIMESTEP_FS * unit.femtoseconds
    total_steps = int(SIM_TIME_NS * 1e6 / TIMESTEP_FS)
    save_steps = int(SAVE_INTERVAL_PS * 1e3 / TIMESTEP_FS)
    log_steps = int(LOG_INTERVAL_PS * 1e3 / TIMESTEP_FS)
    chk_steps = int(CHECKPOINT_INTERVAL_PS * 1e3 / TIMESTEP_FS)

    sim.reporters.append(app.DCDReporter(str(MD_DIR / "production.dcd"), save_steps))
    sim.reporters.append(app.StateDataReporter(
        str(MD_DIR / "production_log.csv"), log_steps,
        time=True, potentialEnergy=True, kineticEnergy=True,
        totalEnergy=True, temperature=True, volume=True,
        density=True, speed=True, separator=","))
    sim.reporters.append(app.CheckpointReporter(str(MD_DIR / "checkpoint.chk"), chk_steps))
    sim.reporters.append(app.StateDataReporter(
        sys.stdout, chk_steps,
        time=True, potentialEnergy=True, temperature=True,
        speed=True, remainingTime=True, totalSteps=total_steps))

    print(f"\nProduction: {SIM_TIME_NS} ns ({total_steps:,} steps, {TIMESTEP_FS} fs timestep)")
    t0 = time.time()
    sim.step(total_steps)
    elapsed = time.time() - t0

    print(f"\nDone! {elapsed/3600:.1f} hours, {SIM_TIME_NS/(elapsed/86400):.1f} ns/day")

    # Save final state
    sim.saveCheckpoint(str(MD_DIR / "final_checkpoint.chk"))
    state = sim.context.getState(getPositions=True)
    with open(MD_DIR / "final_state.pdb", "w") as f:
        app.PDBFile.writeFile(topology, state.getPositions(), f)


if __name__ == "__main__":
    main()
