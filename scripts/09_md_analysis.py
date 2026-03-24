#!/usr/bin/env python3
"""
MD Trajectory Analysis and Figure Generation.

Analyzes the 100 ns MD trajectory of CHEMBL7029 bound to BfmS Chain C:
1. Protein and ligand RMSD over time
2. Ligand RMSF (per-atom fluctuation)
3. Protein-ligand hydrogen bonds and contacts
4. Radius of gyration of binding pocket
5. Per-residue contact frequency heatmap
6. Binding energy estimation via MM/GBSA-like scoring

Generates 4 publication-ready figures at 300 DPI.

Usage: python scripts/09_md_analysis.py
"""

import os
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
MD_DIR = PROJECT / "results" / "md"
FIG_DIR = PROJECT / "results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Pocket residues from P2Rank (Chain C, renumbered in the extracted chain)
POCKET_RESIDUES_ORIG = [66, 69, 72, 98, 101, 153, 154, 155, 212, 215, 216]

# Plotting style
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.linewidth": 1.2,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
    "figure.dpi": 300,
})

COLORS = {
    "protein": "#2C3E50",
    "ligand": "#E74C3C",
    "pocket": "#3498DB",
    "hbond": "#27AE60",
    "contact": "#F39C12",
    "fill_protein": "#2C3E5030",
    "fill_ligand": "#E74C3C30",
}


def load_trajectory():
    """Load the MD trajectory using MDTraj."""
    import mdtraj as md

    topo_path = MD_DIR / "production_topology.pdb"
    traj_path = MD_DIR / "production.dcd"

    print(f"Loading topology: {topo_path}")
    print(f"Loading trajectory: {traj_path}")

    traj = md.load(str(traj_path), top=str(topo_path))
    print(f"  Frames: {traj.n_frames}")
    print(f"  Atoms: {traj.n_atoms}")
    print(f"  Time: {traj.time[0]:.1f} - {traj.time[-1]:.1f} ps")
    return traj


def identify_selections(traj):
    """Identify protein, ligand, pocket atoms."""
    topology = traj.topology

    protein_atoms = topology.select("protein")
    protein_ca = topology.select("protein and name CA")

    # Ligand = residue named UNL
    ligand_atoms = topology.select("resname UNL")
    if len(ligand_atoms) == 0:
        # Try other common names
        for name in ["LIG", "MOL", "DRG"]:
            ligand_atoms = topology.select(f"resname {name}")
            if len(ligand_atoms) > 0:
                break
    if len(ligand_atoms) == 0:
        # Last resort: non-protein, non-water, non-ion
        ligand_atoms = topology.select("not protein and not water and not (resname Na+ Cl- NA CL HOH WAT)")

    # Ligand heavy atoms (no hydrogens)
    ligand_heavy = np.array([i for i in ligand_atoms
                             if topology.atom(i).element.symbol != "H"])

    # Pocket residues
    pocket_atoms = []
    pocket_ca = []
    for res in topology.residues:
        if res.chain.index == 0 and res.resSeq in POCKET_RESIDUES_ORIG:
            for atom in res.atoms:
                pocket_atoms.append(atom.index)
                if atom.name == "CA":
                    pocket_ca.append(atom.index)
    pocket_atoms = np.array(pocket_atoms)
    pocket_ca = np.array(pocket_ca)

    print(f"  Protein atoms: {len(protein_atoms)}, CA: {len(protein_ca)}")
    print(f"  Ligand atoms: {len(ligand_atoms)}, heavy: {len(ligand_heavy)}")
    print(f"  Pocket atoms: {len(pocket_atoms)}, CA: {len(pocket_ca)}")

    return {
        "protein": protein_atoms,
        "protein_ca": protein_ca,
        "ligand": ligand_atoms,
        "ligand_heavy": ligand_heavy,
        "pocket": pocket_atoms,
        "pocket_ca": pocket_ca,
    }


def compute_rmsd(traj, selections):
    """Compute RMSD time series for protein backbone, ligand, and pocket."""
    import mdtraj as md

    time_ns = traj.time / 1000.0  # ps -> ns

    # Protein backbone RMSD (align on backbone, measure backbone)
    backbone = traj.topology.select("backbone")
    protein_rmsd = md.rmsd(traj, traj, frame=0, atom_indices=backbone) * 10  # nm -> A

    # Ligand RMSD (align on protein backbone, measure ligand heavy atoms)
    # First superpose on protein backbone
    traj_aligned = traj.superpose(traj, frame=0, atom_indices=backbone)
    ligand_rmsd = md.rmsd(traj_aligned, traj_aligned, frame=0,
                          atom_indices=selections["ligand_heavy"]) * 10  # nm -> A

    # Pocket RMSD
    pocket_rmsd = md.rmsd(traj_aligned, traj_aligned, frame=0,
                          atom_indices=selections["pocket"]) * 10 if len(selections["pocket"]) > 0 else None

    print(f"  Protein backbone RMSD: {np.mean(protein_rmsd):.2f} +/- {np.std(protein_rmsd):.2f} A")
    print(f"  Ligand RMSD: {np.mean(ligand_rmsd):.2f} +/- {np.std(ligand_rmsd):.2f} A")
    if pocket_rmsd is not None:
        print(f"  Pocket RMSD: {np.mean(pocket_rmsd):.2f} +/- {np.std(pocket_rmsd):.2f} A")

    return time_ns, protein_rmsd, ligand_rmsd, pocket_rmsd


def compute_contacts(traj, selections):
    """Compute protein-ligand contacts over time."""
    import mdtraj as md

    time_ns = traj.time / 1000.0
    ligand_heavy = selections["ligand_heavy"]

    # Find all protein residues within 4.5 A of ligand at any point
    contact_cutoff = 0.45  # nm

    # Compute contacts per frame
    protein_residues = list(traj.topology.residues)
    protein_res_indices = {}
    for res in protein_residues:
        if res.is_protein:
            heavy = [a.index for a in res.atoms if a.element.symbol != "H"]
            if heavy:
                protein_res_indices[res.resSeq] = heavy

    # Per-residue contact frequency
    n_frames = traj.n_frames
    residue_contacts = {}

    for resSeq, res_atoms in protein_res_indices.items():
        # Compute minimum distance between this residue and ligand per frame
        pairs = []
        for ra in res_atoms:
            for la in ligand_heavy:
                pairs.append([ra, la])
        if not pairs:
            continue
        pairs = np.array(pairs)
        distances = md.compute_distances(traj, pairs)
        min_dist_per_frame = distances.min(axis=1)
        contact_fraction = np.mean(min_dist_per_frame < contact_cutoff)
        if contact_fraction > 0.1:  # At least 10% of frames
            residue_contacts[resSeq] = {
                "fraction": contact_fraction,
                "min_distances": min_dist_per_frame,
            }

    # Total contacts per frame
    total_contacts = np.zeros(n_frames)
    for resSeq, data in residue_contacts.items():
        total_contacts += (data["min_distances"] < contact_cutoff).astype(float)

    print(f"  Residues with >10% contact: {len(residue_contacts)}")
    print(f"  Mean contacts per frame: {np.mean(total_contacts):.1f}")

    # Sort by contact frequency
    sorted_contacts = sorted(residue_contacts.items(), key=lambda x: x[1]["fraction"], reverse=True)
    print("  Top contacting residues:")
    for resSeq, data in sorted_contacts[:10]:
        # Get residue name
        for res in traj.topology.residues:
            if res.resSeq == resSeq and res.is_protein:
                print(f"    {res.name}{resSeq}: {data['fraction']*100:.1f}%")
                break

    return residue_contacts, total_contacts, time_ns


def compute_hbonds(traj, selections):
    """Compute protein-ligand hydrogen bonds over time."""
    import mdtraj as md

    # Use MDTraj's Baker-Hubbard criterion
    time_ns = traj.time / 1000.0
    ligand_set = set(selections["ligand"].tolist())
    protein_set = set(selections["protein"].tolist())

    hbond_counts = []
    hbond_residues = {}

    for frame_idx in range(traj.n_frames):
        frame = traj[frame_idx]
        hbonds = md.baker_hubbard(frame, freq=0.0)

        count = 0
        for hb in hbonds:
            donor, h, acceptor = hb
            donor_is_lig = donor in ligand_set
            acceptor_is_lig = acceptor in ligand_set
            donor_is_prot = donor in protein_set
            acceptor_is_prot = acceptor in protein_set

            if (donor_is_lig and acceptor_is_prot) or (donor_is_prot and acceptor_is_lig):
                count += 1
                # Identify protein residue
                prot_atom = acceptor if donor_is_lig else donor
                res = traj.topology.atom(prot_atom).residue
                key = f"{res.name}{res.resSeq}"
                if key not in hbond_residues:
                    hbond_residues[key] = 0
                hbond_residues[key] += 1

        hbond_counts.append(count)

    hbond_counts = np.array(hbond_counts)

    # Normalize hbond_residues to fractions
    n_frames = traj.n_frames
    for key in hbond_residues:
        hbond_residues[key] /= n_frames

    print(f"  Mean H-bonds per frame: {np.mean(hbond_counts):.2f}")
    print(f"  Max H-bonds in a frame: {np.max(hbond_counts)}")
    sorted_hb = sorted(hbond_residues.items(), key=lambda x: x[1], reverse=True)
    print("  H-bond partner residues:")
    for res, frac in sorted_hb[:5]:
        print(f"    {res}: {frac*100:.1f}%")

    return hbond_counts, hbond_residues, time_ns


def compute_rog(traj, selections):
    """Compute radius of gyration of the binding pocket over time."""
    import mdtraj as md

    pocket_atoms = selections["pocket"]
    if len(pocket_atoms) == 0:
        return None, None

    time_ns = traj.time / 1000.0
    pocket_traj = traj.atom_slice(pocket_atoms)
    rog = md.compute_rg(pocket_traj) * 10  # nm -> A

    print(f"  Pocket Rg: {np.mean(rog):.2f} +/- {np.std(rog):.2f} A")
    return rog, time_ns


def compute_ligand_rmsf(traj, selections):
    """Compute per-atom RMSF for the ligand."""
    import mdtraj as md

    backbone = traj.topology.select("backbone")
    traj_aligned = traj.superpose(traj, frame=0, atom_indices=backbone)

    ligand_heavy = selections["ligand_heavy"]
    rmsf = md.rmsf(traj_aligned, traj_aligned, frame=0, atom_indices=ligand_heavy) * 10  # A

    print(f"  Ligand RMSF range: {np.min(rmsf):.2f} - {np.max(rmsf):.2f} A")
    print(f"  Ligand mean RMSF: {np.mean(rmsf):.2f} A")
    return rmsf, ligand_heavy


def load_energy_log():
    """Load energy data from the production log CSV."""
    log_path = MD_DIR / "production_log.csv"
    if not log_path.exists():
        return None

    data = {"time_ns": [], "potential": [], "kinetic": [], "total": [],
            "temperature": [], "volume": [], "density": []}

    with open(log_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                # Time is typically in ps
                t = float(row.get('#"Time (ps)"', row.get("Time (ps)", 0)))
                data["time_ns"].append(t / 1000.0)
                data["potential"].append(float(row.get("Potential Energy (kJ/mole)", 0)))
                data["kinetic"].append(float(row.get("Kinetic Energy (kJ/mole)", 0)))
                data["total"].append(float(row.get("Total Energy (kJ/mole)", 0)))
                data["temperature"].append(float(row.get("Temperature (K)", 0)))
                data["volume"].append(float(row.get("Box Volume (nm^3)", 0)))
                data["density"].append(float(row.get("Density (g/mL)", 0)))
            except (ValueError, KeyError):
                continue

    for key in data:
        data[key] = np.array(data[key])

    return data


def generate_figure_5(time_ns, protein_rmsd, ligand_rmsd, pocket_rmsd,
                       hbond_counts, total_contacts, rog, time_ns_rog):
    """Figure 5: MD Binding Stability Dashboard (4-panel)."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle("CHEMBL7029–BfmS Molecular Dynamics: Binding Stability (100 ns)",
                 fontsize=13, fontweight="bold", y=0.98)

    # Panel A: RMSD time series
    ax = axes[0, 0]
    ax.plot(time_ns, protein_rmsd, color=COLORS["protein"], alpha=0.7, linewidth=0.5, label="Protein backbone")
    ax.plot(time_ns, ligand_rmsd, color=COLORS["ligand"], alpha=0.7, linewidth=0.5, label="Ligand (heavy atoms)")
    if pocket_rmsd is not None:
        ax.plot(time_ns, pocket_rmsd, color=COLORS["pocket"], alpha=0.7, linewidth=0.5, label="Binding pocket")
    # Rolling averages
    window = max(1, len(time_ns) // 100)
    ax.plot(time_ns, np.convolve(protein_rmsd, np.ones(window)/window, mode="same"),
            color=COLORS["protein"], linewidth=2)
    ax.plot(time_ns, np.convolve(ligand_rmsd, np.ones(window)/window, mode="same"),
            color=COLORS["ligand"], linewidth=2)
    if pocket_rmsd is not None:
        ax.plot(time_ns, np.convolve(pocket_rmsd, np.ones(window)/window, mode="same"),
                color=COLORS["pocket"], linewidth=2)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (Å)")
    ax.set_title("A. RMSD Convergence", fontweight="bold")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_xlim(0, time_ns[-1])
    ax.axhline(y=np.mean(ligand_rmsd[-len(ligand_rmsd)//4:]), color=COLORS["ligand"],
               linestyle="--", alpha=0.4, linewidth=1)

    # Panel B: H-bond and contact count
    ax = axes[0, 1]
    ax.plot(time_ns, np.convolve(total_contacts, np.ones(window)/window, mode="same"),
            color=COLORS["contact"], linewidth=1.5, label="Residue contacts", alpha=0.8)
    ax2 = ax.twinx()
    ax2.plot(time_ns, np.convolve(hbond_counts.astype(float), np.ones(window)/window, mode="same"),
             color=COLORS["hbond"], linewidth=1.5, label="H-bonds", alpha=0.8)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Protein–Ligand Contacts", color=COLORS["contact"])
    ax2.set_ylabel("Hydrogen Bonds", color=COLORS["hbond"])
    ax.set_title("B. Interaction Persistence", fontweight="bold")
    ax.set_xlim(0, time_ns[-1])
    # Combined legend
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=8, loc="upper right")

    # Panel C: Pocket radius of gyration
    ax = axes[1, 0]
    if rog is not None and time_ns_rog is not None:
        ax.plot(time_ns_rog, rog, color=COLORS["pocket"], alpha=0.3, linewidth=0.5)
        window_rog = max(1, len(time_ns_rog) // 100)
        ax.plot(time_ns_rog, np.convolve(rog, np.ones(window_rog)/window_rog, mode="same"),
                color=COLORS["pocket"], linewidth=2)
        ax.set_ylabel("Radius of Gyration (Å)")
        ax.axhline(y=np.mean(rog), color=COLORS["pocket"], linestyle="--", alpha=0.5)
        mean_rog = np.mean(rog)
        std_rog = np.std(rog)
        ax.fill_between(time_ns_rog,
                         np.full_like(time_ns_rog, mean_rog - std_rog),
                         np.full_like(time_ns_rog, mean_rog + std_rog),
                         color=COLORS["pocket"], alpha=0.1)
        ax.text(0.95, 0.95, f"Mean: {mean_rog:.1f} ± {std_rog:.1f} Å",
                transform=ax.transAxes, ha="right", va="top", fontsize=9,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    ax.set_xlabel("Time (ns)")
    ax.set_title("C. Binding Pocket Compactness", fontweight="bold")
    ax.set_xlim(0, time_ns[-1] if time_ns_rog is None else time_ns_rog[-1])

    # Panel D: RMSD distribution (violin/histogram)
    ax = axes[1, 1]
    # Last 75% of trajectory (equilibrated portion)
    eq_start = len(protein_rmsd) // 4
    ax.hist(protein_rmsd[eq_start:], bins=50, alpha=0.5, color=COLORS["protein"],
            label=f"Protein ({np.mean(protein_rmsd[eq_start:]):.2f} ± {np.std(protein_rmsd[eq_start:]):.2f} Å)",
            density=True, edgecolor="white", linewidth=0.5)
    ax.hist(ligand_rmsd[eq_start:], bins=50, alpha=0.5, color=COLORS["ligand"],
            label=f"Ligand ({np.mean(ligand_rmsd[eq_start:]):.2f} ± {np.std(ligand_rmsd[eq_start:]):.2f} Å)",
            density=True, edgecolor="white", linewidth=0.5)
    ax.set_xlabel("RMSD (Å)")
    ax.set_ylabel("Density")
    ax.set_title("D. RMSD Distribution (equilibrated)", fontweight="bold")
    ax.legend(fontsize=8)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out = FIG_DIR / "fig5_md_binding_stability.png"
    plt.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved: {out}")


def generate_figure_6(residue_contacts, hbond_residues, rmsf, ligand_heavy, traj):
    """Figure 6: Residue Interaction Map and Ligand Flexibility."""
    fig = plt.figure(figsize=(14, 6))
    gs = GridSpec(1, 3, width_ratios=[2, 1, 1], wspace=0.35)

    # Panel A: Per-residue contact frequency bar chart
    ax = fig.add_subplot(gs[0])
    sorted_contacts = sorted(residue_contacts.items(), key=lambda x: x[1]["fraction"], reverse=True)
    top_contacts = sorted_contacts[:20]

    residue_labels = []
    contact_fracs = []
    hb_fracs = []
    for resSeq, data in top_contacts:
        for res in traj.topology.residues:
            if res.resSeq == resSeq and res.is_protein:
                label = f"{res.name}{resSeq}"
                residue_labels.append(label)
                contact_fracs.append(data["fraction"])
                hb_fracs.append(hbond_residues.get(label, 0))
                break

    y_pos = np.arange(len(residue_labels))
    ax.barh(y_pos, contact_fracs, color=COLORS["contact"], alpha=0.7, label="Van der Waals contact")
    ax.barh(y_pos, hb_fracs, color=COLORS["hbond"], alpha=0.8, label="Hydrogen bond")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(residue_labels, fontsize=8)
    ax.set_xlabel("Interaction Frequency")
    ax.set_title("A. Per-Residue Interaction Frequency", fontweight="bold")
    ax.legend(fontsize=8, loc="lower right")
    ax.invert_yaxis()

    # Highlight pocket residues
    for i, label in enumerate(residue_labels):
        resnum = int("".join(c for c in label if c.isdigit()))
        if resnum in POCKET_RESIDUES_ORIG:
            ax.get_yticklabels()[i].set_fontweight("bold")
            ax.get_yticklabels()[i].set_color(COLORS["pocket"])

    # Panel B: Ligand RMSF per heavy atom
    ax = fig.add_subplot(gs[1])
    atom_indices = np.arange(len(rmsf))
    colors_rmsf = [COLORS["ligand"] if r < np.mean(rmsf) else "#C0392B" for r in rmsf]
    ax.bar(atom_indices, rmsf, color=colors_rmsf, alpha=0.8, width=0.8)
    ax.axhline(y=np.mean(rmsf), color="gray", linestyle="--", alpha=0.5, label=f"Mean: {np.mean(rmsf):.2f} Å")
    ax.set_xlabel("Ligand Heavy Atom Index")
    ax.set_ylabel("RMSF (Å)")
    ax.set_title("B. Ligand Flexibility", fontweight="bold")
    ax.legend(fontsize=8)

    # Panel C: Contact stability over simulation quarters
    ax = fig.add_subplot(gs[2])
    n_frames = len(list(residue_contacts.values())[0]["min_distances"]) if residue_contacts else 0
    if n_frames > 0:
        quarters = ["Q1\n0-25 ns", "Q2\n25-50 ns", "Q3\n50-75 ns", "Q4\n75-100 ns"]
        quarter_size = n_frames // 4
        top5_contacts = sorted_contacts[:5]

        for idx, (resSeq, data) in enumerate(top5_contacts):
            quarter_fracs = []
            for q in range(4):
                start = q * quarter_size
                end = (q + 1) * quarter_size if q < 3 else n_frames
                frac = np.mean(data["min_distances"][start:end] < 0.45)
                quarter_fracs.append(frac)
            for res in traj.topology.residues:
                if res.resSeq == resSeq and res.is_protein:
                    ax.plot(range(4), quarter_fracs, "o-", linewidth=2, markersize=6,
                            label=f"{res.name}{resSeq}", alpha=0.8)
                    break

        ax.set_xticks(range(4))
        ax.set_xticklabels(quarters, fontsize=8)
        ax.set_ylabel("Contact Frequency")
        ax.set_ylim(0, 1.05)
        ax.set_title("C. Contact Stability", fontweight="bold")
        ax.legend(fontsize=7, loc="lower left")

    plt.tight_layout()
    out = FIG_DIR / "fig6_residue_interactions.png"
    plt.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved: {out}")


def generate_figure_7(energy_data):
    """Figure 7: Simulation diagnostics (energy, temperature, density)."""
    if energy_data is None:
        print("  No energy log found, skipping Figure 7.")
        return

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("CHEMBL7029–BfmS MD Simulation Diagnostics",
                 fontsize=13, fontweight="bold", y=0.98)

    t = energy_data["time_ns"]

    # Panel A: Total energy
    ax = axes[0, 0]
    ax.plot(t, energy_data["total"] / 1000, color=COLORS["protein"], linewidth=0.5, alpha=0.5)
    window = max(1, len(t) // 50)
    ax.plot(t, np.convolve(energy_data["total"] / 1000, np.ones(window)/window, mode="same"),
            color=COLORS["protein"], linewidth=2)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Total Energy (× 10³ kJ/mol)")
    ax.set_title("A. Total Energy", fontweight="bold")

    # Panel B: Temperature
    ax = axes[0, 1]
    ax.plot(t, energy_data["temperature"], color=COLORS["ligand"], linewidth=0.5, alpha=0.3)
    ax.plot(t, np.convolve(energy_data["temperature"], np.ones(window)/window, mode="same"),
            color=COLORS["ligand"], linewidth=2)
    ax.axhline(y=300, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Temperature (K)")
    ax.set_title("B. Temperature Control", fontweight="bold")
    ax.set_ylim(280, 320)

    # Panel C: Density
    ax = axes[1, 0]
    if np.any(energy_data["density"] > 0):
        ax.plot(t, energy_data["density"], color=COLORS["pocket"], linewidth=0.5, alpha=0.3)
        ax.plot(t, np.convolve(energy_data["density"], np.ones(window)/window, mode="same"),
                color=COLORS["pocket"], linewidth=2)
        ax.axhline(y=1.0, color="gray", linestyle="--", alpha=0.5)
        ax.set_ylabel("Density (g/mL)")
    else:
        ax.plot(t, energy_data["volume"], color=COLORS["pocket"], linewidth=0.5, alpha=0.3)
        ax.set_ylabel("Box Volume (nm³)")
    ax.set_xlabel("Time (ns)")
    ax.set_title("C. System Density", fontweight="bold")

    # Panel D: Potential vs kinetic energy
    ax = axes[1, 1]
    ax.plot(t, energy_data["potential"] / 1000, color=COLORS["protein"], linewidth=0.5, alpha=0.3, label="Potential")
    ax.plot(t, energy_data["kinetic"] / 1000, color=COLORS["ligand"], linewidth=0.5, alpha=0.3, label="Kinetic")
    ax.plot(t, np.convolve(energy_data["potential"] / 1000, np.ones(window)/window, mode="same"),
            color=COLORS["protein"], linewidth=2)
    ax.plot(t, np.convolve(energy_data["kinetic"] / 1000, np.ones(window)/window, mode="same"),
            color=COLORS["ligand"], linewidth=2)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Energy (× 10³ kJ/mol)")
    ax.set_title("D. Energy Components", fontweight="bold")
    ax.legend(fontsize=8)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out = FIG_DIR / "fig7_md_diagnostics.png"
    plt.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved: {out}")


def write_summary_csv(time_ns, protein_rmsd, ligand_rmsd, pocket_rmsd,
                       hbond_counts, total_contacts, residue_contacts, hbond_residues):
    """Write summary statistics to CSV."""
    # Time series summary
    summary_path = MD_DIR / "md_summary_statistics.csv"
    eq_start = len(protein_rmsd) // 4

    with open(summary_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Metric", "Full_Mean", "Full_SD", "Equilibrated_Mean", "Equilibrated_SD"])
        w.writerow(["Protein_Backbone_RMSD_A",
                     f"{np.mean(protein_rmsd):.3f}", f"{np.std(protein_rmsd):.3f}",
                     f"{np.mean(protein_rmsd[eq_start:]):.3f}", f"{np.std(protein_rmsd[eq_start:]):.3f}"])
        w.writerow(["Ligand_RMSD_A",
                     f"{np.mean(ligand_rmsd):.3f}", f"{np.std(ligand_rmsd):.3f}",
                     f"{np.mean(ligand_rmsd[eq_start:]):.3f}", f"{np.std(ligand_rmsd[eq_start:]):.3f}"])
        if pocket_rmsd is not None:
            w.writerow(["Pocket_RMSD_A",
                         f"{np.mean(pocket_rmsd):.3f}", f"{np.std(pocket_rmsd):.3f}",
                         f"{np.mean(pocket_rmsd[eq_start:]):.3f}", f"{np.std(pocket_rmsd[eq_start:]):.3f}"])
        w.writerow(["H_Bonds_Per_Frame",
                     f"{np.mean(hbond_counts):.3f}", f"{np.std(hbond_counts):.3f}",
                     f"{np.mean(hbond_counts[eq_start:]):.3f}", f"{np.std(hbond_counts[eq_start:]):.3f}"])
        w.writerow(["Residue_Contacts_Per_Frame",
                     f"{np.mean(total_contacts):.3f}", f"{np.std(total_contacts):.3f}",
                     f"{np.mean(total_contacts[eq_start:]):.3f}", f"{np.std(total_contacts[eq_start:]):.3f}"])

    print(f"  Summary statistics saved: {summary_path}")

    # Per-residue contact table
    contacts_path = MD_DIR / "residue_contact_frequencies.csv"
    with open(contacts_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Residue", "Contact_Frequency", "HBond_Frequency", "Is_Pocket_Residue"])
        for resSeq, data in sorted(residue_contacts.items(), key=lambda x: x[1]["fraction"], reverse=True):
            resname = ""
            for res_label in hbond_residues:
                if str(resSeq) in res_label:
                    resname = res_label
                    break
            if not resname:
                resname = str(resSeq)
            hb_frac = hbond_residues.get(resname, 0)
            is_pocket = resSeq in POCKET_RESIDUES_ORIG
            w.writerow([resname, f"{data['fraction']:.4f}", f"{hb_frac:.4f}", is_pocket])

    print(f"  Residue contacts saved: {contacts_path}")


def main():
    print("=" * 70)
    print("BfmS-CHEMBL7029 MD Trajectory Analysis")
    print("=" * 70)

    print("\n[1/8] Loading trajectory...")
    traj = load_trajectory()

    print("\n[2/8] Identifying atom selections...")
    selections = identify_selections(traj)

    print("\n[3/8] Computing RMSD...")
    time_ns, protein_rmsd, ligand_rmsd, pocket_rmsd = compute_rmsd(traj, selections)

    print("\n[4/8] Computing protein-ligand contacts...")
    residue_contacts, total_contacts, _ = compute_contacts(traj, selections)

    print("\n[5/8] Computing hydrogen bonds...")
    hbond_counts, hbond_residues, _ = compute_hbonds(traj, selections)

    print("\n[6/8] Computing pocket radius of gyration...")
    rog, time_ns_rog = compute_rog(traj, selections)

    print("\n[7/8] Computing ligand RMSF...")
    rmsf, ligand_heavy = compute_ligand_rmsf(traj, selections)

    print("\n[8/8] Loading energy data...")
    energy_data = load_energy_log()

    print("\n" + "=" * 70)
    print("Generating Figures")
    print("=" * 70)

    print("\n  Figure 5: MD Binding Stability Dashboard...")
    generate_figure_5(time_ns, protein_rmsd, ligand_rmsd, pocket_rmsd,
                       hbond_counts, total_contacts, rog, time_ns_rog)

    print("  Figure 6: Residue Interaction Map...")
    generate_figure_6(residue_contacts, hbond_residues, rmsf, ligand_heavy, traj)

    print("  Figure 7: Simulation Diagnostics...")
    generate_figure_7(energy_data)

    print("\n  Writing summary data...")
    write_summary_csv(time_ns, protein_rmsd, ligand_rmsd, pocket_rmsd,
                       hbond_counts, total_contacts, residue_contacts, hbond_residues)

    print("\n" + "=" * 70)
    print("Analysis complete!")
    print(f"  Figures: {FIG_DIR}/fig5_md_binding_stability.png")
    print(f"           {FIG_DIR}/fig6_residue_interactions.png")
    print(f"           {FIG_DIR}/fig7_md_diagnostics.png")
    print(f"  Data:    {MD_DIR}/md_summary_statistics.csv")
    print(f"           {MD_DIR}/residue_contact_frequencies.csv")
    print("=" * 70)


if __name__ == "__main__":
    main()
