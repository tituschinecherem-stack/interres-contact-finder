import MDAnalysis as mda
import numpy as np
import pandas as pd
import sys

def interacting_residues(pdb_file, ligand_resname, cutoff):
    """Return residues within cutoff distance as a DataFrame.
    Columns: resname, resid
    """

    u = mda.Universe(pdb_file)

    ligand = u.select_atoms(f"resname {ligand_resname}")
    if len(ligand) == 0:
        raise ValueError(f"Ligand '{ligand_resname}' not found.")

    protein = u.select_atoms(f"not resname {ligand_resname}")

    rows = []

    for res in protein.residues:
        res_atoms = res.atoms.select_atoms("not name H*")
        if len(res_atoms) == 0:
            continue

        dists = mda.lib.distances.distance_array(
            res_atoms.positions,
            ligand.positions
        )

        if np.min(dists) <= cutoff:
            rows.append([res.resname, res.resid])

    return pd.DataFrame(rows, columns=["resname", "resid"])


if __name__ == "__main__":
    # Expect: python Inter-Res.py pdbfile ligand cutoff
    if len(sys.argv) != 4:
        print("Usage: python Inter-Res.py <pdb_file> <ligand_resname> <cutoff>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    ligand = sys.argv[2]
    cutoff = float(sys.argv[3])

    df = interacting_residues(pdb_file, ligand, cutoff)

    df.to_csv("interacting_residues.csv", index=False)
    print("Saved interacting_residues.csv")

