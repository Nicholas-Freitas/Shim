### For creating specific software-specific filetypes from StandardMolecules or ShimStructures
from shim.structure import StandardMolecule
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def write_reduce_db(standard_mol: StandardMolecule,
                      residue_name: str,
                      chemical_name: str,
                      output_path: str | Path) -> None:
    """
    Write a PDB-style CONECT block for an RDKit molecule whose atoms have an
    'atomName' property, listing heavy atoms first and hydrogens last.
    """
    # Get the RDKit molecule from the StandardMolecule
    mol = standard_mol.std_mol

    atom_count = mol.GetNumAtoms()
    formula = rdMolDescriptors.CalcMolFormula(mol)

    heavy_lines: list[str] = []
    hydrogen_lines: list[str] = []

    for atom in mol.GetAtoms():
        name = atom.GetProp("atomName")
        neighbour_names = [n.GetProp("atomName") for n in atom.GetNeighbors()]

        # Build the fixed-width CONECT line
        line = (
            "CONECT      "          # cols 1-11
            f"{name:<6}"            # cols 12-17 (left-aligned)
            f"{len(neighbour_names):>2} "  # cols 18-20
            + "".join(f"{nbr:<5}" for nbr in neighbour_names)  # 21+
        )

        # Route to heavy-atom or hydrogen list
        if atom.GetAtomicNum() == 1 or name.startswith("H"):
            hydrogen_lines.append(line)
        else:
            heavy_lines.append(line)

    conect_lines = heavy_lines + hydrogen_lines

    # ---------- write file ----------
    output_path = Path(output_path)
    with output_path.open("w") as fh:
        fh.write(f"RESIDUE   {residue_name:<3}{atom_count:>7}\n")
        fh.writelines(l + "\n" for l in conect_lines)
        fh.write("END   \n")
        fh.write(f"HET    {residue_name:<5}{atom_count:>13}\n")
        fh.write(f"HETNAM     {residue_name} {chemical_name}\n")
        fh.write(f"FORMUL      {residue_name}    {formula}\n")

    print(f"Wrote {output_path.resolve()} ({atom_count} atoms)")
