from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from pathlib import Path

def get_atom_mapping(mol_A: Chem.Mol, mol_B: Chem.Mol, match_hydrogens=False, silent=False):
    '''
    Given a single-residue RDkit mol from a mol_A file, and a mol_B RDkit mol for that ligand, find the Maximum Common Substructure (MCS) between the two molecules and align the SDF molecule to the mol_A molecule.
    Then, return a mapping of atom indices between the mol_A molecule and the mol_B molecule. Indexing is 0-based.
    ### TODO:
    I want this to be able to match two mols with explicit hydrogens,
    Or get only heavy atoms if one or both are missing hydrogens.
    I've made a patch, but need to test and clean up.
    
    Returns: atom_mapping[mol_A_atom_idx] = mol_B_atom_idx
    '''

    from rdkit.Chem import rdFMCS
    
    if not match_hydrogens:
        # Remove hydrogens to focus on heavy atoms for MCS
        mol_A_noH = Chem.RemoveHs(mol_A)
        mol_B_noH = Chem.RemoveHs(mol_B)
    else:
        # Sneakily, we're still using hydrogens, but calling it noH still.
        mol_A_noH = mol_A
        mol_B_noH = mol_B

    # Find MCS with relaxed bond order comparison
    mcs = rdFMCS.FindMCS(
        [mol_A_noH, mol_B_noH], 
        bondCompare=rdFMCS.BondCompare.CompareAny, 
        ringMatchesRingOnly=True
    )
    if not mcs.smartsString:
        raise ValueError("No MCS found between mol_A and SDF molecules.")

    # Create a query molecule from the MCS
    mcs_query = Chem.MolFromSmarts(mcs.smartsString)
    if mcs_query is None:
        raise ValueError("Could not generate a query molecule from MCS.")

    # Match the MCS query against both molecules
    mol_A_match = mol_A_noH.GetSubstructMatch(mcs_query)
    mol_B_match = mol_B_noH.GetSubstructMatch(mcs_query)

    if not mol_A_match or not mol_B_match:
        raise ValueError("MCS does not match on one of the molecules after all.")

    # mol_A_match and mol_B_match are tuples of atom indices in mol_A_noH and mol_B_noH respectively
    # We need to map back to the original indices in mol_A and mol_B_mol.
    if not match_hydrogens:
        mol_A_heavy_indices = [i for i, a in enumerate(mol_A.GetAtoms()) if a.GetAtomicNum() > 1]
        mol_B_heavy_indices = [i for i, a in enumerate(mol_B.GetAtoms()) if a.GetAtomicNum() > 1]
    else:
        mol_A_heavy_indices = list(range(mol_A.GetNumAtoms()))
        mol_B_heavy_indices = list(range(mol_B.GetNumAtoms()))

    # mol_A_match[i] is the index in mol_A_noH. The atom at mol_A_heavy_indices[mol_A_match[i]] is the original mol_A atom index.
    # Similarly for mol_B_match.
    if len(mol_A_match) != len(mol_B_match):
        raise ValueError("The length of matched MCS atoms does not align between mol_A and mol_B.")

    atom_mapping = {}
    for p_idx, s_idx in zip(mol_A_match, mol_B_match):
        mol_A_atom_idx = mol_A_heavy_indices[p_idx]
        mol_B_atom_idx = mol_B_heavy_indices[s_idx]
        atom_mapping[mol_A_atom_idx] = mol_B_atom_idx

    # Check if the mapping contains all heavy atoms
    num_mol_A_heavy_atoms = mol_A_noH.GetNumAtoms()
    num_mol_B_heavy_atoms = mol_B_noH.GetNumAtoms()

    if len(atom_mapping) != num_mol_A_heavy_atoms or len(atom_mapping) != num_mol_B_heavy_atoms:
        # Which atoms aren't mapped?
        mol_A_mapped_atoms = set(atom_mapping.keys())
        mol_A_unmapped_atoms = set(range(num_mol_A_heavy_atoms)) - mol_A_mapped_atoms
        mol_B_mapped_atoms = set(atom_mapping.values())
        mol_B_unmapped_atoms = set(range(num_mol_B_heavy_atoms)) - mol_B_mapped_atoms

        message = f'The mol_B provided does not contain all heavy atoms from the residue in the mol_A file. \n\
        Unmapped atoms in mol_A: {mol_A_unmapped_atoms}, Unmapped atoms in mol_B: {mol_B_unmapped_atoms}'
        
        if not silent:
            termol.draw(Chem.MolToSmiles(mol_A), name=mol_A_file, three_d=False)
            termol.draw(Chem.MolToSmiles(mol_B), name=mol_B_file, three_d=False)
        raise ValueError(message)
    
    # Align the SDF molecule to the mol_A molecule using this atom mapping
    #AllChem.AlignMol(sdf_mol, mol_A, atomMap=list(atom_mapping.items()))
    
    return atom_mapping

def write_conect_file(mol: Chem.Mol,
                      residue_name: str,
                      chemical_name: str,
                      output_path: str | Path ) -> None:
    """
    Write a PDB-style CONECT block for an RDKit molecule whose atoms have an
    'atomName' property.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Molecule with properly-set bonds and atom property 'atomName'.
    residue_name : str, default 'KP4'
        Three-letter code to appear in RESIDUE / HET lines.
    output_path : str | Path, default 'kp4_conect.txt'
        Path of the file to create.
    """
    # --- gather basic info ----------------------------------------------------
    atom_count = mol.GetNumAtoms()
    formula = rdMolDescriptors.CalcMolFormula(mol)

    # Build “CONECT …” lines
    conect_lines: list[str] = []
    for atom in mol.GetAtoms():
        name = atom.GetProp("atomName")

        # Neighbour names, as stored on the neighbour atoms
        neighbour_names = [n.GetProp("atomName") for n in atom.GetNeighbors()]

        line = ""
        #0-11: CONECT
        #12-18: Atom name
        #19-20: Number of neighbours
        #21+ - Neighbour names, each 5 characters wide

        line += f"CONECT      "
        line += f"{name:<6}"
        line += f"{len(neighbour_names):>2} "
        for neighbour_name in neighbour_names:
            line += f"{neighbour_name:<5}"
        conect_lines.append(line)

    # --- write file -----------------------------------------------------------
    output_path = Path(output_path)
    with output_path.open("w") as fh:
        fh.write(f"RESIDUE   {residue_name:<3}{atom_count:>7}\n")
        fh.writelines(l + "\n" for l in conect_lines)
        fh.write("END   \n")                       # terminate CONECT block
        fh.write(f"HET    {residue_name:<5}{atom_count:>13}\n")
        fh.write(f"HETNAM     {residue_name} {chemical_name}\n")   # customise as needed
        fh.write(f"FORMUL      {residue_name}    {formula}\n")

    print(f"Wrote {output_path.resolve()} ({atom_count} atoms)")


def write_conect_file(mol: Chem.Mol,
                      residue_name: str,
                      chemical_name: str,
                      output_path: str | Path) -> None:
    """
    Write a PDB-style CONECT block for an RDKit molecule whose atoms have an
    'atomName' property, listing heavy atoms first and hydrogens last.
    """
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



