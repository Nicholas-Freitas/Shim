from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from pathlib import Path
import termol

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
            termol.draw(Chem.MolToSmiles(mol_A), name='mol_A', three_d=False)
            termol.draw(Chem.MolToSmiles(mol_B), name='mol_B', three_d=False)
        raise ValueError(message)
    
    # Align the SDF molecule to the mol_A molecule using this atom mapping
    #AllChem.AlignMol(sdf_mol, mol_A, atomMap=list(atom_mapping.items()))
    
    return atom_mapping


