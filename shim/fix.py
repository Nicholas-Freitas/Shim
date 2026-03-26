from shim.structure import StandardMolecule, ShimStructure
from shim.utils import get_atom_mapping, validate_structure_file

from tempfile import NamedTemporaryFile
from Bio.PDB import PDBParser, PDBIO

from pathlib import Path
import numpy as np

def fix_atom_names_in_residue(infile: str, outfile: str, resname, sdf_file, use_hydrogens=False):
    """
    Fixes the atom names in a specific residue of a PDB file using a standard SDF file.
    :param infile: Path to the input PDB/CIF file.
    :param outfile: Path to the output PDB file with fixed atom names.
    :param resname: The name of the residue to fix (e.g., 'LIG').
    :param sdf_file: Path to the SDF file containing the standard ligand.
    """

    # Is the input a pdb or CIF?
    infile_path, infile_type = validate_structure_file(infile)

    # Is the output file a pdb or CIF?
    outfile_path, outfile_type = validate_structure_file(outfile)

    # Get a clean, standardized SDF molecule as our reference:
    standard_shim_mol = StandardMolecule(sdf_file=sdf_file)

    # Extract the residue from the PDB file
    if infile_type == '.pdb':
        in_structure = ShimStructure(pdb_file=infile)
    elif infile_type == '.cif':
        in_structure = ShimStructure(cif_file=infile)

    # Fix each residue of the specified name:
    for chain, index in in_structure.get_residues_by_name(resname):
        # Extract the residue from the structure
        extracted_res_structure = in_structure.extract_residue(chain, index)

        # Convert residue Structure to a Mol object
        extracted_res_mol = extracted_res_structure.to_mol()

        # Rename atoms in the pdb-extracted molecule
        extracted_res_mol = standard_shim_mol.rename_to_standard(extracted_res_mol, use_hydrogens=use_hydrogens)

        # Replace the atom names in the original PDB with the renamed residue:
        in_structure.rename_residue_with_mol(chain, index, extracted_res_mol)

    # Write the modified structure back to a CIF file
    if outfile_type == '.pdb':
        # If the input was a PDB file, write the output as PDB
        in_structure.to_pdb(outfile)
    elif outfile_type == '.cif':
        # If the input was a CIF file, write the output as CIF
        in_structure.to_cif(outfile)

def rename_lig_chain_by_proximity(infile: str, outfile: str, resname: str, chain_mapping=None):
    '''
    For a given input file, rename the chains of all ligands with the specified residue name
    to match the closest protein chain.
    This is for multimer proteins with multiple ligands bound. We want to ensure that for example,
    protein A always is bound to ligand X, and protein B is always bound to ligand Y.
    This mapping is specified with a dictionary chain_mapping = {'A': 'X', 'B': 'Y'}
    :param infile: Path to the input PDB/CIF file.
    :param outfile: Path to the output PDB/CIF file.
    :param resname: The name of the residue to fix (e.g., 'LIG').
    :param chain_mapping: A dictionary mapping protein chains to ligand chains.
    '''

    # Is the input/output a pdb or CIF?
    infile_path, infile_type = validate_structure_file(infile)
    outfile_path, outfile_type = validate_structure_file(outfile)

    # Get the structure"
    if infile_type == '.pdb':
        in_structure = ShimStructure(pdb_file=infile)
    elif infile_type == '.cif':
        in_structure = ShimStructure(cif_file=infile)

    # Mapping from protein chain -> ligand chain:
    if chain_mapping is None:
        chain_mapping = {'A': 'X',
                        'B': 'Y',
                        'C': 'Z'}
    # Get all protein chains (list of strings CHAIN)
    protein_chains = in_structure.get_protein_chains()

    # Get all ligand chains with the specified residue name (list of tuples (chain, resnum))
    ligand_chains_resnums = in_structure.get_residues_by_name(resname)

    # Ligand rename map: We'll save the mapping of old ligand chain to new ligand chain here, then rename at the end:
    ligand_rename_map = {}

    # For each ligand chain, find the closest protein chain:
    for ligand_chain, ligand_resnum in ligand_chains_resnums:

        # Ligand Structure
        ligand = in_structure.extract_residue(ligand_chain, ligand_resnum)

        # Get the coordinates of the ligand residue
        ligand_coords = ligand.get_COG()

        closest_chain = None
        closest_distance = float('inf')

        for protein_chain in protein_chains:
            # Get the coordinates of the protein chain
            current_chain_structure = in_structure.extract_chain(protein_chain)

            current_chain_closest_distance = current_chain_structure.get_minimum_distance(ligand_coords)

            if current_chain_closest_distance < closest_distance:
                closest_distance = current_chain_closest_distance
                closest_chain = protein_chain

        # Use the mapping to determine the new ligand_chain:
        if closest_chain in chain_mapping:
            new_ligand_chain = chain_mapping[closest_chain]
            ligand_rename_map[ligand_chain] = new_ligand_chain
        else:
            mapping_string = '\n'.join([f' Protein {k} -> Ligand {v}' for k, v in chain_mapping.items()])
            raise ValueError(f'No mapping found for protein chain {closest_chain} in chain_mapping.\n \
                Currently mapping: \n \
                {mapping_string}')

    # Rename the ligand chains in the structure:
    for old_chain, new_chain in ligand_rename_map.items():
        in_structure.rename_chain(old_chain, new_chain)

    # Write the modified structure back to a file
    if outfile_type == '.pdb':
        in_structure.to_pdb(outfile)
    elif outfile_type == '.cif':
        in_structure.to_cif(outfile)

    new_structure = ShimStructure(pdb_file=outfile) if outfile_type == '.pdb' else ShimStructure(cif_file=outfile)
    
    for chain in new_structure.structure[0]:
       print(f'Protein chain: {chain}')
    
def rename_atom_by_proximity(infile: str, outfile: str, 
        target_chain: str, 
        target_resid: int,
        target_atom_name: str,
        renaming_chain: str,
        renaming_resid: int,
        closest_atom_name: str,
        ):
    '''
    For a given input file, identify the target atom and the "renaming" residue. Rename the atoms in 
    the "renaming" residue so the atom identified by "closest_atom_name". 
    This is for structures where we want to get the distance from a target atom in a symmetric residue, like ASP or GLU.
    We can rename them so OD1 is always closest to our ligand.
    
    :param infile: Path to the input PDB/CIF file.
    :param outfile: Path to the output PDB/CIF file.
    :param target_chain: The chain of the target atom.
    :param target_resid: The residue number of the target atom.
    :param target_atom_name: The name of the target atom.
    :param renaming_chain: The chain of the renaming residue.
    :param renaming_resid: The residue number of the renaming residue.
    :param closest_atom_name: The name of the atom to rename to.
    
    '''

    SYMMETRIC_ATOMS = {
        "ASP": [["OD1", "OD2"]],
        "GLU": [["OE1", "OE2"]],
        "PHE": [["CD1", "CD2"], ["CE1", "CE2"]],
        "TYR": [["CD1", "CD2"], ["CE1", "CE2"]],
        "ARG": [["NH1", "NH2"]],
        "LEU": [["CD1", "CD2"]],
        "VAL": [["CG1", "CG2"]],
    }
    # Is the input/output a pdb or CIF?
    infile_path, infile_type = validate_structure_file(infile)
    outfile_path, outfile_type = validate_structure_file(outfile)

    # Get the structure"
    if infile_type == '.pdb':
        in_structure = ShimStructure(structure_file=infile)
    elif infile_type == '.cif':
        in_structure = ShimStructure(structure_file=infile)

    # 1. Get the target atom coordinates
    target_atom = None
    for model in in_structure.structure:
        for chain in model:
            if chain.id == target_chain:
                for residue in chain:
                    # Residue ID looks like (" ", residue_number, " ")
                    if residue.get_id()[1] == target_resid:
                        if target_atom_name in residue:
                            target_atom = residue[target_atom_name]
                        break
    
    if target_atom is None:
        raise ValueError(f"Target atom {target_atom_name} not found in residue {target_resid} of chain {target_chain}.")
    
    target_coord = target_atom.coord

    # 2. Get the renaming residue
    renaming_residue = None
    for model in in_structure.structure:
        for chain in model:
            if chain.id == renaming_chain:
                for residue in chain:
                    if residue.get_id()[1] == renaming_resid:
                        renaming_residue = residue
                        break
    
    if renaming_residue is None:
        raise ValueError(f"Renaming residue {renaming_resid} not found in chain {renaming_chain}.")

    resname = renaming_residue.get_resname()
    if resname not in SYMMETRIC_ATOMS:
        print(f"Warning: Residue {resname} is not in the SYMMETRIC_ATOMS list. No action taken.")
        return

    # 3. Determine if we need to flip
    # Find the pair containing closest_atom_name
    primary_pair = None
    for pair in SYMMETRIC_ATOMS[resname]:
        if closest_atom_name in pair:
            primary_pair = pair
            break
    
    if primary_pair is None:
        raise ValueError(f"closest_atom_name {closest_atom_name} not found in SYMMETRIC_ATOMS for {resname}.")

    atom1_name, atom2_name = primary_pair
    if atom1_name not in renaming_residue or atom2_name not in renaming_residue:
        raise ValueError(f"Atoms {atom1_name} or {atom2_name} not found in residue {resname} {renaming_resid}.")

    dist1 = np.linalg.norm(renaming_residue[atom1_name].coord - target_coord)
    dist2 = np.linalg.norm(renaming_residue[atom2_name].coord - target_coord)

    should_flip = (closest_atom_name == atom1_name and dist2 < dist1) or \
                  (closest_atom_name == atom2_name and dist1 < dist2)

    if should_flip:
        print(f"Flipping symmetric atoms in {resname} {renaming_resid} because {closest_atom_name} is further than its symmetric partner.")
        
        swaps = []
        for p in SYMMETRIC_ATOMS[resname]:
            n1, n2 = p[0], p[1]
            swaps.append((n1, n2))
            
            # Hydrogen heuristic: find the differing character (e.g., '1' vs '2')
            # and swap any hydrogens that have that difference at the same position.
            diff_idx = -1
            for i in range(min(len(n1), len(n2))):
                if n1[i] != n2[i]:
                    diff_idx = i
                    break
            
            if diff_idx != -1:
                c1, c2 = n1[diff_idx], n2[diff_idx]
                for atom in renaming_residue:
                    h_name = atom.get_name()
                    if h_name.startswith('H') and len(h_name) > diff_idx:
                        if h_name[diff_idx] == c1:
                            h2_name = h_name[:diff_idx] + c2 + h_name[diff_idx+1:]
                            if h2_name in renaming_residue and (h_name, h2_name) not in swaps:
                                swaps.append((h_name, h2_name))

        # Perform the swaps
        for n1, n2 in swaps:
            a1, a2 = renaming_residue[n1], renaming_residue[n2]
            # Swap name, fullname, id
            a1.name, a2.name = a2.name, a1.name
            a1.fullname, a2.fullname = a2.fullname, a1.fullname
            a1.id, a2.id = a2.id, a1.id
    else:
        print(f"No flip needed for {resname} {renaming_resid}.")

    # Write the modified structure back
    if outfile_type == '.pdb':
        in_structure.to_pdb(outfile)
    elif outfile_type == '.cif':
        in_structure.to_cif(outfile)
    

# ### Test!
# if __name__ == '__main__':
#     fix_atom_names_in_residue('data/cycle_15_top_0_model_0.pdb', 'data/fixed.pdb', 'LIG', 'data/KEMP1_TSA_h.sdf')
#     #rename_lig_chain_by_proximity('data/cycle_15_top_0_model_0.pdb', 'data/renamed.cif', 'LIG')
