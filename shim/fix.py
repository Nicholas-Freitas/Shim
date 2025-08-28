from shim.structure import StandardMolecule, ShimStructure
from shim.utils import get_atom_mapping, validate_structure_file

from tempfile import NamedTemporaryFile
from Bio.PDB import PDBParser, PDBIO

from pathlib import Path

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
    

# ### Test!
# if __name__ == '__main__':
#     fix_atom_names_in_residue('data/cycle_15_top_0_model_0.pdb', 'data/fixed.pdb', 'LIG', 'data/KEMP1_TSA_h.sdf')
#     #rename_lig_chain_by_proximity('data/cycle_15_top_0_model_0.pdb', 'data/renamed.cif', 'LIG')
