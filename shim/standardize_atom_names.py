from data import StandardMolecule, ShimStructure
from utils import get_atom_mapping

from tempfile import NamedTemporaryFile
from Bio.PDB import PDBParser, PDBIO

def fix_atom_names_in_residue(infile: str, outfile: str, chain, index, sdf_file):
    """
    Fixes the atom names in a specific residue of a PDB file using a standard SDF file.
    :param infile: Path to the input PDB file.
    :param outfile: Path to the output PDB file with fixed atom names.
    :param chain: The chain identifier (e.g. 'A').
    :param index: The residue index (an integer).
    :param sdf_file: Path to the SDF file containing the standard ligand.
    """
    # Get a clean, standardized SDF molecule as our reference:
    standard_shim_mol = StandardMolecule(sdf_file=sdf_file)

    # Extract the residue from the PDB file
    in_structure = ShimStructure(pdb_file=infile)
    extracted_res_structure = in_structure.extract_residue(chain, index)

    #temp_residue_file = NamedTemporaryFile(delete=False, suffix='.pdb').name  # temporary file for the residue
    #extract_residue_from_pdb(infile, chain, index, temp_residue_file)
    
    # Read the extracted residue
    ### Convert residue Structure to a Mol object
    ### ! Important! This assumes the atom ordering is preserved when we convert from PDB block to Rdkit Mol.
    ### It seems like this is true, but it needs testing.
    extracted_res_mol = extracted_res_structure.to_mol()
    
    # Rename atoms in the pdb-extracted molecule
    extracted_res_mol = standard_shim_mol.rename_to_standard(extracted_res_mol)

    # Replace the atom names in the original PDB with the renamed residue:
    in_structure.rename_residue_with_mol( chain, index, extracted_res_mol )
    in_structure.to_pdb(outfile)

    for atom in standard_shim_mol.get_mol().GetAtoms():
        print(f"Atom {atom.GetIdx()} {atom.GetSymbol()} -> {atom.GetProp('atomName')}")

    from utils import write_conect_file
    write_conect_file(standard_shim_mol.get_mol(), residue_name="KP4", chemical_name="KEMp-1", output_path="kp4_conect.txt")


### Test!
fix_atom_names_in_residue('data/cycle_15_top_0_model_0.pdb', 
                           'data/cycle_15_top_0_model_0_fixed.pdb',
                           'C', 1, 'data/KEMP1_TSA_h.sdf')