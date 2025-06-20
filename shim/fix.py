from shim.structure import StandardMolecule, ShimStructure
from shim.utils import get_atom_mapping

from tempfile import NamedTemporaryFile
from Bio.PDB import PDBParser, PDBIO

from pathlib import Path

def fix_atom_names_in_residue(infile: str, outfile: str, resname, sdf_file):
    """
    Fixes the atom names in a specific residue of a PDB file using a standard SDF file.
    :param infile: Path to the input PDB/CIF file.
    :param outfile: Path to the output PDB file with fixed atom names.
    :param resname: The name of the residue to fix (e.g., 'LIG').
    :param sdf_file: Path to the SDF file containing the standard ligand.
    """

    # Is the input a pdb or CIF?
    infile_path = Path(infile)
    infile_type = infile_path.suffix.lower()
    if infile_type in ['.pdb', '.cif']:
        # If it's a PDB or CIF file, we can proceed.
        pass
    else:
        raise ValueError(f"Unsupported file type: {infile_path.suffix}. Expected .pdb or .cif.")

    # Is the output file a pdb or CIF?
    outfile_path = Path(outfile)
    outfile_type = outfile_path.suffix.lower()
    if outfile_type not in ['.pdb', '.cif']:
        raise ValueError(f"Unsupported output file type: {outfile_path.suffix}. Expected .pdb or .cif.")

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
        extracted_res_mol = standard_shim_mol.rename_to_standard(extracted_res_mol)

        # Replace the atom names in the original PDB with the renamed residue:
        in_structure.rename_residue_with_mol(chain, index, extracted_res_mol)

    # Write the modified structure back to a CIF file
    if outfile_type == '.pdb':
        # If the input was a PDB file, write the output as PDB
        in_structure.to_pdb(outfile)
    elif outfile_type == '.cif':
        # If the input was a CIF file, write the output as CIF
        in_structure.to_cif(outfile)


# fix_atom_names_in_residue('data/TKSI_KEMP1_dimer.cif', 
#                            'data/TKSI_KEMP1_dimer_fixed.cif',
#                            'LIG', 'data/KEMP1_h.sdf')

# # While we're here, let's write a REDUCE file for the standard molecule:
# from utils import write_conect_file
# write_conect_file(standard_shim_mol.get_mol(), residue_name=resname, chemical_name="KEMp-1", output_path=f"data/{resname}_conect.txt")