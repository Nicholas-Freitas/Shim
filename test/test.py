# from shim.fix import fix_atom_names_in_residue

# fix_atom_names_in_residue('test.cif', 
#                            'test_out.cif',
#                            'LIG3', '../shim/data/KEMP1_TSA_h.sdf')

from shim import ShimStructure, StandardMolecule

structure = ShimStructure('test.cif')

ligand = StandardMolecule(structure_file='../shim/data/KEMP1_TSA_h.sdf')

#print('All hetero residues', structure.get_residues(hetero_only=True))

#all_residues = structure.get_residues()

#print('All residues', all_residues)

#print('Matching residues', structure.match_residues_to_mol(all_residues, ligand))


# Create an initial output:
structure.to_cif('test_out_initial.cif')

# Standardize:

structure.standardize([ligand])

# Write the final output:
structure.to_cif('test_out_final.cif')