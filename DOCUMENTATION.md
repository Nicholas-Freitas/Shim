Shim
====

Overview
--------
Shim is a small Python package that helps standardize and fix common molecular structure formats (PDB, CIF, SDF) and provides lightweight conversion helpers to act as a "shim" between different bioinformatics tools. It is designed for workflows where ligand/residue atom names, chain assignments and connectivity need to be normalized before downstream processing.

Quick install
-------------
Use your preferred environment manager. The package requires Python >= 3.10 and depends on RDKit and Biopython.

Example (conda-like):

- Create and activate an env with RDKit and Biopython installed.
- Install the package locally for development: pip install -e .

Core concepts and objects
-------------------------

- `StandardMolecule` (in `shim.structure`)
  - Purpose: represent a single small molecule (ligand) in a canonical, standardized RDKit form. Useful as a reference for atom names and ordering.
  - When to use: when you have a reference SDF/SMILES/PDB/CIF for a ligand and you want consistent atom names or to rename an extracted ligand from a structure to use the standard names.
  - Key behaviours:
    - Read molecule input from SMILES, PDB, CIF, SDF or an RDKit Mol.
    - Create a canonical RDKit `Mol` (`std_mol`) and a standard SMILES (`std_smiles`).
    - Attach per-atom properties containing a standard `atomName` value.
    - Provide `rename_to_standard(mol, use_hydrogens=False)` to rename atoms on a given RDKit `Mol` to the standard names from the reference.
  - Example usage:

```python
from shim.structure import StandardMolecule

# Create a standard from an SDF file
std = StandardMolecule(sdf_file='data/KEMP1_TSA_h.sdf')

# Given another RDKit mol extracted from a PDB, rename atoms to the standard
renamed = std.rename_to_standard(extracted_mol, use_hydrogens=False)
```

- `ShimStructure` (in `shim.structure`)
  - Purpose: thin wrapper around a Biopython `Structure` object to make common tasks easier: extracting residues or chains, renaming atoms/residues/chains, and converting between PDB/CIF and RDKit `Mol` forms.
  - When to use: when you need to parse a PDB/CIF, extract a ligand/residue, rename atoms, or write back to disk in a consistent form.
  - Typical operations:
    - Construct from file: `ShimStructure(pdb_file='path/to/file.pdb')` or `ShimStructure(cif_file='path/to/file.cif')`.
    - `extract_residue(chain_id, residue_index)` to produce a new `ShimStructure` with only that residue.
    - `rename_residue_with_mol(chain_id, residue_index, mol)` to replace the atom names in the structure with names from an RDKit `Mol`.
    - `to_pdb(outfile)` / `to_cif(outfile)` to write the modified structure to disk.
  - Example usage:

```python
from shim.structure import ShimStructure

s = ShimStructure(pdb_file='data/out1.pdb')
res = s.extract_residue('X', 100)
# do operations on `res`, then write back
s.to_pdb('data/fixed.pdb')
```

Utilities and helpers
---------------------

- `shim.utils.get_atom_mapping(mol_A, mol_B, match_hydrogens=False)`
  - Finds a maximum common substructure (MCS) between two RDKit molecules and returns a mapping of atom indices from `mol_A` -> `mol_B`.
  - Useful for aligning and mapping atom names between an extracted residue (from a PDB) and a canonical SDF.

- `shim.utils.validate_structure_file(path)`
  - Quick validator that returns a `Path` and the suffix type (`.pdb` or `.cif`) or raises for unsupported types.

File conversion and helpers
---------------------------

- `shim.convert.write_reduce_db(standard_mol, residue_name, chemical_name, output_path)`
  - Given a `StandardMolecule` (or any object exposing `std_mol`), write a PDB-style CONECT block and a small Reduce/DB-like header (RESIDUE / HET / HETNAM / FORMUL). Useful when preparing residue dictionaries or small-molecule connectivity files for other tools.

Examples
--------

1) Fixing atom names in a PDB residue using an SDF reference (script-like):

```python
from shim.fix import fix_atom_names_in_residue

fix_atom_names_in_residue(
    infile='data/cycle_15_top_0_model_0.pdb',
    outfile='data/fixed.pdb',
    resname='LIG',
    sdf_file='data/KEMP1_TSA_h.sdf',
    use_hydrogens=False
)
```

2) Renaming ligand chains based on proximity to protein chains:

```python
from shim.fix import rename_lig_chain_by_proximity

rename_lig_chain_by_proximity(
    infile='data/cycle_15_top_0_model_0.pdb',
    outfile='data/renamed.cif',
    resname='LIG',
    chain_mapping={'A':'X','B':'Y','C':'Z'}
)
```

3) Producing a small connectivity file for a standard molecule:

```python
from shim.structure import StandardMolecule
from shim.convert import write_reduce_db

std = StandardMolecule(sdf_file='data/KEMP1_TSA_h.sdf')
write_reduce_db(std, residue_name='LIG', chemical_name='KEMP1', output_path='data/LIG.con')
```

Notes, limitations, and development status
-----------------------------------------
- Several functions in `shim.structure` are thin wrappers around Biopython and RDKit. Implementation details (and some TODOs) exist in the source; behaviour follows the intended API but edge-cases (duplicated atom names, alternate locations, missing hydrogens) may need extra handling in complex PDBs.
- `StandardMolecule.create_standard_mol` prefers inputs in the order: SMILES, CIF, PDB, SDF, then RDKit Mol. If multiple inputs are provided only the first non-empty input is used.
- `utils.get_atom_mapping` uses RDKit's MCS and can raise a `ValueError` if the reference SDF does not contain all heavy atoms present in the extracted residue. Visual debugging is aided by `termol` render calls when not silent.

License
-------
This project is MIT licensed — see `LICENSE` for details.
