from typing import List, Tuple
from tempfile import NamedTemporaryFile

class StandardMolecule:
    def __init__(self, smiles: str = None, 
                        structure_file: str = None,
                        mol = None,
                        use_existing_atom_names=False):
        """
        Initializes a StandardMolecule object.
        TODO: use_existing_atom_names is not implemented
        TODO: During standardization, the renumbering is not tested to be consistent!
        TODO: Not thoroughly tested.
        """

        if use_existing_atom_names:
            raise NotImplementedError("use_existing_atom_names is not implemented yet.")
        
        self.raw_smiles = smiles

        
        self.pdb_block = None
        self.cif_block = None
        self.sdf_block = None
        # Set the pdb, cif, and sdf blocks if provided:
        if structure_file is not None:
            self.read(structure_file)
        self.raw_mol = mol
        
        # We need at least 1 of the above:
        if not (self.raw_smiles or self.pdb_block or self.cif_block or self.sdf_block or self.raw_mol):
            raise ValueError("At least one of smiles, pdb_file, cif_file, sdf_file, or mol must be provided.")

        # Create a standardized representation of the molecule.
        # If the user wants, we'll keep the existing atom names from a CIF or PDB.
        # Otherwise, we'll create standard atom names.
        self.create_standard_mol(use_existing_atom_names=use_existing_atom_names)
        #self.atom_names = self.get_atom_names()

    def read(self, file_path):
        """
        Reads a file and returns its content.
        """

        file_path = str(file_path)
        
        # is it a pdb, cif, or sdf
        if file_path.lower().endswith('.pdb'):
            with open(file_path, 'r') as file:
                self.pdb_block = file.read()
        elif file_path.lower().endswith('.cif'):
            with open(file_path, 'r') as file:
                self.cif_block = file.read()
        elif file_path.lower().endswith('.sdf') or file_path.lower().endswith('.mol'):
            with open(file_path, 'r') as file:
                self.sdf_block = file.read()

    def create_standard_mol(self, use_existing_atom_names=False):
        """
        Creates a standardized representation of the molecule.
        If use_existing_atom_names is True, uses existing atom names from PDB or CIF.
        Otherwise, generates standard atom names.

        We'll only use one possible input, in this order of preference: SMILES, CIF, PDB, SDF, then Mol.
        If multiple inputs are provided, the first non-None preference will be used.
        """
        from rdkit import Chem

        def mol_to_std(raw_mol):
            """
            Converts a raw RDKit molecule to a standardized molecule.
            """
            if not raw_mol:
                raise ValueError("No valid molecule provided.")
            
            # Standardize the molecule by sanitizing it
            Chem.SanitizeMol(raw_mol)
            # Kekulize the molecule to ensure consistent bond orders
            Chem.Kekulize(raw_mol, clearAromaticFlags=True)
            # Convert to canonical smiles for unified format:
            std_smiles = Chem.MolToSmiles(raw_mol, canonical=True)
            # Convert back to a molecule object
            std_mol = Chem.MolFromSmiles(std_smiles, sanitize=False) # for some reason, sanitize=True here removes Hs.
            Chem.SanitizeMol(std_mol)                                # But doesn't here.

            return std_smiles, std_mol

        if self.raw_smiles:
            mol = Chem.MolFromSmiles(self.raw_smiles)
        elif self.pdb_block:
            mol = Chem.MolFromPDBBlock(self.pdb_block, sanitize=False)
        elif self.cif_block:
            mol = Chem.MolFromCIFBlock(self.cif_block, sanitize=False)
        elif self.sdf_block:
            mol = Chem.MolFromMolBlock(self.sdf_block, sanitize=False)
        elif self.raw_mol:
            mol = self.raw_mol
        else:
            raise ValueError("No valid molecule data provided.")

        # Save and set the standard atom names in the atomName value:
        self.std_smiles, self.std_mol = mol_to_std(mol)
        self.set_standard_names()

    def rename_to_standard(self, mol, use_hydrogens=False):
        """
        Renames the atoms in the provided molecule to match the standard names.
        :param mol: RDKit molecule object to rename.
        :return: RDKit molecule with renamed atoms.
        """
        if not hasattr(self, 'std_mol'):
            raise ValueError("Standard molecule has not been created yet.")
        
        from shim.utils import get_atom_mapping

        atom_names = [atom.GetProp('atomName') for atom in self.std_mol.GetAtoms()]
        atom_mapping = get_atom_mapping(mol, self.std_mol, match_hydrogens=use_hydrogens, silent=False)
        
        matched_atoms = []
        unmatched_atoms = []
        for i, atom in enumerate(mol.GetAtoms()):
            if i in atom_mapping:
                # Get the new name from the standard names list
                new_name = atom_names[atom_mapping[i]]
                # Set the atom name property
                atom.SetProp("atomName", new_name)
                matched_atoms.append(f"Atom {i},{atom.GetSymbol()} -> {new_name}")
            else:
                unmatched_atoms.append(f"Atom {i},{atom.GetSymbol()} -> No match found in standard names.")

        # Ensure all atoms have a name property
        if unmatched_atoms:
            raise ValueError(f"Some atoms could not be matched to standard names:\n" + "\n".join(matched_atoms) + "\n"+ "\n".join(unmatched_atoms))
            
        return mol

    def get_mol(self):
        """
        Returns the standardized RDKit molecule.
        """
        if not hasattr(self, 'std_mol'):
            raise ValueError("Standard molecule has not been created yet.")
        return self.std_mol

    def set_standard_names(self):
        """
        Generates standard names for the atoms in the molecule based on their order.
        The first atom will be C1, the second C2, and so on. Oxygen atoms will be O1, O2, etc.
        :param mol: RDKit molecule object.
        :return: None, sets the 'atomName' property for each atom.
        """
        mol = self.std_mol
        if not mol:
            raise ValueError("No molecule provided to generate standard names.")

        standard_names = []
        atom_count = {}  # Extend this dictionary for other elements as needed

        for atom in mol.GetAtoms():
            element = atom.GetSymbol()
            if element in atom_count:
                atom_count[element] += 1
            else:
                atom_count[element] = 1
            name = f"{element}{atom_count[element]}"
            atom.SetProp('atomName', name)

        self.std_mol = mol  # Update the standardized molecule with new names

    def get_smiles(self):
        """
        Returns the standardized SMILES representation of the molecule.
        """
        if not hasattr(self, 'std_smiles'):
            raise ValueError("Standard SMILES has not been created yet.")
        return self.std_smiles
    
    def __repr__(self):
        return f"StandardMolecule(smiles={self.smiles})"

    def show(self):
        """
        Displays in 2D using Termol package.
        """
        pass

class ShimStructure:
    def __init__(self, structure_file=None, structure=None, force=False):
        """
        Thin wrapper around Bio.PDB.Structure.Structure object.
        Created by providing a PDB or CIF file.
        """

        if structure is not None:
            self.structure = structure
        elif structure_file is not None:
            self.structure = self.read_structure(structure_file, force=force)

        # Verify we have a structure:
        if not self.structure:
            raise ValueError("Failed to read structure from provided files.")

    def read_structure(self, structure_file, sanity_check=True, force=False):
        """
        Reads a PDB or CIF file and returns a Biopython structure object.
        """
        from Bio import PDB

        suffix = str(structure_file).split('.')[-1].lower()

        # This allows us to load PDB files with more than 10,000 residues by converting hybrid-36 to decimal.
        # This is an abuse of the PDB format. Ideally, we should use CIF files for such large structures.
        def hybrid36_to_dec(s):
            """ Convert a hybrid-36 encoded string to a decimal integer. """
            # Define the character set for hybrid-36 encoding
            charset = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            base = len(charset)
            
            # Initialize the result
            result = 0
            
            # Iterate over each character in the string
            for i, char in enumerate(reversed(s)):
                if char not in charset:
                    raise ValueError(f"Invalid character '{char}' in hybrid-36 string.")
                value = charset.index(char)
                result += value * (base ** i)
            
            return result

        if suffix == 'pdb':
            parser = PDB.PDBParser(QUIET=False)
            try:
                # Are there more than 10,000 residues?
                too_many_residues = False
                with open(structure_file, 'r') as f:
                    for i, line in enumerate(f):
                        if line.startswith('ATOM') or line.startswith('HETATM'):
                            if hybrid36_to_dec(line[22:26].strip()) > 10000:
                                too_many_residues = True
                                break

                if too_many_residues:
                    if not force:
                        raise ValueError("The PDB file appears to have more than 10,000 residues, which may cause issues with atom naming. Consider using a CIF file instead, or set force=True to proceed anyway.")
                    else:
                        print("Warning: Proceeding with a PDB file that has more than 10,000 residues. This may cause issues with atom naming.")
                        temp_structure = NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False)

                        with open(structure_file, 'r') as original_file:
                            # for each line, replace the hybrid-36 resid with decimal % 10000:
                            for line in original_file:
                                if line.startswith('ATOM') or line.startswith('HETATM'):
                                    resid = line[22:26].strip()
                                    if resid:
                                        dec_resid = hybrid36_to_dec(resid) % 10000
                                        new_line = line[:22] + str(dec_resid).rjust(4) + line[26:]
                                        temp_structure.write(new_line)
                                    else:
                                        temp_structure.write(line)
                                else:
                                    temp_structure.write(line)
                        temp_structure.flush()
                        structure = parser.get_structure('PDB', temp_structure.name)
                else:
                    structure = parser.get_structure('PDB', structure_file)
            
            except Exception as e:
                raise ValueError(f"Failed to parse PDB file: {e}\n. \
                Do you have an improperly formatted PDB file, or more than 10,000 residues? A CIF file may be better.")
        elif suffix == 'cif':
            parser = PDB.MMCIFParser(QUIET=False)
            try:
                structure = parser.get_structure('CIF', structure_file)
            except Exception as e:
                raise ValueError(f"Failed to parse CIF file: {e}")
        else:
            raise ValueError("Either pdb_file or cif_file must be provided.")

        ### I need to figure out how to check when a PDB file has redundant atom names.
        ### Currently, PDBParser just skips duplicate atoms - so we can't check this here.    
        ### We could use a custom StructureBuilder to record duplicates, as drafted above.    
        # if sanity_check:
        #     # Check that each residue doesn't have redundant atom names:
        #     for model in structure:
        #         for chain in model:
        #             for residue in chain:
        #                 atom_names = set()
        #                 print()
        #                 for atom in residue:
        #                     print(atom.get_name(), end=' ')
        #                     if atom.get_name() in atom_names:
        #                         raise ValueError(f"Duplicate atom name {atom.get_name()} found in residue {residue.get_resname()} in chain {chain.id}.")
        #                     atom_names.add(atom.get_name())
        #                 print(atom_names)

        return structure

    ### STANDARDIZE STRUCTURE COMPONENTS ###
    def standardize(self, standard_molecules: List[StandardMolecule], use_hydrogens=False, verbose=False, renumber=False):
        """
        Standardizes the atom names in the structure using a list of StandardMolecule objects.
        :param standard_molecules: List of StandardMolecule objects to use for standardization.
        :return: None, modifies the structure in place.
        """
        # Get all residues in the structure:
        all_residues = self.get_residues()

        # For each standard molecule, find matching residues and standardize them:
        for std_mol in standard_molecules:
            matched_residues = self.match_residues_to_mol(all_residues, std_mol, use_hydrogens=use_hydrogens, verbose=verbose)
            for chain_id, residue_index, resname in matched_residues:
                print(f"Standardizing residue {residue_index} in chain {chain_id} ({resname}) using provided standard molecule.")
                self.standardize_residue(chain_id, residue_index, std_mol, use_hydrogens=use_hydrogens)
        
        # Renumber residues manually:
        if renumber:
            for chain in self.structure.get_chains():
                new_id = 1
                for residue in chain.get_residues():
                    hetflag, old_resseq, icode = residue.id
                    # overwrite residue.id with a new contiguous seq number
                    residue.id = (hetflag, new_id, icode)
                    new_id += 1

    def standardize_residue(self, chain_id, residue_index, standard_molecule: StandardMolecule, use_hydrogens=False):
        """
        Standardizes the atom names in a specific residue of the structure using a StandardMolecule.
        :param chain_id: The chain identifier (e.g. 'A').
        :param residue_index: The residue index (an integer).
        :param standard_molecule: A StandardMolecule object to use for standardization.
        :param use_hydrogens: Whether to consider hydrogens in the matching process.
        :return: None, modifies the structure in place.
        """

        # 1. Extract the residue from the structure
        extracted_res_structure = self.extract_residue(chain_id, residue_index)

        # 2. Convert residue Structure to a Mol object
        extracted_res_mol = extracted_res_structure.to_mol()

        # 3. Rename atoms in the pdb-extracted molecule
        extracted_res_mol = standard_molecule.rename_to_standard(extracted_res_mol, use_hydrogens=use_hydrogens)

        # 4. Replace the atom names in the original Structure with the renamed residue:
        self.rename_residue_with_mol(chain_id, residue_index, extracted_res_mol)

    def match_residues_to_mol(self, residues: List[Tuple[str, int, str]], mol: StandardMolecule, use_hydrogens=False, verbose=False) -> List[Tuple[str, int, str]]:
        """
        Find which residues in the structure match the provided StandardMolecule structure.
        Here we're using the get_atom_mapping function to see if the atoms can be mapped.
        :param residues: List of tuples (chain_id, residue_index) to standardize.
        :param mol: A StandardMolecule object to use for standardization.
        :param use_hydrogens: Whether to consider hydrogens in the matching process.
        :return: List of tuples (chain_id, residue_index) that match the molecule.
        """
        matched_residues = []
        for chain_id, residue_index, resname in residues:
            # 1. Extract the residue from the structure
            extracted_res_structure = self.extract_residue(chain_id, residue_index)

            # 2. Convert residue Structure to a Mol object
            extracted_res_mol = extracted_res_structure.to_mol()

            # 3. Try to rename atoms in the pdb-extracted molecule
            try:
                from shim.utils import get_atom_mapping
                atom_mapping = get_atom_mapping(extracted_res_mol, mol.get_mol(), match_hydrogens=use_hydrogens, silent=True)
                matched_residues.append((chain_id, residue_index, resname))
                if verbose:
                    print(f"SUCCESS! Residue {residue_index} in chain {chain_id} matches the provided molecule with mapping: {atom_mapping}")
            except ValueError as e:
                if verbose:
                    print(f"Residue {residue_index} in chain {chain_id} did not match the provided molecule: {e}")

        return matched_residues

    ### EXTRACT SUBSTRUCTURES ###
    def extract_residue(self, chain_id, residue_index):
        """
        Extracts a specific residue from the structure, and returns a new ShimStructure object
        containing only that residue.
        :param chain_id: The chain identifier (e.g. 'A').
        :param residue_index: The residue index (an integer).
        :return: A new ShimStructure object containing the specified residue.
        """
        
        # 1. Get residue from the Structure object:
        target_residue = None
        for model in self.structure:
            for chain in model:
                if chain.id == chain_id:
                    for res in chain:
                        # Residue ID looks like (" ", residue_number, " ") for standard residues
                        if res.get_id()[1] == residue_index:
                            target_residue = res
                            break
        if target_residue is None:
            raise ValueError(f"Residue {residue_index} not found in chain {chain_id}.")

        # 2. Create a new Structure object with just this residue:
        from Bio.PDB import Structure, Model, Chain, Residue, Atom
        new_structure = Structure.Structure('ExtractedResidue')
        model = Model.Model(0)  # Create a new model
        chain = Chain.Chain(chain_id)  # Create a new chain
        residue = Residue.Residue(target_residue.get_id(), target_residue.get_resname(), target_residue.segid)
        for atom in target_residue:
            # Create a new atom with the same properties
            new_atom = Atom.Atom(atom.get_id(), atom.coord, atom.bfactor, atom.occupancy, atom.altloc,
                                 atom.fullname, atom.serial_number, atom.element)
            residue.add(atom=new_atom)
        chain.add(residue)
        model.add(chain)
        new_structure.add(model)

        # 3. Return a new ShimStructure object with the new structure:
        return ShimStructure(structure=new_structure)

    def extract_chain(self, chain_id):
        """
        Extracts a specific chain from the structure, and returns a new ShimStructure object
        containing only that chain.
        :param chain_id: The chain identifier (e.g. 'A').
        :return: A new ShimStructure object containing the specified chain.
        """
        
        # 1. Get chain from the Structure object:
        target_chain = None
        for model in self.structure:
            for chain in model:
                if chain.id == chain_id:
                    target_chain = chain
                    break
        if target_chain is None:
            raise ValueError(f"Chain {chain_id} not found in structure.")

        # 2. Create a new Structure object with just this chain:
        from Bio.PDB import Structure, Model, Chain, Residue, Atom
        new_structure = Structure.Structure('ExtractedChain')
        model = Model.Model(0)
        chain = Chain.Chain(chain_id)  # Create a new chain
        for residue in target_chain:
            new_residue = Residue.Residue(residue.get_id(), residue.get_resname(), residue.segid)
            for atom in residue:
                # Create a new atom with the same properties
                new_atom = Atom.Atom(atom.get_id(), atom.coord, atom.bfactor, atom.occupancy, atom.altloc,
                                     atom.fullname, atom.serial_number, atom.element)
                new_residue.add(atom=new_atom)
            chain.add(new_residue)
        model.add(chain)
        new_structure.add(model)

        # 3. Return a new ShimStructure object with the new structure:
        return ShimStructure(structure=new_structure)

    ### RENAME STRUCTURE COMPONENTS ###
    def rename_residue_with_mol(self, chain_id, residue_index, mol):
        """
        Renames the atoms in the structure's residue using the provided RDKit molecule.
        :param
        mol: RDKit molecule object with standardized atom names in the 'atomName' field.
        :return: None, modifies the structure in place.

        # TODO: Currently this modifies the structure in place. Should we instead output a copy?
        """

        # 1. Get residue from the Structure object:
        target_residue = None
        for model in self.structure:
            for chain in model:
                if chain.id == chain_id:
                    for res in chain:
                        # Residue ID looks like (" ", residue_number, " ") for standard residues
                        if res.get_id()[1] == residue_index:
                            target_residue = res
                            break
        if target_residue is None:
            raise ValueError(f"Residue {residue_index} not found in chain {chain_id}.")

        # Now, we replace only the names in the atoms of the target residues.
        # To do this, we iteratate over the atoms in both the target residue and the fixed ligand structure - they're ordered the same, so we just swap names.
        if len(target_residue) != len(mol.GetAtoms()):
            raise ValueError(f"Residue {residue_number} in chain {chain_id} of {input_pdb} has {len(target_residue)} atoms, but the ligand PDB has {len(mol.GetAtoms())} atoms.")
        
        for target_atom, ligand_atom in zip(target_residue.get_atoms(), mol.GetAtoms()):
            # Rename the target atom to the ligand atom name
            ligand_atom_name = ligand_atom.GetProp('atomName')

            target_atom.name = ligand_atom_name
            target_atom.fullname = ligand_atom_name#.rjust(4)
            target_atom.id = ligand_atom_name

    def rename_chain(self, old_chain_id, new_chain_id):
        """
        Renames a chain in the structure.
        :param old_chain_id: The current chain identifier (e.g. 'A').
        :param new_chain_id: The new chain identifier (e.g. 'B').
        :return: None, modifies the structure in place.
        """
        for model in self.structure:
            for chain in model:
                if chain.id == old_chain_id:
                    chain.id = new_chain_id
                    return
        raise ValueError(f"Chain {old_chain_id} not found in structure.")

    def rename_residue(self, old_resname, new_resname):
        """
        Renames residues in the structure.
        :param
        old_resname: The current residue name (e.g. 'LIG').
        :param new_resname: The new residue name (e.g. 'LIG1').
        :return: None, modifies the structure in place.
        """       
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() == old_resname:
                        residue.resname = new_resname
        return

    ### OUTPUT STRUCTURE FORMATS ###
    def to_pdb(self, outfile):
        """
        Writes the structure to a PDB file.
        :param outfile: Path to the output PDB file.
        """
        from Bio.PDB import PDBIO

        # truncate atom.id to 3 chars
        for atom in self.structure.get_atoms():
            if len(atom.get_id()) > 3:

                atom.fullname = atom.get_id()[:3]  # Update fullname as well

        # And, truncate the resname to 3 chars
        for residue in self.structure.get_residues():
            if len(residue.get_resname()) > 3:
                residue.resname = residue.get_resname()[:3]  # Update resname to 3 chars

        io = PDBIO()
        io.set_structure(self.structure)
        io.save(outfile, preserve_atom_numbering=True)

    def to_cif(self, outfile):
        """
        Writes the structure to a CIF file.
        :param outfile: Path to the output CIF file.
        """
        from Bio.PDB import MMCIFIO

        io = MMCIFIO()
        io.set_structure(self.structure)
        io.save(outfile, preserve_atom_numbering=True)

    def to_fasta(self, outfile):
        """
        Writes the protein sequences in the structure to a FASTA file.
        :param outfile: Path to the output FASTA file.
        """
        from Bio import SeqIO
        from Bio.PDB import Polypeptide as PP
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        sequences = []
        for model in self.structure:
            for chain in model:
                seq = ''
                for res in chain:
                    if PP.is_aa(res, standard=True):
                        seq += self._resname_to_one(res.get_resname())
                if seq:
                    record = SeqRecord(Seq(seq), id=f"{chain.id}", description="")
                    sequences.append(record)

        SeqIO.write(sequences, outfile, "fasta")

    def to_cif_block(self):
        """
        Returns the CIF block as a string.
        """
        from Bio.PDB import MMCIFIO
        from io import StringIO

        output = StringIO()
        io = MMCIFIO()
        io.set_structure(self.structure)
        io.save(output, preserve_atom_numbering=True)
        return output.getvalue()

    def to_pdb_block(self, suppress_warnings=False):
        """
        Returns the PDB block as a string.
        """
        from Bio.PDB import PDBIO
        from io import StringIO

        # Create a copy of the structure to avoid modifying the original:
        structure_copy = self.structure.copy()

        # Make sure all residue names and atom names are truncated to 3 characters:
        for atom in structure_copy.get_atoms():
            if len(atom.get_id()) > 3:
                atom.fullname = atom.get_id()[:3]
                if not suppress_warnings:
                    print(f"Warning: Atom name {atom.get_id()} truncated to {atom.fullname}.")
        for residue in structure_copy.get_residues():
            if len(residue.get_resname()) > 3:
                residue.resname = residue.get_resname()[:3]
                if not suppress_warnings:
                    print(f"Warning: Residue name {residue.get_resname()} truncated to {residue.resname}.")

        output = StringIO()
        io = PDBIO()
        io.set_structure(structure_copy)
        io.save(output, preserve_atom_numbering=True)
        return output.getvalue()

    def to_fasta_block(self):
        """
        Returns the protein sequences in the structure as a FASTA formatted string.
        """
        from Bio import SeqIO
        from Bio.PDB import Polypeptide as PP
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from io import StringIO

        sequences = []
        for model in self.structure:
            for chain in model:
                seq = ''
                for res in chain:
                    if PP.is_aa(res, standard=True):
                        seq += self._resname_to_one(res.get_resname())
                if seq:
                    record = SeqRecord(Seq(seq), id=f"{chain.id}", description="")
                    sequences.append(record)

        output = StringIO()
        SeqIO.write(sequences, output, "fasta")
        return output.getvalue()
    
    def to_mol(self):
        """
        Converts the structure to an RDKit molecule. (Single residue only)
        :return: An RDKit molecule object.
        """
        from rdkit import Chem

        # Get the number of residues in the structure:
        num_residues = sum(1 for _ in self.structure.get_residues())
        if num_residues > 1:
            raise ValueError("Structure contains multiple residues. Please extract a single residue before converting to RDKit molecule.")

        # # Create a copy of the current ShimStructure as a new ShimStructure object:
        # structure_copy = ShimStructure(structure=self.structure.copy())
        # # Rename the single residue to 'LIG' (or any other name you prefer):
        # structure_copy.rename_residue(old_resname=structure_copy.structure.get_residues().__next__().get_resname(), new_resname='LIG')

        pdb_block = self.to_pdb_block(suppress_warnings=True) ###!
        mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
                        
        if mol is None:
            raise ValueError("Failed to convert structure to RDKit molecule.")
        
        ### Set atomName properties using our saved Structure (the PDB truncates this!):
        bp_atoms = list(self.structure.get_atoms())                 # original order
        for rd_atom, bp_atom in zip(mol.GetAtoms(), bp_atoms):
            rd_atom.SetProp("orig_atom_name", bp_atom.get_name())  # no width limit

        return mol

    ### SET STRUCTURE INFO ###
    def set_b_factor(self, chain_id, residue_index, b_factor):
        """
        Sets the B-factor for all atoms in a specific residue.
        :param chain_id: The chain identifier (e.g. 'A').
        :param residue_index: The residue index (integer).
        :param b_factor: The B-factor value to set (float).
        :return: None, modifies the structure in place.
        """
        for model in self.structure:
            for chain in model:
                if chain.id == chain_id:
                    for res in chain:
                        if res.get_id()[1] == residue_index:
                            for atom in res:
                                atom.set_bfactor(b_factor)
                            return
        raise ValueError(f"Residue {residue_index} not found in chain {chain_id}.")

    ### GET STRUCTURE INFO ###
    def get_residues_by_name(self, resname):
        """
        Returns a list of tuples (chain_id, residue_index) for residues with the given name.
        :param resname: The name of the residue to search for (e.g., 'LIG').
        :return: List of tuples (chain_id, residue_index).
        """
        residues = []
        for model in self.structure:
            for chain in model:
                for res in chain:
                    if res.get_resname() == resname:
                        residues.append((chain.id, res.get_id()[1]))
        return residues
    
    def get_name_by_residue(self, chain_id, residue_index):
        """
        Returns the residue name for a given chain and residue index.
        :param chain_id: The chain identifier (e.g. 'A').
        :param residue_index: The residue index (an integer).
        :return: The residue name (e.g. 'LIG').
        """
        for model in self.structure:
            for chain in model:
                if chain.id == chain_id:
                    for res in chain:
                        if res.get_id()[1] == residue_index:
                            return res.get_resname()
        raise ValueError(f"Residue {residue_index} not found in chain {chain_id}.")

    def get_residues(self, protein_only=False, hetero_only=False):
        """
        Returns a list of tuples (chain_id, residue_index, resname) for all residues.
        :param protein_only: If True, only return protein residues.
        :param hetero_only: If True, only return hetero residues (non-protein).
        :return: List of tuples (chain_id, residue_index, resname).
        """
        from Bio.PDB import Polypeptide as PP

        residues = []
        for model in self.structure:
            for chain in model:
                for res in chain:
                    is_protein = PP.is_aa(res)
                    if (protein_only and not is_protein) or (hetero_only and is_protein):
                        continue
                    residues.append((chain.id, res.get_id()[1], res.get_resname()))
        return residues

    def get_chains(self, protein_only=False, hetero_only=False):
        """
        Return a list of chain IDs in the structure.
        :param protein_only: If True, only return chains that contain protein residues.
        :param hetero_only: If True, only return chains that contain hetero residues (non-protein).
        """

        from Bio.PDB import Polypeptide as PP

        chain_list = []

        for model in self.structure:
            for chain in model:
                has_protein = any(PP.is_aa(res) for res in chain)
                has_hetero = any(not PP.is_aa(res) for res in chain)
                
                if (protein_only and not has_protein) or (hetero_only and not has_hetero):
                    continue
                
                chain_list.append(chain.id)
        return chain_list

    ### UTILITIES ###
    def get_COG(self):
        """ 
        Gets the Center of Geometry of the structure.
        """
        import numpy as np

        coords = []
        for atom in self.structure.get_atoms():
            coords.append(atom.coord)
        coords = np.array(coords)
        cog = np.mean(coords, axis=0)
        return cog

    def get_minimum_distance(self, position: tuple):
        """
        Gets the minimum distance from the given position to any atom in the structure.
        :param position: A tuple of (x, y, z) coordinates.
        :return: The minimum distance.
        """
        import numpy as np

        min_dist = float('inf')
        for atom in self.structure.get_atoms():
            dist = np.linalg.norm(atom.coord - np.array(position))
            if dist < min_dist:
                min_dist = dist
        return min_dist

    def _resname_to_one(self, resname):
        """
        Convert a three-letter residue name to a one-letter code with fallbacks.
        Returns 'X' for unknown residue names.
        """
        # Normalize input
        r = str(resname).strip()
        # Try Biopython Polypeptide converter first
        try:
            from Bio.PDB import Polypeptide as PP
            return PP.three_to_one(r)
        except Exception:
            pass

        # Fallback to Bio.SeqUtils.seq1 which can handle some inputs
        try:
            from Bio.SeqUtils import seq1
            return seq1(r)
        except Exception:
            pass

        # Final fallback: manual mapping for standard amino acids
        mapping = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
            'SEC': 'U',    # selenocysteine
            'PYL': 'O',    # pyrrolysine
            'ASX': 'B',    # Asp or Asn
            'GLX': 'Z',    # Glu or Gln
            'XLE': 'J'     # Leu or Ile (nonstandard code sometimes used)
        }
        key = r[:3].upper()
        return mapping.get(key, 'X')