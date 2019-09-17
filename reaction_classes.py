"""Support classes for the main file."""

from dataclasses import dataclass, field
from typing import List

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFMCS
from rdkit.Chem.Fingerprints import FingerprintMols


@dataclass
class Reaction:
    """Class for representing a reaction in SMILES."""

    reactants: List[str]
    reagents: List[str]
    products: List[str]

    mol_reactants: List[Chem.rdchem.Mol] = field(init=False)
    mol_reagents: List[Chem.rdchem.Mol] = field(init=False)
    mol_products: List[Chem.rdchem.Mol] = field(init=False)

    fp_reactants: List[DataStructs.cDataStructs.ExplicitBitVect] = field(init=False)
    fp_reagents: List[DataStructs.cDataStructs.ExplicitBitVect] = field(init=False)
    fp_products: List[DataStructs.cDataStructs.ExplicitBitVect] = field(init=False)

    def __post_init__(self):
        """Properly initialize some class variable representations (or try)."""

        self.mol_reactants = [Chem.MolFromSmiles(x) for x in self.reactants]
        self.mol_reagents = [Chem.MolFromSmiles(x) for x in self.reagents]
        self.mol_products = [Chem.MolFromSmiles(x) for x in self.products]
        self.fp_reactants = [
            FingerprintMols.FingerprintMol(x) for x in self.mol_reactants
        ]
        self.fp_reagents = [
            FingerprintMols.FingerprintMol(x) for x in self.mol_reagents
        ]
        self.fp_products = [
            FingerprintMols.FingerprintMol(x) for x in self.mol_products
        ]

    def smiles(self):
        """Returns a SMILES representation of this reaction."""

        joined_rxn = (
            ".".join(i) for i in (self.reactants, self.reagents, self.products)
        )
        return ">".join(joined_rxn)


@dataclass
class ReactionDiff:
    """Class for representing a reaction difference."""

    mol_rel_reactant: Chem.rdchem.Mol
    mol_product: Chem.rdchem.Mol
    mol_reactant_dissim: Chem.rdchem.Mol
    mol_product_dissim: Chem.rdchem.Mol
    mol_reactants: Chem.rdchem.Mol

    rel_reactant: str = field(init=False)
    product: str = field(init=False)
    reactant_dissim: str = field(init=False)
    product_dissim: str = field(init=False)
    reactants: str = field(init=False)

    def __post_init__(self):
        """Properly initialize some class variable representations (or try)."""

        self.rel_reactant = Chem.MolToSmiles(self.mol_rel_reactant)
        self.product = Chem.MolToSmiles(self.mol_product)
        self.reactant_dissim = Chem.MolToSmiles(self.mol_reactant_dissim)
        self.product_dissim = Chem.MolToSmiles(self.mol_product_dissim)
        self.reactants = Chem.MolToSmiles(self.mol_reactants)
