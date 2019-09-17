#! /home/paddy/anaconda3/envs/mol_transformer/bin/python
"""
Takes an atom-mapped SMILES, strips it of its atom-mapping, and canonicalizes
it.
"""

from typing import List

import click
from rdkit import Chem
from data_prep_1 import (
    Reaction,
    split_into_rxn_class,
    sep_reactants_and_reagents,
)


def strip_and_canonicalize_species(smi: str) -> str:
    mol = Chem.MolFromSmiles(smi)
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, canonical=True)


def strip_and_canonicalize_rxn(rxn: Reaction) -> Reaction:
    reacts: List[str] = [
        strip_and_canonicalize_species(reactant) for reactant in rxn.reactants
    ]
    reags: List[str] = [
        strip_and_canonicalize_species(reagent) for reagent in rxn.reagents
    ]
    prods: List[str] = [
        strip_and_canonicalize_species(product) for product in rxn.products
    ]

    return Reaction(reactants=reacts, reagents=reags, products=prods)


@click.command()
@click.argument(
    "smiles_file",
    type=click.Path(exists=True, readable=True, resolve_path=True),
)
@click.option(
    "--step1",
    "step_1",
    is_flag=True,
    help="Activate this to also run Step 1 (separate reactants and reagents).",
)
def main(smiles_file: str, step_1: bool):
    """
    This script takes an atom-mapped SMILES, strips it of its atom-mapping, and
    canonicalizes it. Requires an input SMILES file. The SMILES file must be
    atom-mapped in order for reactant/reagent separation to work.
    """

    with open(smiles_file) as sf_file:
        for line in sf_file:
            if line == "":
                continue
            if step_1:
                just_smi: str = line.split()[0]
                prelim_rxn: Reaction = split_into_rxn_class(just_smi)
                rr_sep: Reaction = sep_reactants_and_reagents(prelim_rxn)
            else:
                rr_sep: Reaction = split_into_rxn_class(line.split()[0])
            can_rxn = strip_and_canonicalize_rxn(rr_sep)
            click.echo(can_rxn.smiles())


if __name__ == "__main__":
    main()
