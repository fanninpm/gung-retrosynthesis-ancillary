#! /home/paddy/anaconda3/envs/mol_transformer/bin/python
"""
This script takes a training/testing/validation set and separates the reagents
and reactants.
"""

import re

from dataclasses import dataclass
from typing import List, Set

import click


@dataclass
class Reaction:
    """Class for representing a reaction in SMILES."""

    reactants: List[str]
    reagents: List[str]
    products: List[str]

    def smiles(self):
        """Returns a SMILES representation of this reaction."""
        joined_rxn = (
            ".".join(i) for i in (self.reactants, self.reagents, self.products)
        )
        return ">".join(joined_rxn)


get_atom_map_nos_regex: re.Pattern = re.compile(r":(\d+)\]")


def split_into_rxn_class(smiles: str) -> Reaction:
    """Splits a SMILES string into a Reaction dataclass."""
    smi_list: List[str, str, str] = smiles.split(">")
    assert len(smi_list) == 3  # just a check
    rxn: List[List[str]] = []
    for i in smi_list:
        if i != "":
            rxn.append(i.split("."))
        else:
            rxn.append([])
    return Reaction(reactants=rxn[0], reagents=rxn[1], products=rxn[2])


def get_atom_map_nos(species: str) -> List[int]:
    """Gets a list of atom mapping numbers from a SMILES species."""
    return [int(i.group(1)) for i in get_atom_map_nos_regex.finditer(species)]


def sep_reactants_and_reagents(reaction: Reaction) -> Reaction:
    """Separates atom-mapped species into reactants and reagents."""
    reactants: List[str] = reaction.reactants + reaction.reagents
    new_reactants: List[str] = []
    new_reagents: List[str] = []
    reactants_nos: List[Set[int]] = [
        set(get_atom_map_nos(r)) for r in reactants
    ]
    products_nos: List[Set[int]] = [
        set(get_atom_map_nos(p)) for p in reaction.products
    ]
    for react, reactant_nos in enumerate(reactants_nos):
        for product_nos in products_nos:
            if reactant_nos.intersection(product_nos):
                new_reactants.append(reactants[react])
                break
        else:
            if reactants[react] != "":
                new_reagents.append(reactants[react])
    return Reaction(
        reactants=new_reactants,
        reagents=new_reagents,
        products=reaction.products,
    )


@click.command()
@click.argument(
    "smiles_file",
    type=click.Path(exists=True, readable=True, resolve_path=True),
)
def main(smiles_file: str):
    """
    This script takes a training/testing/validation set and separates the
    reagents and reactants. Requires an input SMILES file. The SMILES file must
    be atom-mapped in order for reactant/reagent separation to work.
    """

    with open(smiles_file) as sf_file:
        for line in sf_file:
            if line == "":
                continue
            just_smi: str = line.split()[0]
            prelim_rxn = split_into_rxn_class(just_smi)
            click.echo(sep_reactants_and_reagents(prelim_rxn).smiles())

            # ~ for i, j in asdict(prelim_rxn).items():
            # ~ click.echo(f"{i}: {j}")


if __name__ == "__main__":
    main()
