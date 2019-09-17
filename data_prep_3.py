#! /home/paddy/anaconda3/envs/mol_transformer/bin/python
"""
Takes a canonicalized SMILES reaction, with reactants and reagents separated,
and tokenizes it.
"""

import re
from typing import List
import click

from data_prep_1 import (
    Reaction,
    split_into_rxn_class,
    sep_reactants_and_reagents,
)
from data_prep_2 import strip_and_canonicalize_rxn

tokenize_pattern: str = r"(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
tokenize_regex: re.Pattern = re.compile(tokenize_pattern)


def tokenize_reaction(rxn: Reaction) -> str:
    """
    Tokenizes the reaction. Reaction tokenizer from Philippe Schwaller et al.
    """
    reactants: str = ".".join(rxn.reactants) + ">"
    reagents: str = ".".join(rxn.reagents)
    products: str = ">" + ".".join(rxn.products)
    react_tokenized: str = " ".join(
        [token for token in tokenize_regex.findall(reactants)]
    )
    reag_tokenized: str = " " + " . ".join(reagents.split(".")) + " "
    prod_tokenized: str = " ".join(
        [token for token in tokenize_regex.findall(products)]
    )

    return react_tokenized + reag_tokenized + prod_tokenized


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
@click.option(
    "--step2",
    "step_2",
    is_flag=True,
    help="Activate this to also run Step 2 (strip atom-mapping and canonicalize reaction).",
)
def main(smiles_file: str, step_1: bool, step_2: bool):
    """
    This script takes a canonicalized SMILES reaction, with reactants and
    reagents separated, and tokenizes it. Requires an input SMILES file. The
    SMILES file must be atom-mapped in order for reactant/reagent separation to
    work (if run with the --step2 flag).
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
            can_rxn: Reaction = strip_and_canonicalize_rxn(
                rr_sep
            ) if step_2 else rr_sep
            tokenized_rxn: str = tokenize_reaction(can_rxn)
            click.echo(tokenized_rxn)


if __name__ == "__main__":
    main()
