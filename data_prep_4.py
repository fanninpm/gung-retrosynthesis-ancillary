#! /home/paddy/anaconda3/envs/mol_transformer/bin/python
"""
Takes a tokenized SMILES reaction and splits it into src and tgt files.
"""
import re

import click

from data_prep_1 import (
    Reaction,
    split_into_rxn_class,
    sep_reactants_and_reagents,
)
from data_prep_2 import strip_and_canonicalize_rxn
from data_prep_3 import tokenize_reaction

split_rxn_regex: re.Pattern = re.compile(r"(?P<src>.*) > ?(?P<tgt>.*)$")


def split_rxn(tok_smi: str, label: str) -> str:
    return split_rxn_regex.match(tok_smi).group(label)


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
@click.option(
    "--step3",
    "step_3",
    is_flag=True,
    help="Activate this to also run Step 3 (tokenize reaction).",
)
@click.option(
    "--label",
    required=True,
    type=click.Choice(["src", "tgt"]),
    help="What type of file to export.",
)
def main(
    smiles_file: str, step_1: bool, step_2: bool, step_3: bool, label: str
):
    """
    This script takes a tokenized SMILES file and splits it into src and tgt
    files. Requires an input SMILES file. The SMILES file must be atom-mapped in
    order for reactant/reagent separation to work (if run with the --step2
    flag).
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
                rr_sep: Reaction = split_into_rxn_class(
                    line.split()[0]
                ) if step_3 else line
            can_rxn: Reaction = strip_and_canonicalize_rxn(
                rr_sep
            ) if step_2 else rr_sep
            tokenized_rxn: str = tokenize_reaction(
                can_rxn
            ) if step_3 else can_rxn
            click.echo(split_rxn(tokenized_rxn, label))


if __name__ == "__main__":
    main()
