"""Splits reaction files (SMILES format, with reaction type annotation)
into src and tgt files for machine learning."""

import os
import re
from typing import List

import click
from tqdm import tqdm


def smi_tokenizer(smi):
    """
    Tokenize a SMILES molecule or reaction.
    Source: https://github.com/pschwllr/MolecularTransformer
    """
    import re

    pattern = r"(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    assert smi == "".join(tokens)
    return " ".join(tokens)


rxn_type_split_regex = re.compile(
    r"(?P<tgt_raw>.*)>(?P<src_products_raw>\S*) (?P<src_rxn_type_raw>.*)$"
)


@click.command()
@click.argument("reaction_file", type=click.Path(exists=True))
def main(reaction_file):
    """
    Splits reaction files (SMILES format, with reaction type annotation)
    into src and tgt files for machine learning.
    
    """
    src_lines: List[str] = []
    tgt_lines: List[str] = []
    head, tail = os.path.split(reaction_file)
    with open(reaction_file) as f:
        lines = f.readlines()
        for line in tqdm(lines):
            rxn_line = line.strip()
            rxn_match = rxn_type_split_regex.match(rxn_line)
            try:
                tgt_line = smi_tokenizer(rxn_match.group("tgt_raw"))
                src_products = smi_tokenizer(
                    rxn_match.group("src_products_raw")
                )
                src_rxn_type = rxn_match.group("src_rxn_type_raw").replace(
                    " ", "_"
                )
            except AttributeError:  # if the regex match returns a NoneType
                continue
            else:
                tgt_lines.append(tgt_line)
                src_lines.append(f"[{src_rxn_type}] {src_products}")
    src_path = os.path.join(head, f"src-{tail}")
    tgt_path = os.path.join(head, f"tgt-{tail}")
    with open(src_path, "w") as src_file, open(tgt_path, "w") as tgt_file:
        src_file.write("\n".join(src_lines))
        tgt_file.write("\n".join(tgt_lines))


if __name__ == "__main__":
    main()
