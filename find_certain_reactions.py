"""Tries to find the differences between products and reactants."""

import operator
from collections import Counter
from itertools import repeat
from typing import List

import click
import tablib
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdFMCS
from tqdm import tqdm

from reaction_classes import Reaction, ReactionDiff


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


def find_reaction_difference(rxn: Reaction) -> List[ReactionDiff]:  # return type?
    """Finds difference between products and reactants."""
    rxn_diff_list: list = []

    # ~ for mol_product in rxn.mol_products:
    # ~ for mol_reactant in rxn.mol_reactants:
    # ~ ms = [mol_reactant, mol_product]
    # ~ mcs = rdFMCS.FindMCS(ms)
    # ~ patt = mcs.smartsString
    # ~ commonMol = Chem.MolFromSmarts(patt)
    # ~ mol_diffs = [Chem.DeleteSubstructs(x,commonMol) for x in ms]
    # ~ diff_list.append([Chem.MolToSmiles(x) for x in mol_diffs])

    for i, fp_product in enumerate(rxn.fp_products):
        tamimoto_scores = [
            DataStructs.FingerprintSimilarity(x, fp_product) for x in rxn.fp_reactants
        ]
        j, max_tamimoto_score = max(
            enumerate(tamimoto_scores), key=operator.itemgetter(1)
        )
        mol_product = rxn.mol_products[i]
        mol_reactant = rxn.mol_reactants[j]
        ms = [mol_reactant, mol_product]
        mcs = rdFMCS.FindMCS(ms, timeout=10)
        if mcs.canceled == False:
            patt = mcs.smartsString
            commonMol = Chem.MolFromSmarts(patt)
            mol_diffs = [Chem.DeleteSubstructs(x, commonMol) for x in ms]
            diffs = ReactionDiff(
                mol_rel_reactant=mol_reactant,
                mol_product=mol_product,
                mol_reactant_dissim=mol_diffs[0],
                mol_product_dissim=mol_diffs[1],
                mol_reactants=Chem.MolFromSmiles(".".join(rxn.reactants)),
            )
            rxn_diff_list.append(diffs)
        else:
            diffs = ReactionDiff(
                mol_rel_reactant=mol_reactant,
                mol_product=mol_product,
                mol_reactant_dissim=Chem.MolFromSmiles(""),
                mol_product_dissim=Chem.MolFromSmiles(""),
                mol_reactants=Chem.MolFromSmiles(".".join(rxn.reactants)),
            )
            rxn_diff_list.append(diffs)
        # ~ diff_list = [
        # ~ rxn.reactants[j],
        # ~ [Chem.MolToSmiles(x) for x in mol_diffs],
        # ~ rxn.products[i],
        # ~ ]
        # ~ rxn_diff_list.append(diff_list)

    return rxn_diff_list


urea_patt = Chem.MolFromSmiles("NC(=O)N")
ester_patt = Chem.MolFromSmiles("CC(=O)OC")
amide_patt = Chem.MolFromSmiles("CC(=O)N")
carbamate_patt = Chem.MolFromSmiles("OC(=O)N")
cooh_patt = Chem.MolFromSmarts("C(=O)[OH]")
weinrab_patt = Chem.MolFromSmiles("CC(=O)N(C)OC")
ketone_patt = Chem.MolFromSmiles("CC(=O)C")
aldehyde_patt = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
alcohol_patt = Chem.MolFromSmarts("[!$([#6](=O)[OH]);$([#6][OH])]")
arom_nit_patt = Chem.MolFromSmarts("c([N+](=O)[O-])")
sulfonamide_patt = Chem.MolFromSmiles("O=S(=O)N")
ether_patt = Chem.MolFromSmarts(
    "[$([OD2]([#6])[#6]);!$([OD2]([#6])[#6]=O);!$([OD2]1[#6][#6]1)]"
)
carbonyl_patt = Chem.MolFromSmiles("C=O")
amine_patt = Chem.MolFromSmarts("[$([#6][NX3,NX4]);!$([#6](=O)N)]")
nitroso_patt = Chem.MolFromSmarts("[#6]N=O")


def classify_reaction(rxn: Reaction, rxn_diffs: ReactionDiff):
    """Currently only classifies an addition of a bromine/chlorine atom to the product."""

    # ~ reactant_br_count = rxn_diff_list[0].count("CBr") + rxn_diff_list[0].count("BrC")
    # ~ product_br_count = rxn_diff_list[2].count("CBr") + rxn_diff_list[2].count("BrC")
    # ~ reactant_cl_count = rxn_diff_list[0].count("CCl") + rxn_diff_list[0].count("ClC")
    # ~ product_cl_count  = rxn_diff_list[2].count("CCl") + rxn_diff_list[2].count("ClC")

    reactant = rxn_diffs.rel_reactant
    mol_reactant = rxn_diffs.mol_rel_reactant
    reactant_list = rxn_diffs.reactants
    mol_reactant_list = rxn_diffs.mol_reactants
    product = rxn_diffs.product
    mol_product = rxn_diffs.mol_product

    # ~ reactant_double_bond_count = (
    # ~ reactant.count("C=C") + reactant.count("=O") + reactant.count("O=")
    # ~ )
    # ~ product_double_bond_count = (
    # ~ product.count("C=C") + product.count("=O") + product.count("O=")
    # ~ )

    reactant_bonds = mol_reactant_list.GetBonds()
    product_bonds = mol_product.GetBonds()
    reactant_bond_types = [bond.GetBondTypeAsDouble() for bond in reactant_bonds]
    product_bond_types = [bond.GetBondTypeAsDouble() for bond in product_bonds]
    reactant_double_bond_count = reactant_bond_types.count(2.0)
    product_double_bond_count = product_bond_types.count(2.0)
    reactant_aromatic_attachments = [
        frozenset(
            (bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum())
        )
        for bond in reactant_bonds
        if bond.GetBondTypeAsDouble() == 1.0
        and operator.xor(
            bond.GetBeginAtom().GetIsAromatic(), bond.GetEndAtom().GetIsAromatic()
        )
    ]
    product_aromatic_attachments = [
        frozenset(
            (bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum())
        )
        for bond in product_bonds
        if bond.GetBondTypeAsDouble() == 1.0
        and operator.xor(
            bond.GetBeginAtom().GetIsAromatic(), bond.GetEndAtom().GetIsAromatic()
        )
    ]
    c_react_arom_attach = Counter(reactant_aromatic_attachments)
    c_prod_arom_attach = Counter(product_aromatic_attachments)
    reactant_arom_arom = [
        frozenset(
            (bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum())
        )
        for bond in reactant_bonds
        if bond.GetBondTypeAsDouble() == 1.0
        and (bond.GetBeginAtom().GetIsAromatic() and bond.GetEndAtom().GetIsAromatic())
    ]
    product_arom_arom = [
        frozenset(
            (bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum())
        )
        for bond in product_bonds
        if bond.GetBondTypeAsDouble() == 1.0
        and (bond.GetBeginAtom().GetIsAromatic() and bond.GetEndAtom().GetIsAromatic())
    ]
    c_react_arom_arom = Counter(reactant_arom_arom)
    c_prod_arom_arom = Counter(product_arom_arom)
    reactant_carbon_carbon = [
        bond
        for bond in reactant_bonds
        if {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()} == {6}
    ]
    product_carbon_carbon = [
        bond
        for bond in product_bonds
        if {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()} == {6}
    ]
    reactant_carbon_nitrogen = [
        bond
        for bond in reactant_bonds
        if {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()}
        == {6, 7}
    ]
    product_carbon_nitrogen = [
        bond
        for bond in product_bonds
        if {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()}
        == {6, 7}
    ]
    reactant_carbon_sulfur = [
        bond
        for bond in reactant_bonds
        if {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()}
        == {6, 16}
    ]
    product_carbon_sulfur = [
        bond
        for bond in product_bonds
        if {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()}
        == {6, 16}
    ]
    reactant_carbon_chlorine = [
        bond
        for bond in reactant_bonds
        if {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()}
        == {6, 17}
    ]
    product_carbon_chlorine = [
        bond
        for bond in product_bonds
        if {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()}
        == {6, 17}
    ]
    reactant_carbon_bromine = [
        bond
        for bond in reactant_bonds
        if {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()}
        == {6, 35}
    ]
    product_carbon_bromine = [
        bond
        for bond in product_bonds
        if {bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()}
        == {6, 35}
    ]

    reactant_ring_count = mol_reactant_list.GetRingInfo().NumRings()
    product_ring_count = mol_product.GetRingInfo().NumRings()
    reactant_heterocycle_count = Descriptors.NumAliphaticHeterocycles(
        mol_reactant_list
    ) + Descriptors.NumAromaticHeterocycles(mol_reactant_list)
    product_heterocycle_count = Descriptors.NumAliphaticHeterocycles(
        mol_product
    ) + Descriptors.NumAromaticHeterocycles(mol_product)

    # ~ reactant_amide_count = reactant.count("NC=O") + reactant.count("C(=O)N")
    # ~ product_amide_count = product.count("NC=O") + product.count("C(=O)N")

    reactant_urea_count = len(mol_reactant_list.GetSubstructMatches(urea_patt))
    product_urea_count = len(mol_product.GetSubstructMatches(urea_patt))

    reactant_ester_count = len(mol_reactant_list.GetSubstructMatches(ester_patt))
    product_ester_count = len(mol_product.GetSubstructMatches(ester_patt))

    reactant_cooh_count = len(mol_reactant_list.GetSubstructMatches(cooh_patt))
    product_cooh_count = len(mol_product.GetSubstructMatches(cooh_patt))

    reactant_amide_count = len(mol_reactant_list.GetSubstructMatches(amide_patt))
    product_amide_count = len(mol_product.GetSubstructMatches(amide_patt))

    reactant_carbamate_count = len(
        mol_reactant_list.GetSubstructMatches(carbamate_patt)
    )
    product_carbamate_count = len(mol_product.GetSubstructMatches(carbamate_patt))

    # ~ reactant_weinrab_count = len(mol_reactant_list.GetSubstructMatches(weinrab_patt))
    # ~ product_weinrab_count = len(mol_product.GetSubstructMatches(weinrab_patt))

    reactant_sulfonamide_count = len(
        mol_reactant_list.GetSubstructMatches(sulfonamide_patt)
    )
    product_sulfonamide_count = len(mol_product.GetSubstructMatches(sulfonamide_patt))

    reactant_ketone_count = len(mol_reactant_list.GetSubstructMatches(ketone_patt))
    product_ketone_count = len(mol_product.GetSubstructMatches(ketone_patt))

    reactant_aldehyde_count = len(mol_reactant_list.GetSubstructMatches(aldehyde_patt))
    product_aldehyde_count = len(mol_product.GetSubstructMatches(aldehyde_patt))

    reactant_alcohol_count = len(mol_reactant_list.GetSubstructMatches(alcohol_patt))
    product_alcohol_count = len(mol_product.GetSubstructMatches(alcohol_patt))

    reactant_arom_nit_count = len(mol_reactant_list.GetSubstructMatches(arom_nit_patt))
    product_arom_nit_count = len(mol_product.GetSubstructMatches(arom_nit_patt))

    reactant_ether_count = len(mol_reactant_list.GetSubstructMatches(ether_patt))
    product_ether_count = len(mol_product.GetSubstructMatches(ether_patt))

    reactant_carbonyl_count = len(mol_reactant_list.GetSubstructMatches(carbonyl_patt))
    product_carbonyl_count = len(mol_product.GetSubstructMatches(carbonyl_patt))

    reactant_amine_count = len(mol_reactant_list.GetSubstructMatches(amine_patt))
    product_amine_count = len(mol_product.GetSubstructMatches(amine_patt))

    reactant_nitroso_count = len(mol_reactant_list.GetSubstructMatches(nitroso_patt))
    product_nitroso_count = len(mol_product.GetSubstructMatches(nitroso_patt))

    # ~ if product_br_count > reactant_br_count or product_cl_count > reactant_cl_count:
    if product_ring_count > reactant_ring_count:
        if product_heterocycle_count > reactant_heterocycle_count:
            return "Heterocycle formation"
        return "Ring formation"
    if len(product_carbon_carbon) > len(reactant_carbon_carbon):
        return "C-C bond formation"
    # if rxn_diffs.product_dissim == "Br" and rxn_diffs.reactant_dissim == "":
    #     return "Addition of Bromine"
    # if rxn_diffs.product_dissim == "Cl" and rxn_diffs.reactant_dissim == "":
    #     return "Addition of Chlorine"
    # ~ if product_ring_count < reactant_ring_count:
    # ~     return "Ring destruction"
    if product_urea_count > reactant_urea_count:
        return "Urea formation"
    if product_ester_count > reactant_ester_count:
        return "Ester formation"
    if product_amide_count > reactant_amide_count:
        return "Amide formation"
    if product_sulfonamide_count > reactant_sulfonamide_count:
        return "Sulfonamide formation"
    if product_ether_count > reactant_ether_count:
        return "Ether formation"
    if product_arom_nit_count > reactant_arom_nit_count:
        return "Electrophilic aromatic substitution"
    if product_ester_count < reactant_ester_count:
        if (
            product_aldehyde_count > reactant_aldehyde_count
            or product_alcohol_count > reactant_alcohol_count
        ):
            return "Reduction"  # temporary
        return "Ester hydrolysis"
    if len(product_carbon_chlorine) > len(reactant_carbon_chlorine) or len(
        product_carbon_bromine
    ) > len(reactant_carbon_bromine):
        return "C-X bond formation"
    if product_amide_count < reactant_amide_count:
        return "Amide hydrolysis"
    if product_ether_count < reactant_ether_count:
        if (
            product_ether_count + product_alcohol_count + product_carbonyl_count
            == reactant_ether_count + reactant_alcohol_count + reactant_carbonyl_count
        ):
            return "Ether cleavage"
    if c_prod_arom_arom[frozenset((6,))] - c_react_arom_arom[frozenset((6,))]:
        return "C-C bond formation"
    if c_prod_arom_attach - c_react_arom_attach:
        return "Nucleophilic aromatic substitution"
    if len(product_carbon_nitrogen) > len(reactant_carbon_nitrogen):
        if product_carbamate_count > reactant_carbamate_count:
            return "Protection"
        if product_amine_count > reactant_amine_count:
            return "Amine alkylation"
        if product_nitroso_count > reactant_nitroso_count:
            return "Nitrosation"
        return "C-N bond formation"
    if len(product_carbon_sulfur) > len(reactant_carbon_sulfur):
        return "C-S bond formation"
    if product_double_bond_count == reactant_double_bond_count - 1:
        if rxn_diffs.product_dissim == rxn_diffs.reactant_dissim:
            return "Reduction"
        elif rxn_diffs.product_dissim == "":
            return "Deprotection"
    if (
        product_cooh_count > reactant_cooh_count
        or product_aldehyde_count > reactant_aldehyde_count
        or product_ketone_count > reactant_ketone_count
    ):
        if (
            product_cooh_count
            + product_aldehyde_count
            + product_alcohol_count
            + product_ketone_count
            == reactant_cooh_count
            + reactant_aldehyde_count
            + reactant_alcohol_count
            + reactant_ketone_count
        ):
            return "Oxidation"

    return "UNSPECIFIED"

product_blacklist = {"Cl", "[Cl-]", "[Na]Cl", "Br", "[Br-]", "[K]", "[K+]", "I", "[I-]", "O"}

@click.command()
@click.argument("src_file", type=click.Path(exists=True))
@click.argument("tgt_file", type=click.Path(exists=True))
# ~ @click.option("-o", "out", type=click.Path(), required=True)
def main(src_file, tgt_file):
    """Main command that click uses to execute everything."""

    with open(src_file) as f_src, open(tgt_file) as f_tgt:
        src_lines = [line.strip() for line in f_src]
        tgt_lines = [line.strip() for line in f_tgt]
    tok_reactions = zip(src_lines, repeat(" > "), tgt_lines)
    reactions: List[str] = ["".join("".join(rxn).split()) for rxn in tok_reactions]
    reaction_types = []
    for index, rxn_string in tqdm(enumerate(reactions)):
        # rxn_line = [rxn_string]
        # ~ click.echo(f"Reaction {i+1}: '{rxn_string}'")
        reaction = split_into_rxn_class(rxn_string)
        # ~ click.echo(type(reaction.mol_reactants[0]))
        # ~ click.echo(type(reaction.fp_reactants[0]))
        reaction_diffs = find_reaction_difference(reaction)
        # ~ click.echo(reaction.smiles())
        # ~ click.echo(reaction_diffs)
        for diff_set in reaction_diffs:
            if diff_set.product in product_blacklist:
                continue
            # ~ click.echo(f"\tRelevant reactant: '{diff_set.rel_reactant}'")
            # ~ click.echo(f"\t\tDissimilarity: '{diff_set.reactant_dissim}'")
            # ~ click.echo(f"\tProduct in question: '{diff_set.product}'")
            # ~ click.echo(f"\t\tDissimilarity: '{diff_set.product_dissim}'")
            rxn_type = classify_reaction(reaction, diff_set)
            reaction_types.append(rxn_type)
            rxn_with_correct_product = ">".join([rxn_string.rsplit(">", maxsplit=1)[0], diff_set.product])
            click.echo(" ".join([rxn_with_correct_product, rxn_type, str(index + 1)]))
            # ~ click.echo(f"\tReaction type: '{rxn_type}'")
        # click.echo(" ".join(rxn_line))
    click.echo()
    most_common_types = Counter(reaction_types).most_common()
    rtypes_dataset = tablib.Dataset(*most_common_types, headers=("Type", "Hits"))
    click.echo(rtypes_dataset)


if __name__ == "__main__":
    main()
