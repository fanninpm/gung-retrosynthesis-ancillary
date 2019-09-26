# gung-retrosynthesis-ancillary
Ancillary code to facilitate prediction of retrosynthesis reactions.

The file `find_certain_reactions.py` currently (as of 24 September 2019) takes a single SMILES file as an input. The command `python3 find_certain_reactions.py FILE` outputs the reaction types next to each reaction that has a non-empty product.

For eventual training with Open_NMT, `src` and `tgt` files can be generated from a file of reactions in SMILES format using the commands `python3 data_prep_4.py --step2 --step3 --label=src FILE` and `python3 data_prep_4.py --step2 --step3 --label=tgt FILE`. Use the `--step1` option with `data_prep_4.py` only if the reactants and reagents are atom-mapped to the products.
