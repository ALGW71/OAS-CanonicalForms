# the template must be the same for both functions
# and from SAbDab high resolution ~1.5A so it is real
import pandas as pd
from pymol import cmd
import os

# put the template in the same directories as the query
# make a base template from the chain to align to:
def template(chain, cdr, length):
    meta = pd.read_csv(f'../resources/csvs/{chain}{cdr}_cdrs.csv')
    out_dir = f'../output/aligned/chain_{chain}/cdr_{cdr}/length_{length}/'
    # make the output directory
    os.makedirs(out_dir, exist_ok=True)

    col1 = f'cdr{cdr}'
    col2 = f'{col1}_len'

    # read in the metadata
    meta = pd.read_csv(f'../resources/csvs/{chain}{cdr}_cdrs.csv')
    meta = meta[meta[col2] == length]
    meta = meta[meta['is_representative'] == 1]
    # cluster must not have * in it
    meta = meta[meta['cluster'].str.contains('\*') == False]
    meta.sort_values(by='resolution', inplace=True)
    print(meta.head(3))

    # # highest resolution & representative structure
    highest_res = meta['pdb'].iloc[0]
    sabdab_chain = meta['dunbrack_chain_id'].iloc[0]
    print(f'XXX - Using {highest_res} as template for chain {sabdab_chain}, CDR{cdr}, Length = {length}')
    if cdr == 1:
        all = '22-43'
    elif cdr == 2:
        all = '51-70'
    elif cdr == 3:
        all = '100-122'

    # make the template
    cmd.reinitialize()
    cmd.load(f'../input/old_sabdab/{highest_res}.pdb', 'temp_in')
    cmd.remove('solvent')
    cmd.remove('hydrogens')
    cmd.remove('sidechain')
    cmd.remove('element O')

    # change back to home
    cmd.select('temp', f'temp_in and chain {sabdab_chain} and resi {all}')
    cmd.alter('temp', f'chain="{chain}"')
    cmd.save(f'{out_dir}/template.pdb', 'temp')
    cmd.reinitialize()
