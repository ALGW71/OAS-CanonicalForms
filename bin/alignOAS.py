import pandas as pd
from pymol import cmd
import os

# this code just works on IDs
# Can you incorporate some of the structures you have built?

def align_oas_loops(query_seq, chain, cdr, root_dir, out_dir):
    # clear
    cmd.reinitialize()
    # load template
    cmd.load(f'{out_dir}/template.pdb', 'temp')
    # load query
    cmd.load(f'{root_dir}{query_seq}.pdb', 'query_input')
    cmd.remove('sidechain')
    cmd.remove('solvent')
    cmd.remove('hydrogens')
    cmd.remove('element O')

    if cdr == 1:
        all = '22-43'
        anchors = '22+23+24+25+26+39+40+41+42+43'
        save = '26-39'
    elif cdr == 2:
        all = '51-70'
        anchors = '51+52+53+54+55+66+67+68+69+70'
        save = '55-66'
    elif cdr == 3:
        all = '100-122'
        anchors = '100+101+102+103+104+118+119+120+121+122'
        save = '104-118'

    # select the region of interest
    cmd.select('query', f'query_input and chain {chain} and resi {all}')

    # align with super
    out = cmd.pair_fit(f'query and chain {chain} and resi {anchors} and n. CA', # query
            f'temp and chain {chain}  and resi {anchors} and n. CA') # target

    # select the region of interest
    cmd.select('query_align', f'query and chain {chain} and resi {save}')
    cmd.remove('solvent')
    cmd.remove('hydrogens')
    cmd.remove('sidechain')
    cmd.remove('element O')

    outname = query_seq.removesuffix(".pdb")
    cmd.save(f'{out_dir}{outname}.pdb', 'query_align')
    cmd.reinitialize()
    return (query_seq + ".pdb", round(out,3))


def oas_align(chain, cdr, length, num):
    # specify input folder
    root_dir = '../input/oas_structures/'
    out_dir = f'../output/aligned/chain_{chain}/cdr_{cdr}/length_{length}/'

    if chain == 'L':
    # define the columns
        col1 = f'cdr{cdr}_aa_light'
        col2 = f'cdr{cdr}l_len'
    elif chain == 'H':
        col1 = f'cdr{cdr}_aa_heavy'
        col2 = f'cdr{cdr}h_len'

    #===========================================================================
    # FILTERING
    # 1. Filter meta by structures that have been run through immune builder
    meta = pd.read_csv('../resources/csvs/paired_info_cdrs.csv.gz')

    # all available structures are those in the root directory...
    available_structures = os.listdir(root_dir)
    available_structures = [x.replace('.pdb', '') for x in available_structures if x.endswith('.pdb')]
    meta = meta[meta['ID'].isin(available_structures)]

    # 2. Filter by length
    meta = meta[['ID', col1, col2]] # cols
    meta = meta[meta[f'{col2}']==length] # length
    meta = meta[meta['ID'].notnull()] # not null

    # You had to swap these round now - get non redundant - then remove already plotted.
    # # 4. Filter out sequences which are repeated more than once.
    meta = meta.groupby(col1).head(1)
    print(f'XXX - Meta length after removing seqs duplicated: {len(meta)}')
    #===========================================================================

    # 3. Remove those files you have already plotted...
    # Read the MDS dist files
    # We may now have multiple samples
    mds_name = f'../output/plotted_data/vis_mds_{chain}_{cdr}_length_{length}.csv'
    print('XXX - ' + mds_name)
    try:
        already_plotted = pd.read_csv(mds_name)['PDB'].values.tolist()
        print("XXX - Meta length before removing files already plotted: " + str(len(meta)))
        meta = meta[~meta['ID'].isin(already_plotted)]
        print("XXX - Meta length after: " + str(len(meta)))
    except FileNotFoundError:
        print(f"XXX - MDS file '{mds_name}' not found. Skipping removal of already plotted files.")

    # # 4. Filter out sequences which are repeated more than once.
    meta = meta.groupby(col1).head(1)
    print(f'XXX - Meta length after removing seqs duplicated: {len(meta)}')
    #===========================================================================

    # Make list of files from the meta
    pdb_files = meta['ID'].tolist()
    pdb_files.sort()
    print(f'XXX - Number of files sent to align: {len(pdb_files[:num])}')
    rmsd_list = []
    for file in pdb_files[:num]:
        try:
                print(f'---------- Running {file} ----------')
                out = align_oas_loops(file, chain, cdr, root_dir, out_dir)
                rmsd_list.append(out)
        except:
                print(f'error with {file}')

    ls = sorted(rmsd_list, key=lambda x: x[1], reverse=True)
    rmsd_df = pd.DataFrame(ls, columns=['File', 'RMSD'])
    rmsd_df.to_csv(out_dir + "oas_" + 'rmsd_summary.csv', index=False)
