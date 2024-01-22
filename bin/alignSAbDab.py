import pandas as pd
from pymol import cmd

def align_sabdab_loops(file, chain, pdb_chain_id, cdr, root_dir, out_dir):
    # clear
    cmd.reinitialize()
    # load template - this must be a new file that is the same for SAbDab and OAS...
    cmd.load(f'{out_dir}/template.pdb', 'temp')
    #load query
    cmd.load(f'{root_dir}{file}.pdb', 'query_input')
    cmd.remove('solvent')
    cmd.remove('hydrogens')
    cmd.remove('sidechain')
    cmd.remove('element O')

    # run pairfit at anchors
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
    # ==========================================================================================

    # select the region of interest
    cmd.select('query', f'query_input and chain {pdb_chain_id} and resi {all}')

    # align with super
    out = cmd.pair_fit(f'query and chain {pdb_chain_id} and resi {anchors} and n. CA', # query
            f'temp and chain {chain}  and resi {anchors} and n. CA') # target

    # select and save the region of interest
    cmd.select('query_align', f'query and chain {pdb_chain_id} and resi {save}')
    cmd.remove('solvent')
    cmd.remove('hydrogens')
    cmd.remove('sidechain')
    cmd.remove('element O')

    cmd.save(f'{out_dir}{file}_chain{pdb_chain_id}.pdb', 'query_align')
    cmd.reinitialize()

    return (file + "_chain" + pdb_chain_id + ".pdb", round(out,3))


def sabdab_align(chain, cdr, length, num):
    # specify i/o folders
    root_dir = '../input/sabdab_chains_renumb/'
    out_dir = f'../output/aligned/chain_{chain}/cdr_{cdr}/length_{length}/'

    # define the columns
    col1 = f'cdr{cdr}'
    col2 = f'{col1}_len'

    # find pdbs files in meta
    meta = pd.read_csv(f'../resources/csvs/{chain}{cdr}_cdrs.csv')
    meta = meta[meta[col2] == length]

    max_sabdab = len(meta[:num])
    rmsd_list = []
    for i in range(max_sabdab):
        file = meta['file'].iloc[i]
        # remove pdb extension
        file = file[:-4]
        pdb_chain_id = meta['chain'].iloc[i]
        try:
            print(f'---------- Running {file}: chain {pdb_chain_id} ----------')
            out = align_sabdab_loops(file, chain, pdb_chain_id, cdr, root_dir, out_dir)
            rmsd_list.append(out)
        except:
            print(f'error with {file}')

    ls = sorted(rmsd_list, key=lambda x: x[1], reverse=True)
    rmsd_df = pd.DataFrame(ls, columns=['File', 'RMSD'])
    rmsd_df.to_csv(out_dir + "sabdab_" + 'rmsd_summary.csv', index=False)