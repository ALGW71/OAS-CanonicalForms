import os
import pandas as pd
import numpy as np
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio import PDB

warnings.simplefilter('ignore', PDBConstructionWarning)

def read_pdb_files(root_dir, files):
    """Read PDB files and return a list of coordinates and atom count."""
    parser = PDB.PDBParser()

    aligned_pdbs = []
    for i in files:
        struct = parser.get_structure('aligned', f'{root_dir}{i}')
        coords = []
        atom_count = 0  # Counter for atoms

        for model in struct:
            for ch in model:
                for residue in ch:
                    for atom in residue:
                        coords.append(
                            np.concatenate([#np.array(residue.get_id()[1]),
                                            atom.get_coord()],
                                           axis=None)
                        )
                        atom_count += 1

        arr = np.stack(coords, axis=0)
        arr = arr.astype(np.double)
        aligned_pdbs.append((arr, atom_count))  # Append coordinates and atom count

    return aligned_pdbs

def calculate_rmsd(coordinates1, coordinates2):
    """Calculate the RMSD between two sets of coordinates."""
    num_atoms = len(coordinates1)
    squared_diff_sum = np.sum((coordinates1 - coordinates2) ** 2)
    rmsd = np.sqrt(squared_diff_sum / num_atoms)
    return rmsd

def calculate_rmsd_matrix(aligned_pdbs, files):
    """Calculate the distance matrix using RMSD."""
    num = len(files)
    distance_matrix_entries = []
    reduced_pdbs = list(aligned_pdbs[0:num])
    reduced_files = list(files[0:num])  # Convert tuple to list
    for i in range(num):
        n = len(reduced_pdbs)
        for j in range(n):
            try:
                rmsd = round(calculate_rmsd(aligned_pdbs[i][0], reduced_pdbs[j][0]),5)
            except Exception as e:
                rmsd = "Missing"
            distance_matrix_entries.append([reduced_files[0], reduced_files[j], rmsd])
        reduced_files.pop(0)
        reduced_pdbs.pop(0)

    return distance_matrix_entries



def save_distance_matrix(distance_matrix_entries, output_file):
    """Save the distance matrix to a CSV file."""
    df = pd.DataFrame(distance_matrix_entries, columns=['pdb1', 'pdb2', 'dist'])
    df.to_csv(output_file, index=False, header=True)

def run_dtw(chain, cdr, length):
    root_dir = f'../output/aligned/chain_{chain}/cdr_{cdr}/length_{length}/'

    # read the csv files and delete those files with more than 0.9A RMSD mismatch
    csvs = os.listdir(root_dir)
    csvs = [f for f in csvs if f.endswith('.csv')]

    files_to_keep = []
    # Loop through each csv file
    for csv in csvs:
        # Read in csv file as pandas dataframe
        df = pd.read_csv(os.path.join(root_dir, csv))

        # We also want to keep only the files that have a RMSD of less than 2.50
        filtered_df = df[df.iloc[:, 1] <= 1.5]
        files_to_keep += filtered_df['File'].values.tolist()

        # Filter dataframe by values greater than 2.25 RMSD in second column
        filtered_df = df[df.iloc[:, 1] > 1.5]
        # Get list of files to delete (from first column of filtered_df)
        files_to_delete = list(filtered_df.iloc[:, 0])
        print("XXX - Files to delete with poor RMSD: ", str(len(files_to_delete)))
        # Loop through files to delete and delete them
        for file in files_to_delete:
            try:
                os.remove(os.path.join(root_dir, file))
            except:
                pass

    files = files_to_keep
    files = [f for f in files if f.endswith('.pdb')]

    print("XXX - Number of files for DTW 1:", len(files))
    print("XXX - Number of files for DTW 2:", len(set(files)))

    files = list(set(files))
    files.sort()

    aligned_pdbs = read_pdb_files(root_dir, files)


    atom_counts = [pdb[1] for pdb in aligned_pdbs]  # Extract atom counts
    median_atom_count = np.median(atom_counts)  # Compute the median atom count
    print("XXX - Median atom count:", median_atom_count)

    print("XXX - Number of aligned PDBS:", len(aligned_pdbs))
    # Zip the files and aligned_pdbs lists together
    zipped = list(zip(files, aligned_pdbs))

    # Filter the pairs based on atom count
    
    # MOD HERE FOR CDRL3 Lengths 8-11====================================================
    filtered_pairs = [pair for pair in zipped if pair[1][1] == median_atom_count]
    # filtered_pairs = zipped
    # MOD HERE FOR CDRL3 Lengths 8-11====================================================

    # Unzip the filtered pairs back into separate lists
    files, aligned_pdbs = zip(*filtered_pairs)
    print("XXX - Number of after drop by median atom counts:", len(aligned_pdbs))

    print("XXX - Files read into nparray.")

    distance_matrix_entries = calculate_rmsd_matrix(aligned_pdbs, files)

    out_dir = f'../output/dist_mats/chain_{chain}/cdr_{cdr}/'
    os.makedirs(out_dir, exist_ok=True)
    out_file = f'{out_dir}/{chain}_{cdr}_length_{length}.csv'

    save_distance_matrix(distance_matrix_entries, out_file)
