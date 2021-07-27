# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv

import glob, os
from os.path import join as join_paths
import pandas as pd
import numpy as np
import src.utilities as utils
from shutil import copy

@click.command()
@click.argument('input_dir', type=click.Path(exists=True))
@click.argument('output_dir', type=click.Path())
def main_commandline(input_dir, output_dir):
    """ Runs data processing scripts to turn interim data (from data/interim/) into
        processed data ready to be analysed (saved in data/processed/).
        Commandline function with Click.
    """
    logger = logging.getLogger(__name__)
    logger.info('making processed data set from interim data')

    config = utils.read_config()
    pdb_codes = config['pdb']['codeList']
    
    # Get paths
    eigenvalues_paths = sorted(glob.glob(os.path.join(input_dir, "*.eigenvalues")))
    # bfactors_paths = sorted(glob.glob(os.path.join(input_dir, "*.mode.m025.bfactors")))
    energy_paths = sorted(glob.glob(os.path.join(input_dir, "*.mode.energy")))
    # frequencies_paths = sorted(glob.glob(os.path.join(input_dir, "*.mode.frequencies")))

    # Load interim data
    eigenvalues = {os.path.basename(path) : load_data(path) for path in eigenvalues_paths}
    # interim_bfactors = {path.replace(input_dir, "") : load_data(path) for path in bfactors_paths}
    energy = {os.path.basename(path) : load_data(path) for path in energy_paths}
    # interim_frequencies = {path.replace(input_dir, "") : load_data(path) for path in frequencies_paths}

    # Restructure dictionary with energy dataframes
    restruct_eigenvalues = {}
    for pdb_code in pdb_codes:
        eigenvalues_dict = {}

        for form_idx in range(3):
            filename = "{}.{}.eigenvalues".format(pdb_code, form_idx)
            eigenvalues_dict[form_idx] = eigenvalues[filename]

        restruct_eigenvalues[pdb_code] = eigenvalues_dict

    restruct_energy = {}
    for pdb_code in pdb_codes:
        energy_dict = {}

        for form_idx in range(3):
            filename = "{}.{}.mode.energy".format(pdb_code, form_idx)
            energy_dict[form_idx] = energy[filename]

        restruct_energy[pdb_code] = energy_dict

    # Process data
    entropy = {}
    cooperativities = {}
    for pdb_code in pdb_codes:
        entropy[pdb_code] = collate_entropy(restruct_energy[pdb_code])
        # Cooperativity
        combined_eigenvalues = collate_eigenvalues(restruct_eigenvalues[pdb_code])

        cooperativity = calculate_cooperativity(combined_eigenvalues)

        cooperativities[pdb_code] = cooperativity

    # Calcualte cooperativity using the classical limit
    # Thomas Rodgers alorithm from 2015 JBC study
    # subprocess.call(['bash', 'src/data/calculate_cooperativity.sh', '1m9a'])

    # Save data
    for pdb_code in pdb_codes:
        save_data(cooperativities[pdb_code], "{}.cooperativity".format(pdb_code), output_dir)
        save_data(entropy[pdb_code], "{}.entropy".format(pdb_code), output_dir)

    # Copy files
    for path in glob.glob(os.path.join(input_dir, "*.CAonly.pdb")):
        copy(path, output_dir)

    for path in glob.glob(os.path.join(input_dir, "*.draw_enm.pml")):
        copy(path, output_dir)

#########################################################################

def main(input_dir, output_dir):
    """ Runs data processing scripts to turn interim data (from data/interim/) into
        processed data ready to be analysed (saved in data/processed/).
    """

    config = utils.read_config()
    pdb_codes = config['pdb']['codeList']
    
    # Get paths
    # Directory path example: "data/raw/-c09.50/-mass-ca-het/0"
    cutoff_paths = glob.glob(join_paths(input_dir, "*"), recursive=True)
    for cutoff_path in cutoff_paths:
        flag_paths = glob.glob(join_paths(cutoff_path, "*"), recursive=True)
        for flag_path in flag_paths:
            cutoff_flag = os.path.basename(cutoff_path)
            other_flags = os.path.basename(flag_path)
            
            # Read martix.eigenfacs files
            idxs = ["0", "1", "2"]
            a_files = [join_paths(flag_path, idx, "matrix.eigenfacs") for idx in idxs]
            a_exist = [f for f in a_files if os.path.isfile(f)]

            if a_files.sort() == a_exist.sort():
                eigenfacs_0 = read_file(join_paths(flag_path, "0", "matrix.eigenfacs"))
                eigenfacs_1 = read_file(join_paths(flag_path, "1", "matrix.eigenfacs"))
                eigenfacs_2 = read_file(join_paths(flag_path, "2", "matrix.eigenfacs"))
            else:
                print("matrix.eigenfacs might be missing")
                return None

            # Create  DataFrames with eigenvalues 
            eigenvals_0 = extract_eigenvals(eigenfacs_0)
            eigenvals_1 = extract_eigenvals(eigenfacs_1)
            eigenvals_2 = extract_eigenvals(eigenfacs_2)

            # Move all eigenvalues into one DataFrame
            eigenvals_all = eigenvals_0.copy()
            eigenvals_all.rename(columns={"eigenvalue": "eigenvalue_0"}, inplace=True)
            
            eigenvals_all['eigenvalue_1'] = eigenvals_1['eigenvalue'][eigenvals_1.index \
                .isin(eigenvals_all.index)]
            eigenvals_all['eigenvalue_2'] = eigenvals_2['eigenvalue'][eigenvals_2.index \
                .isin(eigenvals_all.index)]

            # Calculate dissociation constants and cooperativity
            mode_number = eigenvals_all.index

            diss_consts = pd.DataFrame(index=mode_number)
            diss_consts['K_1'] = (eigenvals_all['eigenvalue_1'] / eigenvals_all['eigenvalue_0']).to_numpy()
            diss_consts['K_2'] = eigenvals_all['eigenvalue_2'] / eigenvals_all['eigenvalue_1']
            
            coop = pd.DataFrame(index=mode_number)
            coop['coop'] = (eigenvals_all['eigenvalue_2'] * eigenvals_all['eigenvalue_0']) / \
                        (eigenvals_all['eigenvalue_1'] ** 2)

            # Calcualte cumulative (total) values
            eigenvals_cum = eigenvals_all.copy()
            eigenvals_cum[:] = np.log(eigenvals_all[:]).cumsum()
            diss_consts_cum = diss_consts.copy()
            diss_consts_cum[:] = np.cumprod(diss_consts_cum[:])
            coop_cum = coop.copy()
            coop_cum[:] = np.cumprod(coop[:])

            # Create results directory
            output_subdir = join_paths(output_dir, cutoff_flag, other_flags)
            os.makedirs(output_subdir, exist_ok=True)

            # Save data
            eigenvals_cum.to_csv(join_paths(output_subdir, "eigenvals.csv"))
            diss_consts_cum.to_csv(join_paths(output_subdir, "diss_consts.csv"))
            coop_cum.to_csv(join_paths(output_subdir, "coop.csv"))
            
    

def read_file(filepath):
    """ Read file line by line.
    """
    with open(filepath) as file: # Use file to refer to the file object
        line_list = file.read().splitlines()
    
    return line_list

def extract_eigenvals(eigenfacs):
    """ Extracts eigenvalues from imported matrix.eigenfacs file
        into DataFrame.
    """
    if eigenfacs[0][1:7] != "VECTOR":
        print("Check matrix.eigenfacs file.")
    else:
        i = 0
        # First line: VECTOR 
        # Second line: ---
        while(eigenfacs[i+2][1:7] != "VECTOR"):
            i += 1
        no_beads = i

    eigenvals = eigenfacs[::no_beads+2]
    mode_numbers = [int(line[8:12]) for line in eigenvals]
    eigenvals = [float(line[-10:]) for line in eigenvals]

    eigenvals = pd.DataFrame(data=eigenvals, index=mode_numbers)
    eigenvals.columns = ['eigenvalue']
    eigenvals.index.name = 'mode_number'

    return eigenvals

def extract_eigenvecs(eigenfacs):
    """ Extracts eigenvectors from imported matrix.eigenfacs file
        into DataFrame.
    """
    if eigenfacs[0][1:7] != "VECTOR":
        print("Check matrix.eigenfacs file.")
    else:
        i = 0
        # First line: VECTOR 
        # Second line: ---
        while(eigenfacs[i+2][1:7] != "VECTOR"):
            i += 1
        no_beads = i

    eigenvals = eigenfacs[::no_beads+2]
    mode_numbers = [int(line[8:12]) for line in eigenvals]
    no_modes = len(mode_numbers)
    eigenvecs = np.split(np.array(eigenfacs), no_modes)
    # Remove VECTOR and --- lines
    eigenvecs = [array[2:] for array in eigenvecs]
    eigenvecs = np.concatenate(eigenvecs)
    eigenvecs = np.loadtxt(eigenvecs)

    mode_col = np.repeat(mode_numbers, no_beads)
    bead_col = np.tile(np.arange(no_beads)+1, no_modes)
    eigenvecs = np.hstack((np.vstack((mode_col, bead_col)).T, eigenvecs))

    return eigenvecs



def load_data(path):
    """ Load interim data into dataframe.
    """
    # Load data
    interim_data = pd.read_csv(path, sep=',', header=0)

    return interim_data

def collate_eigenvalues(input_data, trivial_modes=True):
    """ Creates dataframe with free energy from diffrent structural forms.
        {0:apo , 1:holo1, 2:holo2}
    """
    if trivial_modes:
        output_data = input_data[0]['mode_number'][input_data[0]['mode_number'] > 6].to_frame()
    else:
        output_data = input_data[0]['mode_number'].to_frame()

    for form_idx, df in input_data.items():
        output_data['eigenvalue_{}'.format(form_idx)] = df['eigenvalue'][df['mode_number'].isin(output_data['mode_number'])]
        
    if trivial_modes:
        output_data['mode_number'] = range(1, output_data.shape[0]+1)
        output_data.reset_index(drop=True, inplace=True)
                                                               
    return output_data



def calculate_cooperativity(input_data, cumulative=True, from_eigenvalues=True):
    """ Calculate dissociation constants constants (free eneregy changes) and cooperativity
        from eigenvalues or free energy.
    """
    output_data = input_data.copy()
    if from_eigenvalues:
        output_data.loc[:,'eigenvalue_0':] = output_data.loc[:,'eigenvalue_0':]
        
        output_data['K_1'] = output_data['eigenvalue_1'] / output_data['eigenvalue_0']
        output_data['K_2'] = output_data['eigenvalue_2'] / output_data['eigenvalue_1']
        output_data['cooperativity'] = (output_data['eigenvalue_2'] * output_data['eigenvalue_0']) / (output_data['eigenvalue_1'] ** 2)

        if cumulative:
            output_data.loc[:,'eigenvalue_0':'eigenvalue_2'] = np.cumsum(np.log(output_data.loc[:,'eigenvalue_0':'eigenvalue_2']))
            output_data.loc[:,'K_1':] = np.cumprod(output_data.loc[:,'K_1':])
            
    else:
        # Calculate allostery from free energy (from DDPT FREQEN module)
        if cumulative:
            output_data.loc[:,'G_0':] = np.cumprod(output_data.loc[:,'G_0':])
    
        output_data['dG_1'] = output_data['G_1'] - output_data['G_0']
        output_data['dG_2'] = output_data['G_2'] - output_data['G_1']
        output_data['ddG'] = output_data['dG_2'] - output_data['dG_1']
        output_data['cooperativity'] = np.exp(output_data['ddG'])
    
    return output_data

def save_data(input_data, filename, output_dir):
    """ Save dataframe to a specified .csv file. 
    """
    input_data.to_csv(path_or_buf=os.path.join(output_dir, filename), index=False, float_format='%g')
    
    return None

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()
