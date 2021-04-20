# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv

import os
import glob
import pandas as pd
import numpy as np
import src.utilities as utils

@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())
def main(input_filepath, output_filepath):
    """ Runs data processing scripts to turn interim data (from data/interim/) into
        processed data ready to be analysed (saved in data/processed/).
    """
    logger = logging.getLogger(__name__)
    logger.info('making interim data set from raw data')

    config = utils.read_config()
    pdb_codes = config['pdb']['codeList']
    
    # Get filepaths
    # bfactors_filepaths = sorted(glob.glob(os.path.join(input_filepath, "*.mode.m025.bfactors")))
    energy_filepaths = sorted(glob.glob(os.path.join(input_filepath, "*.mode.energy")))
    # frequencies_filepaths = sorted(glob.glob(os.path.join(input_filepath, "*.mode.frequencies")))

    # Load interim data
    # interim_bfactors = {filepath.replace(input_filepath, "") : load_data(filepath) for filepath in bfactors_filepaths}
    interim_energy = {os.path.basename(filepath) : load_data(filepath) for filepath in energy_filepaths}
    # interim_frequencies = {filepath.replace(input_filepath, "") : load_data(filepath) for filepath in frequencies_filepaths}

    # Restructure dictionary with energy dataframes
    restruct_interim_energy = {}
    for pdb_code in pdb_codes:
        energy_dict = {}

        for form_idx in range(3):
            filename = "{}.{}.mode.energy".format(pdb_code, form_idx)
            energy_dict[form_idx] = interim_energy[filename]

        restruct_interim_energy[pdb_code] = energy_dict

    # Process data
    processed_entropy = {}
    processed_energy = {}
    for pdb_code in pdb_codes:
        processed_entropy[pdb_code] = collate_entropy(restruct_interim_energy[pdb_code])
        free_energy = collate_free_energy(restruct_interim_energy[pdb_code])
        allostery = calculate_allostery(free_energy)
        processed_energy[pdb_code] = allostery

    # Save data
    for pdb_code in pdb_codes:
        save_data(processed_energy[pdb_code], "{}.allostery".format(pdb_code), output_filepath)
        save_data(processed_entropy[pdb_code], "{}.entropy".format(pdb_code), output_filepath)

def load_data(filepath):
    """ Load interim data into dataframe.
    """
    # Load data
    interim_data = pd.read_csv(filepath, sep=',', header=0)

    return interim_data

def collate_free_energy(data, trivial_modes=True):
    """ Creates dataframe with free energy from diffrent structural forms.
        {0:apo , 1:holo1, 2:holo2}
    """
    if trivial_modes:
        output_df = data[0]['mode'][data[0]['mode'] > 6].to_frame()

    else:
        output_df = data[0]['mode'].to_frame()
    
    for form_idx, df in data.items():
        output_df['G_{}'.format(form_idx)] = df['free_energy'][df['mode'].isin(output_df['mode'])]
    
    if trivial_modes:
        output_df['mode'] = range(1, output_df.shape[0]+1)
        output_df.reset_index(drop=True, inplace=True)
                                                               
    return output_df

def collate_entropy(data, trivial_modes=True):
    """ Creates dataframe with entropy from diffrent structural forms.
        {0:apo , 1:holo1, 2:holo2}
    """
    if trivial_modes:
        output_df = data[0]['mode'][data[0]['mode'] > 6].to_frame()

    else:
        output_df = data[0]['mode'].to_frame()
    
    for form_idx, df in data.items():
        output_df['S_{}'.format(form_idx)] = df['entropy'][df['mode'].isin(output_df['mode'])]
    
    if trivial_modes:
        output_df['mode'] = range(1, output_df.shape[0]+1)
        output_df.reset_index(drop=True, inplace=True)
                                                               
    return output_df

def calculate_allostery(data, cumulative=True):
    """ Calculate free energy differences and allostery.
    """
    output_df = data.copy()
    
    if cumulative:
        output_df.loc[:,'G_0':] = np.cumsum(output_df.loc[:,'G_0':])
    
    output_df['dG_1'] = output_df['G_1'] - output_df['G_0']
    output_df['dG_2'] = output_df['G_2'] - output_df['G_1']
    output_df['ddG'] = output_df['dG_2'] - output_df['dG_1']
    output_df['allostery'] = np.exp(output_df['ddG'])
    
    return output_df

def save_data(data, filename, output_filepath):
    """ Save dataframe to a specified .csv file. 
    """
    data.to_csv(path_or_buf=os.path.join(output_filepath, filename), index=False, float_format='%g')
    
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
