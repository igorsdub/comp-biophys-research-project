# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv

import os
import glob
import pandas as pd
import src.utilities as utils

@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())
def main(input_filepath, output_filepath):
    """ Runs data processing scripts to turn raw data (from data/raw/) into
        interim data ready to be processed (saved in data/interim/).
    """
    logger = logging.getLogger(__name__)
    logger.info('making interim data set from raw data')

    # config = utils.read_config()
    # pdb_codes = config['pdb']['codeList']
    
    # Get filepaths
    raw_bfactors_filepaths = sorted(glob.glob(os.path.join(input_filepath, "*.mode.m025.bfactors")))
    raw_energy_filepaths = sorted(glob.glob(os.path.join(input_filepath, "*.mode.energy")))
    raw_frequencies_filepaths = sorted(glob.glob(os.path.join(input_filepath, "*.mode.frequencies")))

    # Load and process raw data
    interim_bfactors = {os.path.basename(filepath) : make_bfactors(filepath) for filepath in raw_bfactors_filepaths}
    interim_energy = {os.path.basename(filepath) : make_energy(filepath) for filepath in raw_energy_filepaths}
    interim_frequencies = {os.path.basename(filepath) : make_frequencies(filepath) for filepath in raw_frequencies_filepaths}

    # Save interim data
    save_data(interim_bfactors, output_filepath)
    save_data(interim_energy, output_filepath)
    save_data(interim_frequencies, output_filepath)



def make_bfactors(filepath):
    """ Process raw `mode.bfactors` file into .csv file.
    """
    column_names = ['record_name', 'atom_number', 'atom_name', 'residue_name', 'chain_id', 'residue_number', \
                'bfactor_predicted', 'bfactor_full_scaled', 'bfactor_scaled', 'bfactor_experimental']
    # Load data
    interim_data = pd.read_csv(filepath, header=None, skiprows=[0,1], delim_whitespace=True)
    # Rename columns
    interim_data.columns = column_names

    return interim_data

def make_energy(filepath):
    """ Process raw `mode.energy` file into .csv file.
    """
    column_names = ['mode', 'partition_function', 'free_energy', 'free_energy_sch', 'entropy', 'entropy_sch']
    # Load data
    interim_data = pd.read_csv(filepath, header=None, skiprows=[0], delim_whitespace=True) 
    # Rename columns
    interim_data.columns = column_names

    return interim_data

def make_frequencies(filepath):
    """ Process raw `mode.frequencies` file into .csv file.
    """
    column_names = ['frequencies']
    # Load data
    interim_data = pd.read_csv(filepath, header=None, skiprows=[0,1], delim_whitespace=True) 
    # Rename columns
    interim_data.columns = column_names

    return interim_data

def save_data(data, output_filepath):
    """ Saves dataframes with a specified filename from a dictionary. 
        Dictionary format : {"filename_1" : df_1, ...}
    """
    for filename, dataframe in data.items():
        dataframe.to_csv(path_or_buf=os.path.join(output_filepath, filename),index=False)
    
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
