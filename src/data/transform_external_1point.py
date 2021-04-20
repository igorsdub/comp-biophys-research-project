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
    """ Runs data processing scripts to turn external 1-point mutational 
        scan data (from data/external/) into interim data ready to be 
        processed (saved in data/interim/).
    """
    logger = logging.getLogger(__name__)
    logger.info('making interim 1-point data from external data')

    # config = utils.read_config()
    # pdb_codes = config['pdb']['codeList']
    spring_strenghts = [0.25, 0.313, 0.375, 0.438, 0.5, 0.563, 0.625, 0.688, 0.75, \
            0.813, 0.875, 0.938, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00, \
            3.25, 3.50, 3.75, 4.00]
    # Get filepaths
    raw_energy_filepaths = sorted(glob.glob(os.path.join(input_filepath, "*.1point.energy.csv")))

    # Load data
    raw_energy_data = {os.path.basename(filepath) : load_data(filepath) for filepath in raw_energy_filepaths}

    # Process data
    interim_energy_data = {}
    for filename, data in raw_energy_data.items():
        splitted_data = split_data(data, spring_strenghts, trivial_modes=False)
        wide_data = process_to_wide(splitted_data)

        interim_energy_data[filename] = wide_data

    # Save data
    for filename, data in interim_energy_data.items():
        data.to_csv(os.path.join(output_filepath, filename), index=False)

# Custom fucntions
def load_data(filepath):
    """ Load 1-point mutational scan concatenated into dataframe.
    """
    
    return  pd.read_csv(filepath, delimiter=",", comment="#", header=0)

def split_data(data, spring_strengths, trivial_modes=True):
    """ Split raw dataframe into a long format dataframes for different.
    """
    
    # The first row contains residue sequnce numbers
    res_nums = data.columns.to_list()

    # Extract energies starting from the second row
    splitted_data = np.split(data.to_numpy(), len(spring_strengths), axis=0)
    
    output = []
    for counter, spring_strength in enumerate(spring_strengths):
        df =  pd.DataFrame(splitted_data[counter])
        df.columns = res_nums

        if trivial_modes:
            # Include only non-trivial modes, i.e. startfing from 7-th mode
            df = df.iloc[6:]
        
        df['mode'] = np.arange(df.shape[0]) + 1
        df['spring_strength'] = df.shape[0] * [spring_strength]

        # Rearrange columns
        column_names = ['mode', 'spring_strength'] + res_nums
        output.append(df[column_names])
    
    return output

def process_to_wide(data):
    """ Process 1-point scan results for different spring strenghts into
        single wide-format dataframe.
    """

    output = pd.concat(data, axis=0, ignore_index=True) # Semi-wide format

    return output

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()
