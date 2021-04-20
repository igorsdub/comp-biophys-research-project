# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv

import os
import glob
from numpy.core.fromnumeric import cumsum
import pandas as pd
import numpy as np
import src.utilities as utils

@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())
def main(input_filepath, output_filepath):
    """ Runs data processing scripts to turn external data (from data/external/) into
        processed data ready to be analysed (saved in data/processed/).
    """
    logger = logging.getLogger(__name__)
    logger.info('making processed 1-point scan data set from external data')

    config = utils.read_config()
    pdb_codes = config['pdb']['codeList']
    pdb_codes = ["1m9a"]

    # Get filepaths
    input_energy_filepaths = sorted(glob.glob(os.path.join(input_filepath, "*.1point.energy.csv")))

    # Load data
    input_energy_data = {os.path.basename(filepath) : load_data(filepath) for filepath in input_energy_filepaths}

    # Apply cumsum and convert into long format
    long_energy_data = {}
    for filename, data in input_energy_data.items():
        cum_data = calculate_cumsum(data)
        long_energy_data[filename] = convert_to_long(cum_data)

    # Process data
    allostery_data = {}
    for pdb_code in pdb_codes:
        data = {form_idx : long_energy_data["{}.{}.1point.energy.csv".format(pdb_code, form_idx)] for form_idx in range(3)}

        collated_energy_data = collate_free_energy(data)
        allostery = calculate_allostery(collated_energy_data)

        
        allostery_data[pdb_code] = allostery

    # Save data
    # for filename, data in long_energy_data.items():
    #     save_data(data.astype({'spring_strength':str}), os.path.join(output_filepath, filename))

    for pdb_code, data in allostery_data.items():
        for mode in range(25,101,25):
            filename = "{}.1point.allostery.m{:03}.csv".format(pdb_code, mode)
            save_data(data[data['mode'] == mode].astype({'spring_strength':str}), os.path.join(output_filepath, filename))

def load_data(filepath):
    """ Load data into dataframe.
    """
    
    return  pd.read_csv(filepath, delimiter=",", comment="#", header=0, dtype={'spring_strength' : float, 'mode' : np.int16})

def calculate_cumsum(input_data):
    """ Performs cumulative sum on number of array sections 
        and returns array of the same size.
    """

    cum_data = input_data.copy()
    column_names = input_data.columns.to_list()
    residue_numbers = column_names[2:]

    split_arrays = np.split(input_data[residue_numbers].to_numpy(), len(input_data['spring_strength'].unique()))
    cum_arrays = np.cumsum(split_arrays, axis=1)
    cum_data = np.concatenate(cum_arrays, axis=0)

    output_data = pd.DataFrame(data=np.concatenate([input_data[['mode', 'spring_strength']].to_numpy(), cum_data], \
        axis=1), columns=column_names)

    return output_data


def convert_to_long(input_data):
    """ Process 1-point scan results for different spring strenghts into
        single long-format dataframe.
    """

    residue_numbers = input_data.columns.to_list()[2:]
    # Convert wide to long
    output = pd.melt(input_data, id_vars=['mode','spring_strength'], \
        value_vars=residue_numbers, var_name='residue_number', value_name='free_energy')
    
    # Ordered and sort dataframe
    output = output[['residue_number', 'spring_strength', 'mode', 'free_energy']]\
        .astype({'spring_strength': np.float32, 'residue_number' : np.int16, 'mode' : np.int16})

    return output.sort_values(by=['residue_number', 'spring_strength', 'mode'])

def collate_free_energy(input_data):
    """ Creates dataframe with free energy from diffrent structural forms.
        {0:apo , 1:holo1, 2:holo2}
    """

    output_data = input_data[0]['residue_number'].to_frame().copy()
    output_data['spring_strength'] = input_data[0]['spring_strength']
    output_data['mode'] = input_data[0]['mode']
    
    for form_idx, df in input_data.items():
        output_data['G_{}'.format(form_idx)] = df['free_energy']
                                                               
    return output_data

def calculate_allostery(input_data):
    """ Calculate free energy differences and allostery.
    """
    output_data = input_data.copy()
    
    output_data['dG_1'] = output_data['G_1'] - output_data['G_0']
    output_data['dG_2'] = output_data['G_2'] - output_data['G_1']
    output_data['ddG'] = output_data['dG_2'] - output_data['dG_1']
    output_data['allostery'] = np.exp(output_data['ddG'])
    
    return output_data

def save_data(input_data, path):
    """ Save dataframe to a specified .csv file. 
    """
    input_data.to_csv(path, index=False, float_format='%.5f')
    
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
