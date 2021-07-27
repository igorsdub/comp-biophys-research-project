# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv

import os
import pandas as pd
from biopandas.pdb import PandasPdb
from copy import deepcopy
import src.utilities as utils

@click.command()
@click.argument('input_dir', type=click.Path(exists=True))
@click.argument('output_dir', type=click.Path())
def main_commandline(input_dir, output_dir):
    """ Make interim PDB structures (from pdb/raw/) 
        for final processing (saved in pdb/processed/).
    """
    logger = logging.getLogger(__name__)
    logger.info('process raw PDB structure')
    main(input_dir, output_dir)

def main(input_dir, output_dir):
    """ Make interim PDB structures (from pdb/raw/) 
        for final processing (saved in pdb/processed/).
    """
    config = utils.read_config()
    pdb_code = config['pdb']['id']

    # Data import
    pdb_struct = load_structure(pdb_code, input_dir, file_extension="pdb1")                    

    # Data processing
    pdb_models = get_models(pdb_struct)

    # Remove water
    pdb_struct.df['HETATM'] = pdb_struct.df['HETATM'][pdb_struct.df['HETATM']['residue_name'] != 'HOH']
    # Remove HEZ
    pdb_struct.df['HETATM'] = pdb_struct.df['HETATM'][pdb_struct.df['HETATM']['residue_name'] != 'HEZ']

    # Select A form
    pdb_struct.df['ATOM'] = pdb_struct.df['ATOM'][pdb_struct.df['ATOM']['residue_name'] != 'BGLU']
    pdb_struct.df['ATOM']['occupancy'] = 1.00
    
    # Rename chains
    pdb_struct.df['ATOM']['chain_id'] = rename_chains(pdb_struct.df['ATOM']['chain_id'], \
        no_protomers=pdb_models)
    pdb_struct.df['HETATM']['chain_id'] = rename_chains(pdb_struct.df['HETATM']['chain_id'], \
        no_protomers=pdb_models)

    # Save data
    save_structure(pdb_struct, pdb_code, output_dir)

    return None

def load_structure(pdb_code, input_dir, file_extension="pdb"):
    """ Loads PDB file inot BioPandas object.
    """
    pdb_filepath = os.path.join(input_dir, "{}.{}".format(pdb_code, file_extension))

    return PandasPdb().read_pdb(pdb_filepath)

def get_models(data):
    """ Extracts number of MODEL records from BioPandas object.
    """
    return data.df['OTHERS'][data.df['OTHERS']['record_name'] == 'MODEL'].shape[0]

def rename_chains(data, no_protomers=2):
    """ Standartizes chain ID labels for homo-multi-mers.
        Note: all chains must be of the equal length. 
    """
    chain_id_labels = ['A', 'B', 'C', 'D', 'E', 'F']

    unique_chain_ids = data.unique()

    if no_protomers == 1 and unique_chain_ids[0] == "A":
        pass

    else:
        atoms_per_protomer = data.shape[0] / no_protomers
        
        for protomer_id in range(no_protomers):
            slice_start = int(protomer_id * atoms_per_protomer)
            slice_end = int((protomer_id + 1) * atoms_per_protomer + 1)
            data.iloc[slice_start:slice_end] = chain_id_labels[protomer_id]
    
    return data

def save_structure(data, pdb_code, output_dir):
    """ Save BioPandas object as a PDB record file.
    """
    output_filename = "{}.pdb".format(pdb_code)

    data.to_pdb(path=os.path.join(output_dir, output_filename),
                    records=['ATOM', 'HETATM'],
                    gz=False,
                    append_newline=True)

    return None
            
if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main_commandline()
