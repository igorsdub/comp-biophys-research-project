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
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())
def main(input_filepath, output_filepath):
    """ Make interim PDB structures (from pdb/raw/) 
        for final processing (saved in pdb/processed/).
    """
    logger = logging.getLogger(__name__)
    logger.info('making interim PDB structures from raw PDB files')

    config = utils.read_config()
    pdb_codes = config['pdb']['codeList']

    # Data import
    pdb_raw_structures = {pdb_code : load_structure(pdb_code, input_filepath, file_extension="pdb1") \
                        for pdb_code in pdb_codes}

    pdb_models = {pdb_code : get_models(PandasPdb_object) \
                    for pdb_code, PandasPdb_object in pdb_raw_structures.items()}

    # Data processing
    pdb_interim_structures = {}
    for pdb_code, PandasPdb_object in pdb_raw_structures.items():
        # Remove water
        PandasPdb_object.df['HETATM'] = PandasPdb_object.df['HETATM'][PandasPdb_object.df['HETATM']['residue_name'] != 'HOH']  
        # Rename chains
        PandasPdb_object.df['ATOM']['chain_id'] = rename_chains(PandasPdb_object.df['ATOM']['chain_id'], \
            number_of_protomers=pdb_models[pdb_code])
        PandasPdb_object.df['HETATM']['chain_id'] = rename_chains(PandasPdb_object.df['HETATM']['chain_id'], \
            number_of_protomers=pdb_models[pdb_code])

        pdb_interim_structures[pdb_code] = PandasPdb_object

    # Save interim data
    for pdb_code, PandasPdb_object in pdb_interim_structures.items():
        save_structure(PandasPdb_object, pdb_code, output_filepath, form_idx=None)

def load_structure(pdb_code, input_filepath, file_extension="pdb"):
    """ Loads PDB file inot BioPandas object.
    """
    pdb_filepath = os.path.join(input_filepath, "{}.{}".format(pdb_code, file_extension))

    return PandasPdb().read_pdb(pdb_filepath)

def get_models(data):
    """ Extracts number of MODEL records from BioPandas object.
    """
    return data.df['OTHERS'][data.df['OTHERS']['record_name'] == 'MODEL'].shape[0]

def rename_chains(data, number_of_protomers=2):
    """ Standartizes chain ID labels for homo-multi-mers.
        Note: all chains must be of the equal length. 
    """
    chain_id_labels = ['A', 'B', 'C', 'D', 'E', 'F']

    unique_chain_ids = data.unique()

    if number_of_protomers == 1 and unique_chain_ids[0] == "A":
        pass

    else:
        atoms_per_protomer = data.shape[0] / number_of_protomers
        
        for protomer_id in range(number_of_protomers):
            slice_start = int(protomer_id * atoms_per_protomer)
            slice_end = int((protomer_id + 1) * atoms_per_protomer + 1)
            data.iloc[slice_start:slice_end] = chain_id_labels[protomer_id]
    
    return data

def save_structure(data, pdb_code, output_filepath, form_idx=None):
    """ Save BioPandas object as a PDB record file.
    """
    if form_idx == None:
        output_filename = "{}.pdb".format(pdb_code)
    
    else:
        output_filename = "{}.{}.pdb".format(pdb_code, form_idx)

    data.to_pdb(path=os.path.join(output_filepath, output_filename),
                    records=['ATOM', 'HETATM', 'OTHERS'],
                    gz=False,
                    append_newline=True)
            
if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()
