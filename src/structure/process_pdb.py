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
    """ Proccesses interim PDB structures (from pdb/interim/) and creates PDB 
        structural forms for simulations (saved in pdb/processed/).
    """
    logger = logging.getLogger(__name__)
    logger.info('making processed PDB forms from interim PDB structures')
    main(input_dir, output_dir)

def main(input_dir, output_dir):
    """ Proccesses interim PDB structures (from pdb/interim/) and creates PDB 
        structural forms for simulations (saved in pdb/processed/).
    """
    config = utils.read_config()
    pdb_code = config['pdb']['id']

    # Data import
    pdb_struct = load_structure(pdb_code, input_dir, file_extension="pdb") 

    # Data processing
    # Delete residues 1 and 216
    # pdb_struct.df['ATOM'] = pdb_struct.df['ATOM'][(pdb_struct.df['ATOM']['residue_number'] != 1) \
    #     & (pdb_struct.df['ATOM']['residue_number'] != 216)]

    # Create structural forms
    pdb_0 = create_form(pdb_struct, form_idx=0)
    pdb_1 = create_form(pdb_struct, form_idx=1)
    pdb_2 = create_form(pdb_struct, form_idx=2)

    # Save processed data
    save_structure(pdb_0, 0, output_dir)
    save_structure(pdb_1, 1, output_dir)
    save_structure(pdb_2, 2, output_dir)



def load_structure(pdb_code, input_dir, file_extension="pdb"):
    """ Loads PDB file inot BioPandas object.
    """
    pdb_filepath = os.path.join(input_dir, "{}.{}".format(pdb_code, file_extension))
    return PandasPdb().read_pdb(pdb_filepath)

def create_form(data, form_idx=0):
    """ Creates PDB structure forms.
        form_idx = 0 is apo; 1 - holo1; and 2 - holo2
        Note: Only works for homodimers.
    """
    # Make a deep copy of BioPandas object to make changes
    data_out = deepcopy(data)
    
    # If form_idx == 2 that's holo2 already
    if form_idx == 1:
        hetatm_record_len = data_out.df['HETATM'].shape[0]
        # Keep only one ligand
        data_out.df['HETATM'] = data_out.df['HETATM'][:int(hetatm_record_len/2)]
        
    elif form_idx == 0:
        # Delete all 'HETATM' records
        data_out.df['HETATM'] = pd.DataFrame(columns=data_out.df['HETATM'].columns)

    return data_out

def save_structure(data, form_idx, output_dir):
    """ Save BioPandas object as a PDB record file.
    """
    output_filename = "{}.pdb".format(form_idx)

    data.to_pdb(path=os.path.join(output_dir, output_filename),
                    records=['ATOM', 'HETATM'],
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

    main_commandline()
