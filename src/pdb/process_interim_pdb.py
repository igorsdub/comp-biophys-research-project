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
    """ Proccesses interim PDB structures (from pdb/interim/) and creates PDB 
        structural forms for simulations (saved in pdb/processed/).
    """
    logger = logging.getLogger(__name__)
    logger.info('making processed PDB forms from interim PDB structures')

    config = utils.read_config()
    pdb_codes = config['pdb']['codeList']

    # Data import
    pdb_interim_structures = {pdb_code : load_structure(pdb_code, input_filepath, file_extension="pdb") \
                        for pdb_code in pdb_codes}

    # Data processing
    pdb_processed_structures = {}
    for pdb_code, PandasPDB_object in pdb_interim_structures.items():
        # Delete the first residue (1) and the last residues (216) from both chains
        PandasPDB_object.df['ATOM'] = PandasPDB_object.df['ATOM'][(PandasPDB_object.df['ATOM']['residue_number'] != 1) \
            & (PandasPDB_object.df['ATOM']['residue_number'] != 216)]
        # Create structural forms
        pdb_forms = {form_idx : create_form(PandasPDB_object, form_idx=form_idx) for form_idx in range(3)}

        pdb_processed_structures[pdb_code] = pdb_forms

    # Save processed data
    for pdb_code, pdb_forms in pdb_processed_structures.items():
        for form_idx, pdb_form in pdb_forms.items():
            save_structure(pdb_form, pdb_code, output_filepath, form_idx=form_idx, )


def load_structure(pdb_code, input_filepath, file_extension="pdb"):
    """ Loads PDB file inot BioPandas object.
    """
    pdb_filepath = os.path.join(input_filepath, "{}.{}".format(pdb_code, file_extension))
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
