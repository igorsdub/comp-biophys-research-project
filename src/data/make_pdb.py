# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv

import os
import glob
from biopandas.pdb import PandasPdb

@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())

def main(input_filepath, output_filepath):
    """ Runs data processing scripts to turn raw data from (../raw) into
        cleaned data ready to be analyzed (saved in ../processed).
    """

    pdb_capsid_filename = "Polyomaviridae_Large_72N_nogroups.pdb"

    pdb_capsid = PandasPdb().read_pdb(os.path.join(input_filepath, pdb_capsid_filename))

    pdb_capsid.df['ATOM'] = pdb_capsid.df['HETATM']

    pdb_capsid.df['ATOM']['record_name'].replace(to_replace="HETATM", value="ATOM", inplace=True)
    pdb_capsid.df['ATOM']['atom_name'].replace(to_replace="C", value="CA", inplace=True)
    pdb_capsid.df['ATOM']['chain_id'].replace(to_replace="", value="A", inplace=True)
    pdb_capsid.df['ATOM']['residue_number'] = pdb_capsid.df['ATOM']['atom_number']

    logger = logging.getLogger(__name__)
    logger.info('making final PDB file from PDB raw file')

    pdb_capsid.to_pdb(path=os.path.join(output_filepath, "viral_capsid.72n.pdb"), 
            records=['ATOM'], 
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
