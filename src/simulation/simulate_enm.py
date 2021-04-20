# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv

import subprocess
import os
import glob
import src.utilities as utils

@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())
def main(input_filepath, output_filepath):
    """ Runs simualtion scripts for processed PDB data (from pdb/processed/) 
        to generate raw data ready to be processed (saved in data/raw/).
    """
    logger = logging.getLogger(__name__)
    logger.info('making simulation data set from processed PDB structures')

    # config = utils.read_config()
    # pdb_codes = config['pdb']['codeList']

    # Get PDB files in input_filepath
    pdb_filepaths = sorted(glob.glob(os.path.join(input_filepath, "*.pdb")))

    # Simulate
    for pdb_filepath in pdb_filepaths:
        scan_pdb(pdb_filepath, output_filepath, cutoff=8.00)


def scan_pdb(pdb_filepath, output_filepath, cutoff=8.00):
    """ Executes Shell script with essential DDPT routines.
        For inputs see scan_pdb.sh
    """
    # Usage: scan_pdb.sh <pdb-filepath> <results-filepath> <cutoff>
    subprocess.call(['bash', 'src/simulation/scan_pdb.sh', pdb_filepath, output_filepath, str(cutoff)]) 

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()
