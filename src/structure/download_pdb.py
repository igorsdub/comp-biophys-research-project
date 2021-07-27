# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv

import os
import src.utilities as utils
from urllib.request import urlretrieve

@click.command()
@click.argument('output_dir', type=click.Path())
def main_commandline(output_dir):
    """ Downloads raw PDB structures listed in YAML file
        from PDB website (saved in pdb/raw).
    """
    logger = logging.getLogger(__name__)
    logger.info('download raw PDB files from PDB website')

    main(output_dir)
    

def main(output_dir):
    """ Downloads raw PDB structures listed in YAML file
        from PDB website (saved in pdb/raw).
    """
    config = utils.read_config()
    pdb_code = config['pdb']['id']
    download_pdb(pdb_code, output_dir, biounit = True, compressed = False)

    return None


def download_pdb(pdb_code, output_dir, biounit = True, compressed = False):
    """ Downloads raw PDB files form a list of PDB IDs.
        Authored by Chris Swain (http://www.macinchem.org)
        Modified by Igors Dubanevics (https://github.com/igordub)
        Copyright CC-BY
    """
    # Add .pdb extension and remove ':1' suffix in entities
    filename = "{:4s}.pdb".format(pdb_code[:4])
    
    # Add '1' if biounit
    if biounit:
        filename = "{}1".format(filename)
    # Add .gz extenison if compressed
    elif compressed:
        filename = "{}.gz".format(filename)
    
    url = os.path.join("https://files.rcsb.org/download/", filename.lower())
    destination_file = os.path.join(output_dir, filename)
    # Download file
    urlretrieve(url, destination_file)

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
