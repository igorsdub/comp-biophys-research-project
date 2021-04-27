# -*- coding: utf-8 -*-
from typing import Counter
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv

import subprocess
import os
import glob
import src.utilities as utils
import numpy as np
import itertools
import pandas as pd


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
    pdb_filepaths = sorted(glob.glob(os.path.join(input_filepath, "1m9a.2.pdb")))
    # pdb_filepaths = sorted(glob.glob(os.path.join(input_filepath, "*.pdb")))

    # Simulate ENM
    for pdb_filepath in pdb_filepaths:
        cutoff_radius = find_smallest_cutoff_radius(pdb_filepath, output_filepath, start_cutoff_radius = 7.0)
        print("The smallest cutoff radius for 1M9A is {:.2f} angstroms".format(cutoff_radius))


def scan_pdb(pdb_filepath, output_filepath, flag_combination="-ca -het -c 8.00"):
    """ Executes Shell script with essential DDPT routines.
        For inputs see scan_pdb.sh
    """
    # Usage: scan_pdb.sh <pdb-filepath> <results-filepath> <cutoff>
    subprocess.call(['bash', 'src/simulation/scan_pdb.sh', pdb_filepath, output_filepath, flag_combination]) 

def find_smallest_cutoff_radius(pdb_filepath, output_filepath, start_cutoff_radius = 5.0, step=0.5):
    """ Finds minimal cutoff radius values for the ENM which 
        avoids floppy modes due-to underconnected EN.

        Starts with 5.0 angstrom cutoff radius and proceeds
        in 0.5 angstrom step up.
    """

    for cutoff_radius in np.arange(start_cutoff_radius, 15.5, 0.5):
        flag_combination = "-ca -c {}".format(cutoff_radius)
        subprocess.call(['bash', 'src/simulation/find_smallest_cutoff_radius.sh', pdb_filepath, output_filepath, flag_combination])

        eigenvalues = np.loadtxt(os.path.join(output_filepath, "eigenvalues"), max_rows=7)
        eigenvals_sum_6 = np.sum(eigenvalues[0:6])
        eigenvals_sum_7 = np.sum(eigenvalues[0:7])

        print("""
    Cutoff radius:          {:.2f}
    First 6 eigenvalue sum: {:.3e}
    First 7 eigenvalue sum: {:.3e}
        """.format(cutoff_radius, eigenvals_sum_6, eigenvals_sum_7))

        if (eigenvals_sum_6 < 1e-7) and (eigenvals_sum_7 > 1e-7):
            # Cutoff condition is met
            return cutoff_radius

        # # Create distance matrix
        # dist_data = np.loadtxt(os.path.join(output_filepath, "dist.dat"), \
        #     dtype={'names': ('res_num_i', 'res_num_j', 'dist'), 'formats': ('i4', 'i4', 'f4')})

        # counter = 0
        # while dist_data[counter][0] == dist_data[0][0]:
        #     counter += 1

        # total_res_num = counter + 1
        # dist_matrix = np.zeros((total_res_num, total_res_num))

        # # Populate matrix
        # for i in range(total_res_num):
        #     start_idx = i * total_res_num
        #     end_idx = (i+1) * total_res_num
        #     dist_matrix[i] = np.array([row[-1] for row in dist_data[start_idx:end_idx]])    # Use list comperhensiion

        # Next, create connectivity matrix by boolean condition
        # dist_matrix <= cutoff_radius
        # Substitute True with 1, and Flase with 0
        # Then, check connectivity of each residue
    return None
    

# def make_distance_matrix(path):


def brutally_scan_pdb(pdb_filepath, output_filepath, start_cutoff_radius=5.0):
    """ Performs simualtion for all possible GENENMM module useful flag
        combination. 
    """
    # DDPT flags in the ordr of apperas in GENENMM sourcecode
    mass_flag   = ['', '-mass']
    ca_flag     = ['-ca']   # Always present
    het_flag    = ['het']   # Always present
    lig1_flag   = ['', '-lig1']
    res_flag    = ['', '-res']
    # Cutoff radius
    cutoff_radius_values = np.arange(start_cutoff_radius, 15.5, 0.5)
    c_flag   = ["-c {}".format(cutoff_radius_value) for cutoff_radius_value in cutoff_radius_values] # Cutoff radius

    # Combine all flag lists
    all_flags = [mass_flag, ca_flag, het_flag, lig1_flag, res_flag, c_flag]

    # Create tuple with all flag permutaions (non-reapeating)
    all_flag_permutations = list(itertools.product(*all_flags))
    # Convert list of tuples into list lists
    all_flag_permutations = list(map(list, all_flag_permutations))

    for flag_combination in all_flag_permutations:

        # Usage: scan_pdb.sh <pdb-filepath> <results-filepath> <cutoff>
        subprocess.call(['bash', 'src/simulation/scan_pdb.sh', pdb_filepath, output_filepath, flag_combination]) 

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()
