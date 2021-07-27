# -*- coding: utf-8 -*-
from posixpath import join
from typing import Counter
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv

import subprocess
import os
from os.path import join as join_paths
import glob
import src.utilities as utils
import numpy as np
import itertools
import pandas as pd
from biopandas.pdb import PandasPdb
from scipy.spatial.distance import pdist, squareform


@click.command()
@click.argument('input_dir', type=click.Path(exists=True))
@click.argument('output_dir', type=click.Path())
def main_comandline(input_dir, output_dir):
    """ Runs simualtion scripts for processed PDB data (from pdb/processed/) 
        to generate raw data ready to be processed (saved in data/raw/).
    """
    logger = logging.getLogger(__name__)
    logger.info('making simulation data set from processed PDB structures')

    # config = utils.read_config()
    # pdb_codes = config['pdb']['codeList']

    # Get PDB files in input_dir
    pdb_filepaths = sorted(glob.glob(join_paths(input_dir, "*.pdb")))

    apo_pdb_path = join_paths(input_dir, "0.pdb")
    
    dist = create_distance_matrix(apo_pdb_path)

    # Find the smallest cutoff for all EN beads
    # to have at least three springs in the ENM with 
    # isotropic cutoff radius
    cutoff_radius_3springs = np.amax(np.sort(dist, axis=0)[:, 2])

    # Find smallest non-floppy ENM cutoff radius
    # cutoff_radius_nonfloppy = find_smallest_cutoff_radius(apo_pdb_path, output_dir)

    # Simulate ENM
    for pdb_filepath in pdb_filepaths:
        run_enm(pdb_filepath, output_dir, \
            flag_combo="8.00")

def main(input_dir, output_dir):
    """ Runs simualtion scripts for processed PDB data (from pdb/processed/) 
        to generate raw data ready to be processed (saved in data/raw/).
    """
    logger = logging.getLogger(__name__)
    logger.info('making simulation data set from processed PDB structures')

    # config = utils.read_config()
    # pdb_codes = config['pdb']['codeList']

    # Get PDB files in input directory
    pdb_filepaths = sorted(glob.glob(join_paths(input_dir, "*.pdb")))

    apo_pdb_path = join_paths(input_dir, "0.pdb")

    dist = create_distance_matrix(apo_pdb_path)
    np.savetxt(join_paths(output_dir, "dist.csv"), dist, delimiter=',', fmt='%.3f')

    # Find the smallest cutoff for all EN beads
    # to have at least three springs in the ENM with 
    # isotropic cutoff radius
    sorted_dist = np.sort(dist, axis=0)
    cutoff_radius_3springs = np.amax(sorted_dist[2])
    print("Minimum 3 springs per EN bead\nCutoff radius = {:.3f}".\
        format(cutoff_radius_3springs))

    # Find smallest non-floppy ENM cutoff radius
    # cutoff_radius_nonfloppy = find_smallest_cutoff_radius(apo_pdb_path, output_dir)
    cutoff_radius_nonfloppy = 7.5
    
    # Brute-force ENM scan
    for pdb_filepath in pdb_filepaths:
        brute_force_scan(pdb_filepath, output_dir, start_cutoff_radius=cutoff_radius_nonfloppy)
    
    # Simulate ENM
    # for pdb_filepath in pdb_filepaths:
    #     run_enm(pdb_filepath, output_dir, \
    #         flag_combo="8.00")


def run_enm(pdb_filepath, output_dir, flag_combo="-ca -het -c 8.00"):
    """ Executes Shell script with essential DDPT routines.
        For inputs see run_enm.sh
    """
    # Usage: run_enm.sh <pdb-filepath> <results-filepath> <cutoff>
    subprocess.call(['bash', 'src/simulation/run_enm.sh', pdb_filepath, output_dir, flag_combo]) 

    return None

def create_distance_matrix(pdb_filepath):
    """ Creates distance matrix for PDB carbon alpha coordinates.
    """
    ppdb = PandasPdb().read_pdb(pdb_filepath)
    ca_records = ppdb.df['ATOM'][ppdb.df['ATOM']['atom_name'] == 'CA']
    ca_coord = ca_records[['x_coord', 'y_coord', 'z_coord']].to_numpy()

    # Calculate distance matrix
    dist = pdist(ca_coord)
    dist = squareform(dist)

    return dist

def find_smallest_cutoff_radius(pdb_filepath, output_dir, start_cutoff_radius = 5.0, step=0.5):
    """ Finds minimal cutoff radius values for the ENM which 
        avoids floppy modes due-to underconnected EN.
        Choose starting cutoff radius (default == 5.0) and scan up to 
        15.0 angstroms in 0.5 angstrom steps until non-floppy ENM is found.
    """

    for cutoff_radius in np.arange(start_cutoff_radius, 15.5, 0.5):
        flag_combo = "-c {} -ca".format(cutoff_radius)
        run_enm(pdb_filepath, output_dir, flag_combo=flag_combo)

        eigenvalues = np.loadtxt(join_paths(output_dir, "eigenvalues"), max_rows=7)
        eigenvals_sum_6 = np.sum(eigenvalues[0:6])
        eigenvals_sum_7 = np.sum(eigenvalues[0:7])

        print("""
    Cutoff radius:          {:.2f}
    First 6 eigenvalue sum: {:.3e}
    First 7 eigenvalue sum: {:.3e}
        """.format(cutoff_radius, eigenvals_sum_6, eigenvals_sum_7))

        if (eigenvals_sum_6 < 1e-7) and (eigenvals_sum_7 > 1e-7):
            # Cutoff condition is met
            print("The ENM has no floppy modes!")
            return cutoff_radius
        elif cutoff_radius < 15.0:
            print("The ENM is underconnected. Continue the scan.\n".format(cutoff_radius))
            continue
        else:
            print("15.0 angstrom cutoff radius is not enough\nCheck the PDB file".format(cutoff_radius))

    return None

def brute_force_scan(pdb_filepath, output_dir, start_cutoff_radius=5.0):
    """ Brute-force ENM scan to find an optimal ENM.
    """
    # DDPT flags in the ordr of apperas in GENENMM sourcecode
    mass_flag   = ['', '-mass']
    ca_flag     = ['-ca']   # Always present
    het_flag    = ['-het']   # Always present
    lig1_flag   = ['', '-lig1'] # Has no effect on apo form
    res_flag    = ['', '-res']

    cutoff_radii = np.arange(start_cutoff_radius, 15.5, 0.5)

    # Combine all flags without cutoff lists
    flags = [mass_flag, ca_flag, het_flag, lig1_flag, res_flag]

    # Create tuple with all flag permutaions (non-reapeating)
    flag_combos = list(itertools.product(*flags))
    # Convert list of tuples into list lists
    flag_combos = list(map(list, flag_combos))

    pdb_filename = os.path.splitext(os.path.basename(pdb_filepath))[0]

    # ANM (with cutoff radius)
    for cutoff_radius in cutoff_radii:
        for flag_combo in flag_combos:
              
            cutoff_flag = "-c {}".format(cutoff_radius)
            cutoff_flag_lbl = "-c{:05.2f}".format(cutoff_radius)

            output_subdir = join_paths(output_dir, cutoff_flag_lbl, \
                "".join(flag_combo).replace(" ", ""), pdb_filename)
            os.makedirs(output_subdir, exist_ok=True)

            appended_flag_combo = flag_combo.copy()
            appended_flag_combo.append(cutoff_flag)

            with open(join_paths(output_subdir, "main.log"), 'w') as log_file:
                # Usage: run_enm.sh <pdb-filepath> <results-filepath> <GENENMM-flags>
                subprocess.call(['bash', 'src/simulation/run_enm.sh', pdb_filepath, output_subdir, \
                    " ".join(appended_flag_combo)], stdout=log_file)
    
    # pfENM
    for flag_combo in flag_combos:
        pf_flag = "-pf"

        output_subdir = join_paths(output_dir, pf_flag, \
            "".join(flag_combo).replace(" ", ""), pdb_filename)
        os.makedirs(output_subdir, exist_ok=True)

        appended_flag_combo = flag_combo.copy()
        appended_flag_combo.append(pf_flag)

        with open(join_paths(output_subdir, "main.log"), 'w') as log_file:
            # Usage: run_enm.sh <pdb-filepath> <results-filepath> <GENENMM-flags>
            subprocess.call(['bash', 'src/simulation/run_enm.sh', pdb_filepath, output_subdir, \
                " ".join(appended_flag_combo)], stdout=log_file)


    return None

def write_cfile(input_data):
    """ Writes cfile that contains custom cutoff radii
        for different atom names, e.g. 'CA ', 'C  ', 'O  '.
    """
    "{:4s} {:7.3f}".format(atom_name, cutoff_radius)
    

def write_ffile(input_data):
    """ Writes ffile that contains custom residue-residue
        intercations.
    """

    " {:4d} {:1s} {:4d} {:1s} {:8.3f}".format(res_1, chain_1, res_2, chain_2, k_cust)

def write_spfile(input_data):
    """ Writes spfile that contains custom residue-residue
        springs regardless of the global cutoff radius.
    """

    " {:4d} {:1s} {:4d} {:1s}".format(res_1, chain_1, res_2, chain_2)

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()
