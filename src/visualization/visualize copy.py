# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv

import subprocess
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import numpy as np
from pymol import cmd
import src.utilities as utils
from src.visualization.modevectors import modevectors
import src.visualization.viz_1point as viz_1point


@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())
def main(input_filepath, output_filepath):
    """ Runs data visualization scripts to turn processed data (from data/processed)
        into plots (saved in scratch/).
    """
    logger = logging.getLogger(__name__)
    logger.info('making plots from processed data')

    config = utils.read_config()
    pdb_codes = config['pdb']['codeList']

    mpl.rcParams.update(mpl.rcParamsDefault)  # VS Code plots not black
    plt.style.use(config['viz'])
    
    # Get filepaths
    entropy_wt_filepaths = sorted(glob.glob(os.path.join(input_filepath, "*.entropy")))
    allostery_wt_filepaths = sorted(glob.glob(os.path.join(input_filepath, "*.allostery")))
    allostery_1point_filepaths = sorted(glob.glob(os.path.join(input_filepath, "*.1point.allostery.m???.csv")))

    # Load data
    entropy_wt_data = {os.path.basename(filepath).replace(".entropy", "") : load_data(filepath) for filepath in entropy_wt_filepaths}
    allostery_wt_data = {os.path.basename(filepath).replace(".allostery", "") : load_data(filepath) for filepath in allostery_wt_filepaths}
    allostery_1point_data = {os.path.basename(filepath).replace(".csv", "") : load_data(filepath) for filepath in allostery_1point_filepaths}

    # Plot wild-type data
    # for pdb_code in pdb_codes:
    #     entropy_wt_plot = entropy_wt_data[pdb_code][:100]
    #     allostery_wt_plot = allostery_wt_data[pdb_code][:100]

    #     # Entropy
    #     _, ax = plt.subplots()

    #     x_lbl = "Mode (non-trivial)"
    #     y_lbl = "$-S/(kT)$"
    #     ttl = "Entropy | {}".format(pdb_code.upper())

    #     ax.set_xlabel(x_lbl)
    #     ax.set_ylabel(y_lbl)
    #     ax.set_title(ttl)

    #     for column_name in ['S_0', 'S_1', 'S_2']: 
    #         sns.lineplot(data=entropy_wt_plot[column_name], label="${}$".format(column_name), ax=ax)

    #     ax.legend()

    #     plt.savefig(os.path.join(output_filepath, "{}.wt.entropy.png".format(pdb_code)))

    #     # Free energy
    #     _, ax = plt.subplots()

    #     x_lbl = "Mode (non-trivial)"
    #     y_lbl = "G/(kT)"
    #     ttl = "Free energy | {}".format(pdb_code.upper())

    #     ax.set_xlabel(x_lbl)
    #     ax.set_ylabel(y_lbl)
    #     ax.set_title(ttl)

    #     for column_name in ['G_0', 'G_1', 'G_2']: 
    #         sns.lineplot(data=allostery_wt_plot[column_name], label="${}$".format(column_name), ax=ax)

    #     ax.legend()

    #     plt.savefig(os.path.join(output_filepath, "{}.wt.free_energy.png".format(pdb_code)))

    #     # Free energy change
    #     _, ax = plt.subplots()

    #     x_lbl = "Mode (non-trivial)"
    #     y_lbl = r"$\Delta G/(kT)$"
    #     ttl = "Free energy change | {}".format(pdb_code.upper())

    #     ax.set_xlabel(x_lbl)
    #     ax.set_ylabel(y_lbl)
    #     ax.set_title(ttl)

    #     for column_name in ['dG_1', 'dG_2']: 
    #         sns.lineplot(data=allostery_wt_plot[column_name], label="${}$".format(column_name), ax=ax)

    #     ax.legend()

    #     plt.savefig(os.path.join(output_filepath, "{}.wt.free_energy_change.png".format(pdb_code)))

    #     # Allostery
    #     _, ax1 = plt.subplots(figsize=(12,6))

    #     x_lbl = "Mode (non-trivial)"
    #     y1_lbl = "$K_{2}/K_{1}$"
    #     ttl = "Allostery | {}".format(pdb_code.upper())

    #     ax1.set_xlabel(x_lbl)
    #     ax1.set_ylabel(y1_lbl)
    #     ax1.set_title(ttl)

    #     # Show non-cooperative region
    #     ax1.axhline(y=1.0, color='black', linestyle=':', linewidth=1)

    #     sns.lineplot(data=allostery_wt_plot['allostery'], ax=ax1)
        
    #     # Plot ddG
    #     ax2 = ax1.twinx()
    #     y2_lbl = r"$\Delta \Delta G/(kT)$"
    #     ax2.set_ylabel(y2_lbl)

    #     sns.lineplot(data=allostery_wt_plot['ddG'], ax=ax2)

    #     plt.savefig(os.path.join(output_filepath, "{}.wt.allostery.png".format(pdb_code)))

    # plt.close('all')

    # Plot 1-point mut scan heatmap
    for filename, data in allostery_1point_data.items():
        pdb_code = filename[:4]
        no_modes = int(filename[-3:].lstrip("0"))

        # Free energy
        for form_idx in range(3):
            column_name = "G_{}".format(form_idx)
            # Convert to wide format
            free_energy_1point_plot = data.pivot(index='spring_strength', columns='residue_number', values=column_name)

            _, ax = plt.subplots()
            
            ttl = "${}$ | {} | {} modes".format(column_name, pdb_code.upper(), no_modes)
            ax.set_title(ttl)
            viz_1point.plot_heatmap(free_energy_1point_plot, cbar_lbl="$G/(kT)$".format(column_name),axis=ax)

            plt.savefig(os.path.join(output_filepath, "{}.1point.free_energy_{}.m{:03d}.png".format(pdb_code, form_idx ,no_modes)))

        plt.close('all')
        # Free energy change
        for form_idx in range(1,3):
            column_name = "dG_{}".format(form_idx)
            # Convert to wide format
            free_energy_1point_plot = data.pivot(index='spring_strength', columns='residue_number', values=column_name)

            _, ax = plt.subplots()
            
            ttl = "${}$ | {} | {} modes".format(column_name.replace("d", "\Delta "), pdb_code.upper(), no_modes)
            ax.set_title(ttl)
            viz_1point.plot_heatmap(free_energy_1point_plot, cbar_lbl="$\Delta G/(kT)$",axis=ax)

            plt.savefig(os.path.join(output_filepath, "{}.1point.free_energy_change_{}.m{:03d}.png".format(pdb_code, form_idx, no_modes)))
        
        plt.close('all')
        # Cooperativity
        # Convert to wide format
        allostery_1point_plot = data.pivot(index='spring_strength', columns='residue_number', values='allostery')

        _, ax = plt.subplots()
        
        ttl = "1-point scan | {} | {} modes".format(pdb_code.upper(), no_modes)
        ax.set_title(ttl)
        viz_1point.plot_heatmap(allostery_1point_plot, axis=ax)

        plt.savefig(os.path.join(output_filepath, "{}.1point.allostery.m{:03d}.png".format(pdb_code, no_modes)))

    # Plot heatmap in real-space
    filename = "1m9a.1point.allostery.m025"
    data = allostery_1point_data[filename]
        
    pdb_code = filename[:4]
    no_modes = int(filename[-3:].lstrip("0"))

    # Free energy for apo structure
    for form_idx in range(1):
        column_name = "G_{}".format(form_idx)
        # Convert to wide format
        selected_data = data.pivot(index='spring_strength', columns='residue_number', values=column_name)
        for spring_strength in [0.25, 4.00]:
            selected_kcust_data= selected_data.loc[spring_strength, :]
            vmin = selected_data.min().min()
            vmax = selected_data.max().max()
            vcentre = selected_data.loc[1.00, :].iloc[0]
            # print("vmin = {}\nvcentre = {}\nvmax = {}".format(vmin, vcentre, vmax))

            colour_data, _ = viz_1point.code_heatmap(selected_kcust_data, vmin=vmin, vmax=vmax, vcenter=vcentre)
            path = os.path.join(output_filepath, "{}.1point.free_energy.m{:03d}.k{:06.3f}".format(pdb_code, no_modes, spring_strength))
            cmd.delete('all')
            viz_1point.colour_by_heatmap(colour_data, structure_path="pdb/processed/1m9a.0.pdb", molecule_name="1m9a", output_path=path)

    # Cooperativity
    column_name = "allostery"
    # Convert to wide format
    selected_data = data.pivot(index='spring_strength', columns='residue_number', values=column_name)
    for spring_strength in [0.25, 4.00]:
        selected_kcust_data= selected_data.loc[spring_strength, :]
        vmin = selected_data.min().min()
        vmax = selected_data.max().max()
        vcentre = selected_data.loc[1.00, :].iloc[0]
        
        colour_data, _ = viz_1point.code_heatmap(selected_kcust_data, vmin=vmin, vmax=vmax, vcenter=vcentre)

        path = os.path.join(output_filepath, "{}.1point.allostery.m{:03d}.k{:06.3f}".format(pdb_code, no_modes, spring_strength))
        cmd.delete('all')
        viz_1point.colour_by_heatmap(colour_data, structure_path="pdb/processed/1m9a.2.pdb", molecule_name="1m9a", output_path=path)


    # Draw ENM
    # for pdb_code in pdb_codes:
    #     for form_idx in range(3):
    #         script_filepath = os.path.join(input_filepath, "{}.{}.draw_enm.pml".format(pdb_code, form_idx))
    #         structure_filepath = os.path.join(input_filepath, "{}.{}.CAonly.pdb".format(pdb_code, form_idx))
    #         cmd.reinitialize()
    #         structure_name ="{}.{}.enm".format(pdb_code, form_idx)
    #         draw_ENM(structure_filepath, script_filepath, structure_name=structure_name, \
    #             output_filepath=output_filepath, view=None)
    
    # # Draw eigenvectors
    # cmd.reinitialize()
    # for mode in range(7,31):
    #     first_structure = os.path.join(input_filepath, "0.CAonly.pdb")
    #     last_structure = os.path.join(input_filepath, "0.Mode_{:03}.pdb".format(mode))
    #     output_name = "mode_{}".format(mode)
    #     draw_eigenvecotrs(first_structure, last_structure, output_name)
    #     cmd.save(os.path.join(output_filepath, "0.modes.pse"))

def load_data(data_filepath):
    """ Load interim data into dataframe.
    """
    # Load data
    interim_data = pd.read_csv(data_filepath, sep=',', header=0)

    return interim_data

def draw_ENM(structure_filepath, script_filepath, structure_name="enm", output_filepath='.', view=None):
    """ Draws elastic network model of a structure and save image.
    """

    cmd.delete('all')
    cmd.load(structure_filepath, "enm")
    cmd.run(script_filepath)

    # Set name
    if structure_name == "enm":
        pass
    else:
        cmd.set_name("enm", "enm")

    # Set view
    if view == None:
        cmd.orient()
    else:
        cmd.set_view(view)
    
    cmd.viewport(width=1200, height=1200)
    cmd.zoom(complete=1)

    png_filepath = os.path.join(output_filepath, structure_name) + ".png"
    pse_filepath = os.path.join(output_filepath, structure_name) + ".pse"

    cmd.save(pse_filepath)
    cmd.set('ray_opaque_background', 0)
    cmd.png(png_filepath, width=1200, height=1200, ray=1)

    return (pse_filepath, png_filepath)

def draw_eigenvecotrs(first_structure, last_structure, output_name):
    """ Draws eginvectors based on two PDB structure coordinate difference.
    """
    cmd.delete("first_obj_frame")
    cmd.delete("last_obj_frame")
    cmd.delete(output_name)
    cmd.load(first_structure, object="first_obj_frame")
    cmd.load(last_structure, object="last_obj_frame")
    modevectors("first_obj_frame", "last_obj_frame", outname=output_name,cutoff=0, cut=0, factor=5)

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()
