# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv

import subprocess
import glob, os
from os.path import join as join_paths
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
@click.argument('input_dir', type=click.Path(exists=True))
@click.argument('output_dir', type=click.Path())
def main_commandline(input_dir, output_dir):
    """ Runs data visualization scripts to turn processed data (from data/processed)
        into plots (saved in scratch/).
    """
    logger = logging.getLogger(__name__)
    logger.info('making plots from processed data')

    main(input_dir, output_dir)
    
###################################################################################

def main(input_dir, output_dir):
    """ Runs data visualization scripts to turn processed data (from data/processed)
        into plots (saved in scratch/).
    """
    config = utils.read_config()
    pdb_codes = config['pdb']['codeList']

    mpl.rcParams.update(mpl.rcParamsDefault)  # VS Code plots not black
    plt.style.use(config['viz']['jupyter'])

    # Get paths
    # Directory path example: "data/processed/-c09.50/-mass-ca-het/0"
    cutoff_paths = sorted(glob.glob(join_paths(input_dir, "*"), recursive=True))

    # Plot parameters
    eigenvals_ylims = [-350, -100]
    diss_consts_ylims = [0, 15]
    coop_ylims = [0.9, 1.1]
    plot_no_modes = 300

    for cutoff_path in cutoff_paths:
        cutoff_flag = os.path.basename(cutoff_path)
        flag_paths = sorted(glob.glob(join_paths(cutoff_path, "*"), recursive=True))
        figure_path = "allo.{}.pdf".format(cutoff_flag)

        rows = 3
        cols = len(flag_paths)
        fig_side = 3 # inches
        fig_len = fig_side * cols
        fig_wid = fig_side * rows

        fig1, axs = plt.subplots(rows, cols, figsize=(fig_len, fig_wid))

        fig1.suptitle("{}".format(cutoff_flag))

        for idx, flag_path in enumerate(flag_paths):
            other_flags = os.path.basename(flag_path)

            # Load data
            eigenvals = pd.read_csv(join_paths(flag_path, "eigenvals.csv"), index_col='mode_number')
            diss_consts = pd.read_csv(join_paths(flag_path, "diss_consts.csv"), index_col='mode_number')
            coop = pd.read_csv(join_paths(flag_path, "coop.csv"), index_col='mode_number')


            # Plot data
            ttl = "{}".format(other_flags)
            ax1 = axs[0][idx]
            ax2 = axs[1][idx]
            ax3 = axs[2][idx]

            # EIGENVALS
            sns.scatterplot(data=eigenvals[7:plot_no_modes+7], ax=ax1)
            ax1.set_title(ttl, pad=15)
            ax1.set_xlabel("")
            ax1.set_xticklabels([])
            ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=False)

            # DISS_CONSTS
            sns.scatterplot(data=diss_consts[7:plot_no_modes+7], ax=ax2)
            ax2.set_xlabel("")
            ax2.set_xticklabels([])
            ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=False)

            # COOP
            sns.scatterplot(data=coop[7:plot_no_modes+7], ax=ax3)
            # Show non-cooperativity
            ax3.axhline(y=1.0, color='black', linestyle=':')
            ax3.get_legend().remove()
            ax3.set_ylim(coop_ylims)

            if idx == 0:
                ax1.set_ylabel("$\log(\lambda_{n})_{total}$")
                ax2.set_ylabel("$K$")
                ax3.set_xlabel("Mode number")
                ax3.set_ylabel("$K_{2}/K_{1}$")
            else:
                ax1.set_xlabel("")
                ax1.get_legend().remove()
                ax2.set_xlabel("")
                ax2.get_legend().remove()
                ax3.set_xlabel("")

            # Subplots' axes aspect ratio
            ax1.set_box_aspect(1)
            ax2.set_box_aspect(1)
            ax3.set_box_aspect(1)

        fig1.tight_layout(w_pad=1)
        
        plt.savefig(join_paths(output_dir, figure_path), bbox_inches='tight')
        plt.close(fig1)

    
    # Draw ENM
    # for pdb_code in pdb_codes:
    #     for form_idx in range(3):
    #         script_filepath = os.path.join(input_dir, "{}.{}.draw_enm.pml".format(pdb_code, form_idx))
    #         structure_filepath = os.path.join(input_dir, "{}.{}.CAonly.pdb".format(pdb_code, form_idx))
    #         cmd.reinitialize()
    #         structure_name ="{}.{}.enm".format(pdb_code, form_idx)
    #         draw_ENM(structure_filepath, script_filepath, structure_name=structure_name, \
    #             output_dir=output_dir, view=None)
    
    # # Draw eigenvectors
    # cmd.reinitialize()
    # for mode in range(7,31):
    #     first_structure = os.path.join(input_dir, "0.CAonly.pdb")
    #     last_structure = os.path.join(input_dir, "0.Mode_{:03}.pdb".format(mode))
    #     output_name = "mode_{}".format(mode)
    #     draw_eigenvecotrs(first_structure, last_structure, output_name)
    #     cmd.save(os.path.join(output_dir, "0.modes.pse"))

    # Plot 1-point mutational scan heatmap
    # for filename, data in allostery_1point_data.items():
    #     pdb_code = filename[:4]
    #     no_modes = int(filename[-3:].lstrip("0"))

    #     # Free energy
    #     for form_idx in range(3):
    #         column_name = "G_{}".format(form_idx)
    #         # Convert to wide format
    #         free_energy_1point_plot = data.pivot(index='spring_strength', columns='residue_number', values=column_name)

    #         _, ax = plt.subplots()
            
    #         ttl = "${}$ | {} | {} modes".format(column_name, pdb_code.upper(), no_modes)
    #         ax.set_title(ttl)
    #         viz_1point.plot_heatmap(free_energy_1point_plot, cbar_lbl="$G/(kT)$".format(column_name),axis=ax)

    #         plt.savefig(os.path.join(output_dir, "{}.1point.free_energy_{}.m{:03d}.png".format(pdb_code, form_idx ,no_modes)))

    #     plt.close('all')
    #     # Free energy change
    #     for form_idx in range(1,3):
    #         column_name = "dG_{}".format(form_idx)
    #         # Convert to wide format
    #         free_energy_1point_plot = data.pivot(index='spring_strength', columns='residue_number', values=column_name)

    #         _, ax = plt.subplots()
            
    #         ttl = "${}$ | {} | {} modes".format(column_name.replace("d", "\Delta "), pdb_code.upper(), no_modes)
    #         ax.set_title(ttl)
    #         viz_1point.plot_heatmap(free_energy_1point_plot, cbar_lbl="$\Delta G/(kT)$",axis=ax)

    #         plt.savefig(os.path.join(output_dir, "{}.1point.free_energy_change_{}.m{:03d}.png".format(pdb_code, form_idx, no_modes)))
        
    #     plt.close('all')
    #     # Cooperativity
    #     # Convert to wide format
    #     allostery_1point_plot = data.pivot(index='spring_strength', columns='residue_number', values='allostery')

    #     _, ax = plt.subplots()
        
    #     ttl = "1-point scan | {} | {} modes".format(pdb_code.upper(), no_modes)
    #     ax.set_title(ttl)
    #     viz_1point.plot_heatmap(allostery_1point_plot, axis=ax)

    #     plt.savefig(os.path.join(output_dir, "{}.1point.allostery.m{:03d}.png".format(pdb_code, no_modes)))

    # # Plot heatmap in real-space
    # filename = "1m9a.1point.allostery.m025"
    # data = allostery_1point_data[filename]
        
    # pdb_code = filename[:4]
    # no_modes = int(filename[-3:].lstrip("0"))

    # # Free energy for apo structure
    # for form_idx in range(1):
    #     column_name = "G_{}".format(form_idx)
    #     # Convert to wide format
    #     selected_data = data.pivot(index='spring_strength', columns='residue_number', values=column_name)
    #     for spring_strength in [0.25, 4.00]:
    #         selected_kcust_data= selected_data.loc[spring_strength, :]
    #         vmin = selected_data.min().min()
    #         vmax = selected_data.max().max()
    #         vcentre = selected_data.loc[1.00, :].iloc[0]
    #         # print("vmin = {}\nvcentre = {}\nvmax = {}".format(vmin, vcentre, vmax))

    #         colour_data, _ = viz_1point.code_heatmap(selected_kcust_data, vmin=vmin, vmax=vmax, vcenter=vcentre)
    #         path = os.path.join(output_dir, "{}.1point.free_energy.m{:03d}.k{:06.3f}".format(pdb_code, no_modes, spring_strength))
    #         cmd.delete('all')
    #         viz_1point.colour_by_heatmap(colour_data, structure_path="pdb/processed/1m9a.0.pdb", molecule_name="1m9a", output_path=path)

    # # Cooperativity
    # column_name = "allostery"
    # # Convert to wide format
    # selected_data = data.pivot(index='spring_strength', columns='residue_number', values=column_name)
    # for spring_strength in [0.25, 4.00]:
    #     selected_kcust_data= selected_data.loc[spring_strength, :]
    #     vmin = selected_data.min().min()
    #     vmax = selected_data.max().max()
    #     vcentre = selected_data.loc[1.00, :].iloc[0]
        
    #     colour_data, _ = viz_1point.code_heatmap(selected_kcust_data, vmin=vmin, vmax=vmax, vcenter=vcentre)

    #     path = os.path.join(output_dir, "{}.1point.allostery.m{:03d}.k{:06.3f}".format(pdb_code, no_modes, spring_strength))
    #     cmd.delete('all')
    #     viz_1point.colour_by_heatmap(colour_data, structure_path="pdb/processed/1m9a.2.pdb", molecule_name="1m9a", output_path=path)


def load_data(data_filepath):
    """ Load interim data into dataframe.
    """
    # Load data
    interim_data = pd.read_csv(data_filepath, sep=',', header=0)

    return interim_data

def draw_ENM(structure_filepath, script_filepath, structure_name="CAonly", output_dir='.', view=None):
    """ Draws elastic network model of a structure and saves image.
    """

    cmd.delete('all')
    cmd.load(structure_filepath, "CAonly")
    cmd.run(script_filepath)

    # Set name
    if structure_name == "CAonly":
        pass
    else:
        cmd.set_name("CAonly", structure_name)

    # Set view
    if view == None:
        cmd.orient()
    else:
        cmd.set_view(view)
    
    cmd.viewport(width=1200, height=1200)
    cmd.zoom(complete=1)

    png_filepath = os.path.join(output_dir, structure_name) + ".png"
    pse_filepath = os.path.join(output_dir, structure_name) + ".pse"

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

def save_figure(path, extensions=['pdf']):
    """ Saves a Matplolib figure in supplied extesnions.
    """
    for extension in extensions:
        plt.savefig("{}.{}".format(path, extension))

    return None

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()
