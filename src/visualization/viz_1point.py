# Data manipulation
import numpy as np
import pandas as pd

# Data visualization
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os
# %matplotlib inline
import seaborn as sns
from pymol import cmd

def plot_heatmap(data, cbar_lbl='$K_{2}/K_{1}$', axis=None):
    '''Plots 1-point mutational scan dataframe. 
    >>> plot_heatmap(df_1, '$K_{2}/K_{1}$')
    '''

    xlbl = 'Amino Acid Number'
    ylbl = '$k_R/k$'
    
    if data.shape[0] == 25:
        ypositions = np.array([0, 4, 8, 12, 16, 20, 24]) + 0.5
        ylabels = ("0.25", "0.50", "0.75", "1.00", "2.00", "3.00", "4.00")
    else:
        ypositions = np.arange(data.shape[0]) + 0.5
        ylabels = [str(label) for label in data.index.values]
        
    if axis == None:
        sns.heatmap(data, center=data.loc[1.00].iloc[0], cmap = plt.cm.RdBu_r, xticklabels = 50, 
                     cbar_kws={'label': cbar_lbl})
        ax = plt.gca()
    else:
        sns.heatmap(data, center=data.loc[1.00].iloc[0], cmap = plt.cm.RdBu_r, xticklabels = 50, 
                     cbar_kws={'label': cbar_lbl}, ax=axis)
        ax = axis

    ax.invert_yaxis()
    ax.set_ylim(0, data.shape[0])
    
    ax.set_yticklabels(ax.get_yticks(), rotation = 0)
    ax.set(yticks=ypositions, yticklabels=ylabels)
    
    ax.set_xlim(0, data.shape[1])
    xlabels = np.arange(50, data.columns[-1] + 1, 50)
    # Shift x-ticks to cell's middle
    xpositions = xlabels + 0.5 - data.columns[0]
    ax.set(xticks=xpositions, xticklabels=xlabels)

    # Add frame around plot
    for _, spine in ax.spines.items():
        spine.set_visible(True)

    ax.set_xlabel(xlbl)
    ax.set_ylabel(ylbl)

    return ax

    # Set second axis
#     ax2 = ax.twiny()
#     # Active residues
#     # Move ticks to the center of a cell
#     res_act = np.concatenate((np.arange(41,42), np.arange(49,50), np.arange(143,146), np.arange(163,168), np.arange(187,193), np.arange(214,215), np.arange(284,287)))+0.5
#     # Shift residue positions
#     x = tuple(res_act-1)
#     # Set x-axis limit
#     ax2.set_xlim(0, 306-1)
#     # Set tick style
#     plt.style.use('seaborn')
#     forest = (0.2,0.6,0.2)
#     ax2.tick_params(axis="x", direction="out", length=5*f, width=f/2, color=forest)
#     # Hide tick labels
#     ax2.set_xticks(x)
#     ax2.set_xticklabels([])

def code_heatmap(data, vmin=-1, vmax=1, vcenter=0, cmap=plt.cm.RdBu_r):
    ''' Converts passed numerical value from coloured heatmap into hex colour code.
        Returns data frame of PyMOL hex colour value for each residue and colour for the center value.
        >>> convert_heatmap_to_colourcode(df_1c3b_wide.loc[4.00, :], center=df_1c3b_wide.loc[1.00, :][0])
        (colour_data, center_hex)
    '''
    
    norm = mpl.colors.TwoSlopeNorm(vmin=vmin, vmax=vmax, vcenter=vcenter)
    
    center_rgba = cmap(vcenter)
    center_hex = mpl.colors.to_hex(center_rgba).replace("#", "0x")
    
    if isinstance(data, np.ndarray):
        rgba_values = cmap(norm(data))
    else:
        rgba_values = cmap(norm(data.to_numpy()))
        
    hex_values = []
    
    # Convert RGB values into PyMOL hex format
    for i in range(np.shape(rgba_values)[0]):
            hex_values.append(mpl.colors.to_hex(rgba_values[i,:]).replace("#", "0x")) # PyMOL syntax
    
    if isinstance(data, pd.core.series.Series) or isinstance(data, pd.core.frame.DataFrame):
        colour_data = pd.DataFrame(hex_values, index=data.index, columns=[data.name]).transpose()
        colour_data.index.name = 'spring_strength'
    else:
        colour_data = pd.DataFrame(hex_values).transpose()
        
    return (colour_data, center_hex)

def colour_by_heatmap(colour_data, structure_path, molecule_name="protein", output_path="colour_by_heatmap", view=None):
    '''
DESCRIPTION

    Colours PDB structure by colour map data.
    output_filepath
    >>> colour_by_heatmap(df, structure_filepath="xxxx.pdb", mol_name="protein", out_dir='.')
    '''
    
    # if not (isinstance(colour_data, pd.core.series.Series) or isinstance(colour_data, dict)):
    #     print('''
    #             Passed data must be either dictionary or Pandas Series object.
    #             Key = residue number
    #             Value = PyMOL hex code
    #             ''')
    #     return None

    cmd.load(structure_path, object=molecule_name)
    
    # Set view
    if view == None:
        cmd.reset()
        cmd.orient()
    else:
        cmd.set_view(view)
    
    cmd.viewport(width=1200, height=1200)
    cmd.zoom(complete=1)

    cmd.set('cartoon_discrete_colors', 1)
    cmd.set('sphere_scale', 1)

    cmd.show_as('cartoon', molecule_name)
    cmd.color('white', molecule_name)

    # Iterate over the alpha-carbons
    # residue_numbers = []
    # cmd.iterate('{} and name CA'.format(molecule_name) , 'residue_numbers.append(resi)')

    # Colour the structure
    for residue_number in colour_data.columns:
        # print(colour_data[residue_number].item())
        cmd.color(colour_data[residue_number].item(), '{0} and resi {1}'.format(molecule_name, residue_number))

    png_out_path = output_path + ".png"
    pse_out_path = output_path + ".pse"

    cmd.save(pse_out_path)

    cmd.set('ray_opaque_background', 0)
    cmd.png(png_out_path, width=1200, height=1200, ray=1, quiet=0)
