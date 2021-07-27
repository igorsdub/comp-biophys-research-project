#!/usr/bin/env python
""" This is the master script for recreating the results.

    It imports each of the key module fucntions 
    and runs them one by one.

    Run the whole thing from the root directory 
    to replicate all the results:
    
    $ python -m src.main
"""

import src.utilities as utils
import src.structure as struct
import src.simulation as sim
import src.data as data
import src.visualization as viz

config = utils.read_config()
# utils.clean()
struct.download_pdb.main(config['pdb']['rawFilePath'])
# mkpdb.main(config['pdb']['rawFilePath'])
# simenm.main(config['pdb']['proFilePath'], config['data']['rawFilePath'])
# prowt.main(config['data']['rawFilePath'], config['data']['proFilePath'])
# viz.main(config['data']['proFilePath'], config['data']['outPathScratch'])