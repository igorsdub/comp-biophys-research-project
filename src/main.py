#!/usr/bin/env python
""" This is the master script for recreating the results

    It imports each of the key other scripts and
    runs them one by one.

    Run the whole thing from the root directory 
    to replicate all of the python analysis
    $ python -m src.main
"""

# import src.pdb.download_pdb as dl_pdb
# import src.pdb.make_interim_pdb as mk_int_pdb
# import src.pdb.process_pdb as pro_pdb
# import src.simulation.simulate_enm as sim_enm
# import src.data.make_interim_dataset as mk_int_dataset
# import src.data.make_dataset as pro_data
# import src.visualization.visualize as viz

# dl_pdb.main(config['pdb']['rawFilePath'])
# mk_int_pdb.main(config['pdb']['rawFilePath'], config['pdb']['intFilePath'])
# pro_pdb.main(config['pdb']['intFilePath'], config['pdb']['proFilePath'])
# sim_enm.main(config['pdb']['proFilePath'], config['data']['rawFilePath'])
# mk_int_dataset.main(config['data']['rawFilePath'], config['data']['intFilePath'])
# pro_data.main(config['data']['intFilePath'], config['data']['proFilePath'])
# viz.main(config['data']['proFilePath'], config['data']['outPathScratch'])

import subprocess
import src.utilities as utils

config = utils.read_config()

subprocess.call(['python', '-m', 'src.pdb.download_pdb', config['pdb']['rawFilePath']])
subprocess.call(['python', '-m', 'src.pdb.make_interim_pdb', config['pdb']['rawFilePath'], config['pdb']['intFilePath']])
subprocess.call(['python', '-m', 'src.pdb.process_pdb', config['pdb']['intFilePath'], config['pdb']['proFilePath']])
subprocess.call(['python', '-m', 'src.simulation.simulate_enm', config['pdb']['proFilePath'], config['data']['rawFilePath']])
subprocess.call(['python', '-m', 'src.data.transform_raw_wt', config['data']['rawFilePath'], config['data']['intFilePath']])
subprocess.call(['python', '-m', 'src.data.process_interim_wt', config['data']['intFilePath'], config['data']['proFilePath']])
subprocess.call(['python', '-m', 'src.visualization.visualize', config['data']['proFilePath'], config['data']['outPathScratch']])


