#!/usr/bin/env python
"""
This script provides useful funcs to all other scripts
"""
import yaml
import os

def read_config():
    # Read in config file
    with open("config.yaml") as yaml_file:
        # YAML loads a list of dictionaries
        config_list = yaml.full_load(yaml_file)
        # Convert list into dict
        config_dict = {key: value for dict in config_list for key, value in dict.items()}
    return config_dict
    