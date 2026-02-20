#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:04:44 2024

@author: u6044586
"""

import os
import yaml
from planeflight_utils import _parse_geoschem_config, _get_tracer_name_num_mapping


def create_num_inputs(in_dir, out_dir, logfile_pth, config_yaml):
    # Parse the config file & get # of advectoed species:
    adv_species, simtype = _parse_geoschem_config(config_yaml)

    # Parse log file and get dict of tracer name/number pairs:
    name_map = _get_tracer_name_num_mapping(logfile_pth, len(adv_species), keys_are_nums=False)

    files = [f for f in os.listdir(in_dir) if os.path.isfile(os.path.join(in_dir, f))]

    for f in files:
        # Open file and read in all lines to a list.
        fin = open(os.path.join(in_dir, f), 'r')
        lines = fin.readlines()

        new_lns = []
        for line in lines:
            got_it = False
            for sp in list(name_map.keys()):
                if line.replace(' ', '').replace('\n', '') == sp and not got_it:
                    new_lns.append(name_map[sp] + '\n')
                    got_it = True
            if not got_it:
                new_lns.append(line)

        # Create new file with
        outF = open(os.path.join(out_dir, f), 'w')
        for line in new_lns:
            outF.write(line)

        del fin, lines, outF, new_lns

    return


def create_num_outputs(in_dir, out_dir, logfile_pth, config_yaml):

    def update_headers(header_line, name_map):
        # We are also checking for leading and trailing spaces
        headers = header_line.split()
        out_chars = list(header_line)
        for header in headers:
            if header in name_map:
                key = header
                value = name_map[key]
                key_index = header_line.find(key)

                while key_index != -1:
                    # Ensure that the key is matched as a whole word surrounded by spaces
                    is_full_match = (
                        (key_index == 0 or header_line[key_index - 1] == ' ') and
                        (key_index + len(key) == len(header_line) or header_line[key_index + len(key)] == ' ')
                    )

                    if is_full_match:
                        # Replace the found key with the corresponding value in out_chars
                        for i in range(len(key)):
                            out_chars[key_index + i] = ' '  # Clear original key

                        for i in range(len(value)):
                            if i < len(value):
                                out_chars[key_index + i] = value[i]  # Insert the new value

                    # Look for the next occurrence
                    key_index = header_line.find(key, key_index + 1)

        # Join the list back into a string
        out_str = ''.join(out_chars)

        return out_str

    # Parse the config file & get # of advected species:
    adv_species, simtype = _parse_geoschem_config(config_yaml)

    # Parse log file and get dict of tracer name/number pairs:
    name_map = _get_tracer_name_num_mapping(logfile_pth, len(adv_species), keys_are_nums=False)

    files = [f for f in os.listdir(in_dir) if os.path.isfile(os.path.join(in_dir, f))]

    for f in files:
        with open(os.path.join(in_dir, f), 'r') as infile, open(os.path.join(out_dir, f), 'w') as outfile:

            # Open file and read in all lines to a list.
            header_line1 = infile.readline()
            header_line2 = infile.readline()
            data_lines = infile.readlines()

            new_header1 = update_headers(header_line1, name_map)
            new_header2 = update_headers(header_line2, name_map)

            # Create new file with
            outfile.write(new_header1)
            outfile.write(new_header2)
            for line in data_lines:
                outfile.write(line)
    return


_run_dir = (
    '/uufs/chpc.utah.edu/common/home/haskins-group1/users/szhao'
    '/GEOS_CHEM/GC_RunDirs/gc_2x25_merra2_fullchem_base'
)
config_yaml = _run_dir + '/geoschem_config_jan_apr.yml'
spdb_path = _run_dir + '/species_database.yml'
logfile_pth = _run_dir + '/geoschem_2017_jan_apr.log'


_test_data = (
    '/uufs/chpc.utah.edu/common/home/u6044586/python_scripts'
    '/modules/gcpy_campaigns/test_plane_data'
)
i_names = _test_data + '/inputs_names/'
i_nums = _test_data + '/inputs_nums/'
o_names = _test_data + '/outputs_names/'
o_nums = _test_data + '/outputs_nums/'


# Create input fils with tracer #s from tracer name inputs:
# create_num_inputs(i_names, i_nums, logfile_pth, config_yaml)

# Create output files with tracer #s from tracer name outputs:
# create_num_outputs(o_names, o_nums, logfile_pth, config_yaml):

with open(spdb_path, 'r') as f:
    spdb = yaml.load(f, Loader=yaml.FullLoader)
