#!/usr/bin/env python3

#
# Copyright (C) 2020 Cambridge Astronomical Survey Unit
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <https://www.gnu.org/licenses/>.
#


import numpy as np
from astropy.table import Table
import os
import sys
import argparse
import logging

def test_path(*path):
    full_path = os.path.join(*path)
    if not os.path.exists(full_path):
        raise IOError('Sorry, but "{}" does not exist!'.format(full_path))
    return full_path


root_path = os.path.dirname(os.path.realpath(__file__))  # Path to this file
ifu_directory = test_path(root_path, '..', '..')  # Path to the 'ifu' directory, containing the whole workflow, etc.
workflow_directory = test_path(ifu_directory, 'workflow')  # within the `ifu_directory`

if ifu_directory not in sys.path:
    sys.path.append(ifu_directory)
    # print(sys.path)

from workflow.stage1 import create_ifu_driver_cat


def get_params_dict(filename,template_filename=None):
    """Read basic parameters and defaults from CSV file"""

    table = Table.read(filename, format='ascii.csv', comment='#')
    params = {}

    template = None
    if template_filename:
        template = Table.read(template_filename)
        template['PROGTEMP'] = template['PROGTEMP'].astype(str)

    
    allow_null = ['TACALLOC', 'IFU_PA_REQUEST'] # float values
    allow_blank = ['TARGPROG', 'cc_report']
    non_zero = ['TARGPRIO','HA_RANGE']
    
    for row in table:
        is_null = row['value'] in ['None','nan']
        try:
            floatval = float(row['value'])
            is_zero = floatval == 0
        except ValueError:
            is_zero = False

        try:
            is_blank = row['value'].mask
            is_blank = True
        except AttributeError:
            is_blank = False

        if (is_blank) and (not row['parameter'] in allow_blank):
             raise AssertionError("Quantity {} in params input CSV may not be defined as a blank-string value".format(row['parameter']))
            
        if (is_null) and (not row['parameter'] in allow_null):
            raise AssertionError("Quantity {} in params input CSV may not be defined as a null value".format(row['parameter']))
        if (is_zero) and (row['parameter'] in non_zero):
            raise AssertionError("Quantity {} in params input CSV may not be set to zero".format(row['parameter']))

        
        if (row['parameter'] == 'TACALLOC') and (row['value']) == 'None':
            params[row['parameter']] = np.str_("")
            continue

        if (row['parameter'] in allow_blank) and (is_blank):
            params[row['parameter']] = np.str_("")
            continue
        
        if row['value'] in ['None', 'nan']:
            # set as nan
            params[row['parameter']] = np.float32(np.nan)
            continue

        params[row['parameter']] = row['value']
        if (template is not None) and (row['parameter'] in template.columns):
            try:
                cast_val = row['value'].astype(template[row['parameter']].dtype)
            except ValueError:
                raise AssertionError("Quantity {} in params input CSV ('{}') could not be converted to datatype {}".format(row['parameter'],row['value'],str(template[row['parameter']].dtype)))
        
    return params


def get_data_dict(filename, params_dict, template_filename):
    """Read individual target information from CSV file, fill the rest with defaults"""


    targets = Table.read(filename, comment='#',converters={'PROGTEMP':str,'TARGPROG':str},format='ascii.csv')
    n_targets = len(targets)
    template = Table.read(template_filename)
    template['PROGTEMP'] = template['PROGTEMP'].astype(str)
    data_dict = {}
    
    for name in template.colnames:
        if name in targets.colnames:
            if hasattr(targets[name],'mask'):
                data = targets[name].astype(template[name].dtype)
                mask = targets[name].mask
                data_dict[name] = [_data if not _mask else None for _data, _mask in zip(data,mask)]                 
            else:
                data_dict[name] = targets[name].astype(template[name].dtype).tolist()
                
            if None in data_dict[name]:
                # replace instances of None with the default value (assuming it is in the params_dict)
                assert name in params_dict, "User specified a null value for '{}', but this quantity is not in the _params CSV file".format(name)
                data_dict[name] = [x if x is not None else params_dict[name] for x in data_dict[name]]
                # now cast to the required dtype
                data_dict[name] = list(np.array(data_dict[name],dtype=template[name].dtype))
                
            if name in params_dict:
                logging.info('Default value {0}={1} ignored'.format(name,params_dict[name]))
        else:
            if name in params_dict:
                data_dict[name] = n_targets * [params_dict[name].astype(template[name].dtype)]
                logging.info('Setting default value {0}={1}'.format(name,params_dict[name]))
            else:
                logging.error('ERROR: information about {} not found'.format(name))

    return data_dict


# ==============================================================================
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Create IFU driver catalogue from <base_name>_params.csv` and `<base_name>_target_list.csv`""")

    parser.add_argument('base_name',
                        help="""Root for all input and output files, conforming to the WEAVE standard:\n
                        '<TARGSRVY>_<trimester>' for Science Team catalogues\n
                        'W<S/V><trimester>-<xxx>' for SV and Open Time""")

    parser.add_argument('--input_directory', nargs='?', default=os.getcwd(),
                        help="if different from current directory")

    parser.add_argument('--output_directory',
                        help="if different from '<input_directory>/<base_name>'")

    parser.add_argument('--overwrite', dest='overwrite', action='store_true',
                        help='overwrite the output files')
    
    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')
    
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()))
    
    ############################################################################
    # Read input files to fill params_dict and data_dict

    input_directory = test_path(args.input_directory)
    ifu_driver_template = os.path.join(workflow_directory, 'aux', 'ifu_driver_template.fits')    
    logging.info('Reading {} input files from "{}":'.format(args.base_name,input_directory))
    param_file = os.path.join(input_directory, '{}_params.csv'.format(args.base_name))
    params_dict = get_params_dict(param_file, template_filename=ifu_driver_template)
    logging.info("Parameter file: {}".format(param_file))
    
    
    trimester = params_dict['trimester']
    author_name = params_dict['catauthor_name']
    author_email = params_dict['catauthor_email']
    tacalloc = params_dict['TACALLOC']
    report_verbosity = int(params_dict['report_verbosity'])
    cc_report = params_dict['cc_report']
    assert len(cc_report) <= 80, "Maximum permitted field size for cc_report field is 80 char"
        
    target_file = os.path.join(input_directory, f'{args.base_name}_target_list.csv')
    logging.info("Input target definition file: {}".format(target_file))
    
    logging.info("Base driver template file: {}".format(ifu_driver_template))

    data_dict = get_data_dict(target_file, params_dict, ifu_driver_template)

    if (trimester in ['2023A1','2022A2']) and (not 'HA_REQUEST' in data_dict):
        logging.warning("No HA_REQUEST in trimester {} CSV file - adding as exception".format(trimester))
        data_dict['HA_REQUEST'] = [0.0 for val in data_dict['TARGSRVY']]
    
    ############################################################################
    # Create the output directory if it does not exist

    if args.output_directory is None:
        output_directory = os.path.join(input_directory, f'{args.base_name}')
    else:
        output_directory = args.output_directory

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    output_filename = os.path.join(output_directory, f'{args.base_name}-ifu_driver_cat.fits')

    ############################################################################
    # Create the IFU driver catalogue
    
    overwrite = args.overwrite
    
    output = create_ifu_driver_cat(ifu_driver_template, data_dict, output_filename,
                                   trimester, tacalloc, author_name, author_email,
                                   report_verbosity=report_verbosity,
                                   overwrite=overwrite,cc_report=cc_report)
    logging.info("Driver FITS catalogue file output for {}: {}".format(args.base_name,output))
