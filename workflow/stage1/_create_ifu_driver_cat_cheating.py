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


import glob
import os

import numpy as np

from workflow.stage1 import create_ifu_driver_cat
from workflow.utils.get_data_from_xmls import \
    get_spa_data_of_target_fibres_from_xmls, get_trimester_from_xmls, \
    get_author_from_xmls, get_cc_report_from_xmls, \
    get_report_verbosity_from_xmls
from workflow.utils.get_progtemp_info import get_obsmode_from_progtemp


def get_keywords_info(xml_filename_list):

    trimester = get_trimester_from_xmls(xml_filename_list)
    author = get_author_from_xmls(xml_filename_list)
    report_verbosity = get_report_verbosity_from_xmls(xml_filename_list)
    cc_report = get_cc_report_from_xmls(xml_filename_list)
    
    return trimester, author, report_verbosity, cc_report


def get_data_dict(xml_filename_list):
    
    data_dict = get_spa_data_of_target_fibres_from_xmls(xml_filename_list,
                                                        post_configure=False)
    
    # Add a column with the IFU_PA_REQUEST values. It will be 0 for LIFU and
    # NULL for mIFU
    
    data_dict['IFU_PA_REQUEST'] = [
        0. if get_obsmode_from_progtemp(progtemp) == 'LIFU'
           else np.nan
        for progtemp in data_dict['PROGTEMP']
    ]
    
    return data_dict

    
if __name__ == '__main__':

    ############################################################################
    # Set the the list of XML files used for cheating
    
    xml_files_pattern = os.path.join('..', 'stage4', 'input', '*.xml')
    
    xml_filename_list = glob.glob(xml_files_pattern)
    xml_filename_list.sort()
    
    ############################################################################
    # Set the location of the template and the output file and directory
    
    ifu_driver_template = os.path.join('aux', 'ifu_driver_template.fits')
    
    output_dir = 'output'
    output_filename = os.path.join(output_dir,
                                   'WC_2020A1-ifu_driver_cat-cheating.fits')
    
    ############################################################################
    # Create the output directory if it does not exist
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    ############################################################################
    # Get the dictionary with the data which will populate the template from the
    # XMLs
    
    data_dict = get_data_dict(xml_filename_list)
    
    ############################################################################
    # Set the needed information to populate some keywords of the primary header

    trimester, author, report_verbosity, cc_report = \
        get_keywords_info(xml_filename_list)
    
    ############################################################################
    # Create the IFU driver catalogue
    
    create_ifu_driver_cat(ifu_driver_template, data_dict, output_filename,
                          trimester, author, report_verbosity=report_verbosity,
                          cc_report=cc_report)

