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


import os

import numpy as np

from workflow.stage1 import create_ifu_driver_cat


def set_keywords_info():

    trimester = '2020A1'
    author = 'a@domain.com'
    report_verbosity = 1
    cc_report = 'b@domain.com,c@domain.com'
    
    return trimester, author, report_verbosity, cc_report


def get_data_dict():
    
    data_dict = {}
    
    data_dict['TARGSRVY'] = \
        ['WC', 'WC', 'WC', 'WC', 'WC', 'WC', 'WC', 'WC', 'WC', 'WC',
         'WC', 'WC', 'WC', 'WC', 'WC', 'WC', 'WC', 'WC', 'WC', 'WC',
         'WC', 'WC', 'WC', 'WC', 'WC', 'WC', 'WC', 'WC', 'WC', 'WC',
         'WC', 'WC', 'WC', 'WC', 'WC', 'WC', 'WC', 'WC', 'WC', 'WC',
         'WC', 'WC', 'WC', 'WC', 'WC', 'WC', 'WC', 'WC']
    
    data_dict['TARGPROG'] = \
        ['LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR',
         'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR',
         'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR',
         'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR',
         'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR', 'LR']
    
    data_dict['TARGID'] = \
        ['PSZ1_125A', 'PSZ1_125A', 'PSZ1_125A', 'PSZ1_125A', 'PGC_91697',
           'PGC_65950', 'PGC_1051614', 'PGC_1048686', 'PGC_1058689', 'PGC_1044404',
         'PGC_1045221', 'PGC_1054814', 'PGC_1041786', 'PGC_1048019', 'PGC_1065694',
           'PGC_1051735', 'PGC_1053791', 'PGC_1060024', 'PGC_1059931', 'PGC_1063030',
         'PGC_1053595', 'PGC_1055624', 'PGC_1049654', 'PGC_1060631', 'PGC_1058049',
           'PGC_1056295', 'PGC_91697', 'PGC_65950', 'PGC_1051614', 'PGC_1048686',
         'PGC_1058689', 'PGC_1044404', 'PGC_1045221', 'PGC_1054814', 'PGC_1041786',
           'PGC_1048019', 'PGC_1065694', 'PGC_1051735', 'PGC_1053791', 'PGC_1060024',
         'PGC_1059931', 'PGC_1063030', 'PGC_1053595', 'PGC_1055624', 'PGC_1049654',
           'PGC_1060631', 'PGC_1058049', 'PGC_1056295']
    
    data_dict['TARGNAME'] = \
        ['PSZ1_125A', 'PSZ1_125A', 'PSZ1_125A', 'PSZ1_125A', 'CLUSTER1',
           'CLUSTER1', 'CLUSTER1', 'CLUSTER1', 'CLUSTER1', 'CLUSTER1',
         'CLUSTER1', 'CLUSTER1', 'CLUSTER1', 'CLUSTER1', 'CLUSTER1',
           'CLUSTER1', 'CLUSTER1', 'CLUSTER1', 'CLUSTER1', 'CLUSTER1',
         'CLUSTER1', 'CLUSTER1', 'CLUSTER1', 'CLUSTER1', 'CLUSTER1',
           'CLUSTER1', 'CLUSTER1', 'CLUSTER1', 'CLUSTER1', 'CLUSTER1',
         'CLUSTER1', 'CLUSTER1', 'CLUSTER1', 'CLUSTER1', 'CLUSTER1',
           'CLUSTER1', 'CLUSTER1', 'CLUSTER1', 'CLUSTER1', 'CLUSTER1',
         'CLUSTER1', 'CLUSTER1', 'CLUSTER1', 'CLUSTER1', 'CLUSTER1',
           'CLUSTER1', 'CLUSTER1', 'CLUSTER1']
    
    data_dict['TARGPRIO'] = \
        [10.0, 10.0, 10.0, 10.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
         1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
         1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
         1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
         1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    
    data_dict['PROGTEMP'] = \
        ['41331', '41331', '41331', '41331', '71331',
           '71331', '71331', '71331', '71331', '71331',
         '71331', '71331', '71331', '71331', '71331',
           '71331', '71331', '71331', '71331', '71331',
         '71331', '71331', '71331', '71331', '71331',
           '71331', '71661', '71661', '71661', '71661',
         '71661', '71661', '71661', '71661', '71661',
           '71661', '71661', '71661', '71661', '71661',
         '71661', '71661', '71661', '71661', '71661',
           '71661', '71661', '71661']
    
    data_dict['OBSTEMP'] = \
        ['IAEEB', 'IAEEB', 'IAEEB', 'IAEEB', 'IAEEB',
           'IAEEB', 'IAEEB', 'IAEEB', 'IAEEB', 'IAEEB',
         'IAEEB', 'IAEEB', 'IAEEB', 'IAEEB', 'IAEEB',
           'IAEEB', 'IAEEB', 'IAEEB', 'IAEEB', 'IAEEB',
         'IAEEB', 'IAEEB', 'IAEEB', 'IAEEB', 'IAEEB',
           'IAEEB', 'IAEEB', 'IAEEB', 'IAEEB', 'IAEEB',
         'IAEEB', 'IAEEB', 'IAEEB', 'IAEEB', 'IAEEB',
           'IAEEB', 'IAEEB', 'IAEEB', 'IAEEB', 'IAEEB',
         'IAEEB', 'IAEEB', 'IAEEB', 'IAEEB', 'IAEEB',
           'IAEEB', 'IAEEB', 'IAEEB']
    
    data_dict['GAIA_RA'] = \
        [316.369609537, 316.369609537, 316.369058752, 316.370143237, 316.9558,
           315.5608, 315.9079, 315.8575, 315.67, 316.6396,
         316.2858, 316.7646, 316.3958, 316.705, 316.5483,
           316.7125, 316.4896, 316.1362, 316.9229, 315.6946,
         316.9179, 316.0612, 316.7354, 316.2367, 316.5967,
           316.9825, 316.9558, 315.5608, 315.9079, 315.8575,
         315.67, 316.6396, 316.2858, 316.7646, 316.3958,
           316.705, 316.5483, 316.7125, 316.4896, 316.1362,
         316.9229, 315.6946, 316.9179, 316.0612, 316.7354,
           316.2367, 316.5967, 316.9825]
    
    data_dict['GAIA_DEC'] = \
        [-4.71060356792, -4.71060356792, -4.70967250147, -4.70966298315, -5.3525,
           -5.1131, -4.9278, -5.1475, -4.4014, -5.4644,
         -5.4, -4.6878, -5.6722, -5.1933, -3.8808,
           -4.9186, -4.7664, -4.3089, -4.315, -4.0897,
         -4.7797, -4.6247, -5.0761, -4.2647, -4.4497,
           -4.5789, -5.3525, -5.1131, -4.9278, -5.1475,
         -4.4014, -5.4644, -5.4, -4.6878, -5.6722,
           -5.1933, -3.8808, -4.9186, -4.7664, -4.3089,
         -4.315, -4.0897, -4.7797, -4.6247, -5.0761,
           -4.2647, -4.4497, -4.5789]
    
    data_dict['GAIA_EPOCH'] = \
        [2015.5, 2015.5, 2015.5, 2015.5, 2015.5,
           2015.5, 2015.5, 2015.5, 2015.5, 2015.5,
         2015.5, 2015.5, 2015.5, 2015.5, 2015.5,
           2015.5, 2015.5, 2015.5, 2015.5, 2015.5,
         2015.5, 2015.5, 2015.5, 2015.5, 2015.5,
           2015.5, 2015.5, 2015.5, 2015.5, 2015.5,
         2015.5, 2015.5, 2015.5, 2015.5, 2015.5,
           2015.5, 2015.5, 2015.5, 2015.5, 2015.5,
         2015.5, 2015.5, 2015.5, 2015.5, 2015.5,
           2015.5, 2015.5, 2015.5]
    
    data_dict['GAIA_PMRA'] = \
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    
    data_dict['GAIA_PMDEC'] = \
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    
    data_dict['GAIA_PARAL'] = \
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    
    data_dict['IFU_PA_REQUEST'] = \
        [0.0, 0.0, 0.0, 0.0, np.nan,
           np.nan, np.nan, np.nan, np.nan, np.nan,
         np.nan, np.nan, np.nan, np.nan, np.nan,
           np.nan, np.nan, np.nan, np.nan, np.nan,
         np.nan, np.nan, np.nan, np.nan, np.nan,
           np.nan, np.nan, np.nan, np.nan, np.nan,
         np.nan, np.nan, np.nan, np.nan, np.nan,
           np.nan, np.nan, np.nan, np.nan, np.nan,
         np.nan, np.nan, np.nan, np.nan, np.nan,
           np.nan, np.nan, np.nan]
    
    data_dict['IFU_DITHER'] = \
        [3, -1, -1, -1, 3, 3, 3, 3, 3, 3,
         3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
         3, 3, 3, 3, 3, 3, 6, 6, 6, 6,
         6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
         6, 6, 6, 6, 6, 6, 6, 6]

    data_dict['HA_RANGE'] = \
        [2.75, 4.0, 3.75, 4.0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
         np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
         np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
         np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
         np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]


    
    return data_dict

    
if __name__ == '__main__':
    
    ############################################################################
    # Set the location of the template and the output file and directory
    
    ifu_driver_template = os.path.join('aux', 'ifu_driver_template.fits')
    
    output_dir = 'output'
    output_filename = os.path.join(output_dir, 'WC_2020A1-ifu_driver_cat.fits')
    
    ############################################################################
    # Create the output directory if it does not exist
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    ############################################################################
    # Get a dictionary with the data which will populate the template
    
    # NOTE: See the above function get_data_dict to understand the structure of
    # the dictionary
    
    data_dict = get_data_dict()
    
    ############################################################################
    # Set the needed information to populate some keywords of the primary header

    trimester, author, report_verbosity, cc_report = set_keywords_info()
    
    ############################################################################
    # Create the IFU driver catalogue
    
    create_ifu_driver_cat(ifu_driver_template, data_dict, output_filename,
                          trimester, author, report_verbosity=report_verbosity,
                          cc_report=cc_report)

