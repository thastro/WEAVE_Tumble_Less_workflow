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


import argparse
import collections
import logging
import os.path

from workflow.utils import create_sub_template
from workflow.utils.get_resources import get_master_cat


def create_ifu_driver_template(catalogue_template, output_filename,
                               update_datetime=False, overwrite=False):
    """
    Create the base template for the IFU driver catalogues.

    Warning: do not overwrite the base template unless you know what
    you are doing! Use this function at your own risk. It will replace
    the aux/ifu_driver_template.fits file that is used to build the 
    survey-specific driver files.

    Parameters
    ----------
    catalogue_template : str
        Any catalogue template containing the SPA columns.
    output_filename : str
        The name of the output file for the new IFU driver template.
    update_datetime : bool, optional
        Update DATETIME keyword in the new IFU driver template.
    overwrite : bool, optional
        Overwrite the output FITS file containing the IFU driver template.
    """

    # Set the list of columns which will be added to the table
    
    col_list = ['TARGSRVY', 'TARGPROG', 'TARGID', 'TARGNAME', 'TARGPRIO',
                'PROGTEMP', 'OBSTEMP', 'GAIA_RA', 'GAIA_DEC', 'GAIA_EPOCH',
                'GAIA_PMRA', 'GAIA_PMDEC', 'GAIA_PARAL', 'IFU_PA', 'IFU_DITHER']

    non_spa_list = {}
    non_spa_list['HA_RANGE'] = {'TFORM':'E',
                                'TUNIT':'h',
                                'TDISP':'F4.1',
                                'TLMIN':0.0,
                                'TLMAX':12.0,
                                'TUCD':'pos.eq.ha',
                                'COMMENT':{
                                    'TTYPE':"Hour angle range",
                                    'TFORM':"data format of field: 4-byte REAL",
                                    'TUNIT':"physical unit of field",
                                    'TDISP':"Display format for column",
                                    'TLMIN':"Minimum value expected for field",
                                    'TLMAX':"Maximum value expected for field",
                                    'TUCD':"UCD for column"}
                                }

    non_spa_list['HA_REQUEST'] = {'TFORM':'E',
                                'TUNIT':'h',
                                'TDISP':'F4.1',
                                'TLMIN':0.0,
                                'TLMAX':12.0,
                                'TUCD':'pos.eq.ha',
                                'COMMENT':{
                                    'TTYPE':"Requested hour angle of observation",
                                    'TFORM':"data format of field: 4-byte REAL",
                                    'TUNIT':"physical unit of field",
                                    'TDISP':"Display format for column",
                                    'TLMIN':"Minimum value expected for field",
                                    'TLMAX':"Maximum value expected for field",
                                    'TUCD':"UCD for column"}
                                }

    
    non_spa_list['TAC_ID'] = {'TFORM':'10A',
                              'TDISP':'A10',
                              'COMMENT':{
                                  'TTYPE':"TAC ID code",
                                  'TFORM':"data format of field: ASCII Character",
                                  'TDISP':"Display format for column"}
                              }

    
    
    # Set the extension name

    extname = 'IFU DRIVER CATALOGUE'

    # Choose the keywords which will be inherited from the catalogue template

    inherited_kwds = ['DATAMVER', 'TRIMESTE', 'DATETIME']

    # Set the information of the new keywords that will be added

    new_primary_kwds = collections.OrderedDict()
    new_primary_kwds['VERBOSE'] = \
        (1, 'Attribute "report_verbosity" of the XML files')
    new_primary_kwds['TACALLOC'] = ('', 'Proposal identifier')
    new_primary_kwds['CAT_NME1'] = ('', 'Submitter forename')
    new_primary_kwds['CAT_NME2'] = ('', 'Submitter surname(s)')
    new_primary_kwds['CAT_MAIL'] = ('', 'Submitter email')
    new_primary_kwds['CAT_CC'] = ('', 'Report cc list')
    
    # Create a dictionary to renaming the column IFU_PA
    
    rename_col_dict = {'IFU_PA': 'IFU_PA_REQUEST'}

    # Create the sub-template with the above paremeters

    create_sub_template(catalogue_template, output_filename, col_list,
                        extname=extname, inherit_primary_kwds=False,
                        inherited_kwds=inherited_kwds,
                        new_primary_kwds=new_primary_kwds,
                        rename_col_dict=rename_col_dict,
                        non_spacols=non_spa_list,
                        update_datetime=update_datetime, overwrite=overwrite)

    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Create a template for the IFU driver catalogues')

    parser.add_argument('--in', dest='catalogue_template',
                        default=os.path.join('aux',
                                             'Master_CatalogueTemplate.fits'),
                        help="""name of catalogue template containing the SPA
                        columns""")

    parser.add_argument('--out', dest='ifu_driver_template',
                        default=os.path.join('aux', 'ifu_driver_template.fits'),
                        help="""name for the output file which will contain the
                        new template for the IFU driver catalogues""")

    parser.add_argument('--update_datetime', action='store_true',
                        help="""update DATETIME keyword from the catalogue
                        template""")

    parser.add_argument('--overwrite', action='store_true',
                        help='overwrite the output file')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()
    
    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    if not os.path.exists(args.catalogue_template):
        logging.info('Downloading the master catalogue template')
        get_master_cat(file_path=args.catalogue_template)

    create_ifu_driver_template(args.catalogue_template,
                               args.ifu_driver_template,
                               update_datetime=args.update_datetime,
                               overwrite=args.overwrite)

