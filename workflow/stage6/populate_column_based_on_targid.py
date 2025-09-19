#!/usr/bin/env python3

#
# Copyright (C) 2021 Scott C. Trager
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

"""
Populate a column depending on the corresponding TARGID values.

Populate a column depending on the corresponding TARGID values in an IFU FITS
catalogue. Assumes each target has a single set of parameters. Ignores all
fibres marked as TARGUSE="S".

Author: sctrager@astro.rug.nl
"""

import argparse
import logging

import numpy as np
from astropy.io import fits


def populate_column_based_on_targid(ifu_cat, nested_dict):
    """
    Populate a column depending on the corresponding TARGID values.
    
    Populate a column depending on the corresponding TARGID values in an IFU
    FITS catalogue. Assumes each target has a single set of parameters. Ignores
    all fibres marked as TARGUSE="S".

    Parameters
    ----------
    ifu_cat : str
        Name of a FITS file containing a IFU catalogue.
    nested_dict : dict
        A dictionary containing {col : dict} flags, where dict is a dictionary
        of {targid : value} pairs
    """
    logging.info('Adding APS_ parameters to {}'.format(ifu_cat))

    with fits.open(ifu_cat, mode='update') as hdu_list:

        # get the data from the catalogue

        for galaxy in np.unique(hdu_list[1].data['TARGID']):
            rows = np.where(hdu_list[1].data['TARGID'] == galaxy)
            logging.info('Modifying APS_ values for {}'.format(galaxy))
            try:
                galaxy=int(galaxy)
            except:
                pass

            targclass=hdu_list[1].data['TARGCLASS'][rows]
            for aps_key in nested_dict.keys():
                try:
                    logging.debug('nested_dict[{}][{}] = {}'.format(aps_key,galaxy,nested_dict[aps_key][galaxy]))
                    if ('_DR' in aps_key) or ('_ID' in aps_key):
                        aps_default=np.where(np.char.equal(targclass,'GALAXY'),
                                             nested_dict[aps_key][galaxy],'')
                    else:
                        aps_default=np.where(np.char.equal(targclass,'GALAXY'),
                                             nested_dict[aps_key][galaxy],None)
                    hdu_list[1].data[aps_key][rows]=aps_default
                except:
                    pass

        logging.info('Saving {}'.format(ifu_cat))
        hdu_list.flush()

