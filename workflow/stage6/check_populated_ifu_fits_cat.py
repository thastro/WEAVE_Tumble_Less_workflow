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
import logging

from astropy.io import fits


def check_populated_ifu_fits_cat(ifu_cat_from_xmls, ifu_cat):
    """
    Check an IFU FITS catalogue to detect undesired changes.

    Parameters
    ----------
    ifu_cat_from_xmls : str
        A non-populated IFU FITS catalogue.
    ifu_cat : str
        A populated IFU FITS catalogue.

    Returns
    -------
    result : bool
        True if the file passes the check, otherwise False.
    """
    
    result = True

    # Set the parameters to be ignored in the comparison
    
    ignore_hdus = ['SIM']
    
    ignore_keywords = [
        'CHECKSUM', 'DATASUM', 'DATETIME', 'MAG_G_CM', 'MAG_R_CM', 'MAG_I_CM']
    
    non_ignore_field = [
        'CNAME', 'TARGSRVY', 'TARGPROG', 'TARGCAT', 'TARGID', 'TARGNAME',
        'TARGPRIO', 'PROGTEMP', 'OBSTEMP', 'GAIA_ID', 'GAIA_DR', 'GAIA_RA',
        'GAIA_DEC', 'GAIA_EPOCH', 'GAIA_PMRA', 'GAIA_PMRA_ERR', 'GAIA_PMDEC',
        'GAIA_PMDEC_ERR', 'GAIA_PARAL', 'GAIA_PARAL_ERR', 'HEALPIX',
        'IFU_SPAXEL', 'IFU_PA', 'IFU_DITHER']

    with fits.open(ifu_cat_from_xmls) as hdu_list:
        ignore_fields = [colname
                          for colname in hdu_list[1].data.names
                              if colname not in non_ignore_field]
    
    # Compare the FITS files and set the result
    
    diff = fits.diff.FITSDiff(ifu_cat_from_xmls, ifu_cat,
                              ignore_hdus=ignore_hdus,
                              ignore_keywords=ignore_keywords,
                              ignore_fields=ignore_fields)

    result = diff.identical
    
    # Print the report
    
    report = diff.report()
    
    for report_line in report.split('\n'):
        logging.info(report_line)

    # Write a message

    if result == True:
        logging.info(
            '{} has passed the check against {}'.format(
                ifu_cat, ifu_cat_from_xmls))
    else:
        logging.error(
            '{} has NOT passed the check against {}'.format(
                ifu_cat, ifu_cat_from_xmls))

    return result


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        description='Check an IFU FITS catalogue to detect undesired changes')
    
    parser.add_argument('ifu_cat_from_xmls',
                        help='a non-populated IFU FITS catalogue')
    
    parser.add_argument('ifu_cat',
                        help='a populated IFU FITS catalogue')
    
    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')
    
    args = parser.parse_args()
    
    logging.basicConfig(level=getattr(logging, args.log_level.upper()))
    
    check_populated_ifu_fits_cat(args.ifu_cat_from_xmls, args.ifu_cat)

