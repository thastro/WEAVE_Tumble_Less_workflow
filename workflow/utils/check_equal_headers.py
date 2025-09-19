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


import logging
import re as _re

from astropy.io import fits as _fits


def _get_num_ext_in_fits_file(filename):
    
    with _fits.open(filename) as hdu_list:
        num_ext = len(hdu_list)
    
    return num_ext


def _match_in_regex_list(string, regex_list):
    
    result = False
    
    for regex in regex_list:
        if _re.match(regex, string):
            result = True
            break
    
    return result


def _check_equal_one_header(fits_filename1, fits_filename2, ext_i,
                            ignore_values=[]):
    
    # Get the list of extra keywords which will be ignored in the comparison
    
    extra_kwd_list = ['COMMENT', 'HISTORY']
    
    # Get the list of the integrity keywords (which will be ignored in the
    # comparison of the values)
    
    integrity_kwd_list = ['CHECKSUM', 'DATASUM']
    
    # Get the list of regular expresions for the basic keywords which values
    # will not be compared
    
    basic_kwd_regex_list =['SIMPLE', 'BITPIX', 'NAXIS', 'EXTEND', 
                           'XTENSION', 'PCOUNT', 'GCOUNT', 'NAXIS([0-9]+)']

    # Get the headers of each file in the requested extension

    hdr1 = _fits.getheader(fits_filename1, ext_i)
    hdr2 = _fits.getheader(fits_filename2, ext_i)
    
    # Get the list of keywords in the headers of the requested extension
    
    kwd_list1 = list(hdr1.keys())
    kwd_list2 = list(hdr2.keys())
    
    # Guess if both files contain the same keywords in the requested extensions
    
    result = True
    
    for kwd in kwd_list1:
        if (kwd not in kwd_list2) and (kwd not in extra_kwd_list):
            logging.error('keyword {} of {} is missing in {}'.format(
                          kwd, fits_filename1, fits_filename2))
            result = False
    
    for kwd in kwd_list2:
        if (kwd not in kwd_list2) and (kwd not in extra_kwd_list):
            logging.error('keyword {} of {} is missing in {}'.format(
                          kwd, fits_filename2, fits_filename1))
            result = False
    
    # If both files have the same keywords in the requested extension, check
    # whether their values are equal
    
    if result == True:
    
        for kwd in kwd_list1:
        
            if ((not _match_in_regex_list(kwd, basic_kwd_regex_list)) and
                (kwd not in extra_kwd_list)):
            
                if ((kwd not in ignore_values) and
                    (kwd not in integrity_kwd_list)):
                    
                    if hdr1[kwd] != hdr2[kwd]:
                        logging.error(
                          'keyword {} has different values in {} and {}'.format(
                             kwd, fits_filename1, fits_filename2))
                        result = False
    
    return result


def check_equal_headers(fits_filename1, fits_filename2, ignore_values=[]):
    """
    Check that the headers of two FITS files are equal.

    Parameters
    ----------
    fits_filename1 : str
        The name of the first FITS file.
    fits_filename2 : str
        The name of the second FITS file.
    ignore_values : list of str
        A list of the keywords in the header which values could be different
        (in addition to the basic ones, like NAXISi).
    
    Returns
    -------
    result : bool
        A boolean value indicating whether the headers are equal (True) or
        different (False).
    """
    
    # Check whether both files have the same number of extensions
    
    num_ext1 = _get_num_ext_in_fits_file(fits_filename1)
    num_ext2 = _get_num_ext_in_fits_file(fits_filename2)
    
    if num_ext1 == num_ext2:
        result = True
    else:
        logging.error('{} and {} have different number of extensions'.format(
                      fits_filename1, fits_filename2))
        result = False
    
    # If both files have the same number of extensions, check whether the header
    # of each extension are equal
    
    if result == True:
    
        for ext_i in range(num_ext1):
        
            result_i = _check_equal_one_header(fits_filename1, fits_filename2,
                                               ext_i,
                                               ignore_values=ignore_values)
            
            if result_i == False:
                logging.error(
                    'headers of {} and {} are different in extension {}'.format(
                        fits_filename1, fits_filename2, ext_i))
                result = False
    
    return result

