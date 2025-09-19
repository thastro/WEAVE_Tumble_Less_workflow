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


import datetime as _datetime

from astropy.io import fits as _fits
import astropy.time
import numpy as _np


def _try_to_copy_nonspa_comment(kwd_type, ttype, hdu, i, non_spacols):

    kwd_i = '{}{}'.format(kwd_type, i + 1)

    if kwd_type in non_spacols[ttype]['COMMENT']:

        hdu.header.comments[kwd_i] = non_spacols[ttype]['COMMENT'][kwd_type]


def _try_to_copy_comment(kwd_type, hdu, i, template_hdu, j):

    kwd_i = '{}{}'.format(kwd_type, i + 1)
    kwd_j = '{}{}'.format(kwd_type, j + 1)

    if kwd_i in hdu.header.keys():

        hdu.header.comments[kwd_i] = template_hdu.header.comments[kwd_j]


def _try_to_copy_nonspa_keyword(kwd_type, ttype, hdu, i, non_spacols):

    kwd_i = '{}{}'.format(kwd_type, i + 1) 

    if kwd_type in non_spacols:

        # The requested keyword is in the non_spacols dictionary - OK to add the kwd_type

        next_ttype_kwd = 'TTYPE{}'.format(i + 2)

        # Copy the keyword before the next TTYPE kwd if possible, otherwise at
        # the end of the header
        
        if next_ttype_kwd in hdu.header.keys():
            hdu.header.insert(next_ttype_kwd,
                              (kwd_i,non_spacols[ttype][kwd_type]))
            
def _try_to_copy_keyword(kwd_type, hdu, i, template_hdu, j):

    kwd_i = '{}{}'.format(kwd_type, i + 1)
    kwd_j = '{}{}'.format(kwd_type, j + 1)

    if kwd_j in template_hdu.header.keys():

        # Copy the keyword before the next TTYPE kwd if possible, otherwise at
        # the end of the header

        next_ttype_kwd = 'TTYPE{}'.format(i + 2)

        if next_ttype_kwd in hdu.header.keys():
            hdu.header.insert(next_ttype_kwd,
                              (kwd_i, template_hdu.header[kwd_j]))
        else:
            hdu.header[kwd_i] = template_hdu.header[kwd_j]

        # Copy the comment of the keyword

        hdu.header.comments[kwd_i] = template_hdu.header.comments[kwd_j]


def create_sub_template(catalogue_template, output_filename, col_list,
                        extname=None, inherit_primary_kwds=True,
                        inherited_kwds=[], new_primary_kwds={},
                        rename_col_dict={}, non_spacols={}, update_datetime=True,
                        checksum=True, overwrite=False):
    """
    Create a template from a subset of columns from other template.

    Parameters
    ----------
    catalogue_template : str
        The name of the catalogue template.
    output_filename : str
        The name of the output file for the new sub-template.
    col_list : list of str
        A list containing all the columns which will be included in the new
        sub-template.
    extname : str, optional
        Name to be used in the table extension.
    inherit_primary_kwds : bool, optional
        A boolean value to indicate whether all keywords in the primary header
        will be inherited.
    inherited_kwds : list of str, optional
        A list of keywords which will be inherited when not all the keywords are
        inherited.
    new_primary_kwds: dict, optional
        A dictionary with the keywords, their values and their comments which
        should be added to the primary header. The keywords are specified via
        the dictionary keys, while the values of the dictionary should contain
        a tuple with the value and comment of each keyword (None is used to
        indicate no comment). It can be used to overwrite inherited keywords
        (in that case, if the comment is None, the original comment will be
        preserved). An ordered dictionary could be used to indicate the order of
        these keywords.
    rename_col_dict : dict, optional
        A dictionary used to rename the columns. Its keys should be the name of
        the columns in the input catalogue template, while its values should be
        the new name which will be used in the new sub-template.
    update_datetime : bool, optional
        Update DATETIME keyword in the new sub-template.
    checksum : bool, optional
        Add CHECKSUM and DATASUM keywords in the new sub-template.
    overwrite : bool, optional
        Overwrite the output FITS file containing the sub-template.
    """

    # Lists with the type of keywords that will be copied from the catalogue
    # template

    basic_kwd_type_list = ['TTYPE', 'TFORM', 'TDISP', 'TUNIT', 'TNULL']
    extra_kwd_type_list = ['TLMIN', 'TLMAX', 'TUCD', 'TPROP']

    # Read the catalogue template

    template_hdulist = _fits.open(catalogue_template)
    template_primary_hdr = template_hdulist[0].header
    template_hdu = template_hdulist[1]

    # Check that all the requested columns exist in the catalogue template and
    # that there are not repeated columns

    template_column_names = [col.name for col in template_hdu.columns]

    for col_name in col_list:
        assert col_name in template_column_names

    assert len(col_list) == len(set(col_list))

    # Create the column list and save the mapping between the columns of the
    # new sub-template and the input catalogue template

    column_list = []
    column_mapping = {}

    col_counter = 0

    for j, col in enumerate(template_hdu.columns):

        if col.name in col_list:

            # Save the mapping and increase the counter of created columns

            column_mapping[col_counter] = j
            col_counter += 1

            # Get the properties of the column

            if col.name in rename_col_dict.keys():
                col_name = rename_col_dict[col.name]
            else:
                col_name = col.name
            
            col_format = col.format
            
            col_disp = col.disp
            col_unit = col.unit
            col_null = col.null

            # Create the column and add it to the column list

            column = _fits.Column(name=col_name, format=col_format,
                                  disp=col_disp, unit=col_unit, null=col_null)

            column_list.append(column)

    # Now add any non-SPA columns in here:

    for _col in non_spacols:
        j += 1
        column_mapping[col_counter] = j
        col_counter += 1
        
        col_name = _col

        col_unit = None
        col_null = None

        col_format = non_spacols[_col]['TFORM']
        col_disp = non_spacols[_col]['TDISP']
        if 'TUNIT' in non_spacols[_col]:
            col_unit = non_spacols[_col]['TUNIT']
        
        if 'TNULL' in non_spacols[_col]:
            col_null = non_spacols[_col]['TNULL']

        # Create the column and add it to the column list
        
        column = _fits.Column(name=col_name, format=col_format,
                              disp=col_disp, unit=col_unit, null=col_null)
        
        column_list.append(column)
        col_list.append(col_name)

        
    # Create the HDU from the column list
    
    coldefs = _fits.ColDefs(column_list)
    hdu = _fits.BinTableHDU.from_columns(coldefs)

    # Copy the comments of the original catalogue for the created keywords

    for i in range(len(col_list)):

        for kwd_type in basic_kwd_type_list:


            if col_list[i] in non_spacols: 
                _try_to_copy_nonspa_comment(kwd_type, col_list[i], hdu, i, non_spacols)                
            else:
                j = column_mapping[i]
                _try_to_copy_comment(kwd_type, hdu, i, template_hdu, j)


    
    # Add other keywords related with each column which could not be added
    # in the defitions of the columns

    for i in range(len(col_list)):

        for kwd_type in extra_kwd_type_list:

            if col_list[i] in non_spacols:
                _try_to_copy_nonspa_keyword(kwd_type, col_list[i], hdu, i, non_spacols)
            else:
                j = column_mapping[i]
                _try_to_copy_keyword(kwd_type, hdu, i, template_hdu, j)

    # Give a name to the HDU

    if extname is not None:
        hdu.name = extname
    else:
        hdu.name = 'CATALOGUE'
    
    # Create the primary header copying its keywords from the template
    
    primary_hdr = _fits.Header()
    
    basic_kwd_list =['SIMPLE', 'BITPIX', 'NAXIS', 'EXTEND', 'COMMENT']
    
    for kwd in template_primary_hdr.keys():
        if kwd not in basic_kwd_list:
            if (inherit_primary_kwds is True) or (kwd in inherited_kwds):
                primary_hdr[kwd] = template_primary_hdr[kwd]
                primary_hdr.comments[kwd] = template_primary_hdr.comments[kwd]
    
    # Add/Overwrite the requested keywords (before DATETIME please)
    
    for kwd in new_primary_kwds.keys():
        
        value, comment = new_primary_kwds[kwd]
        
        
        if kwd in primary_hdr:
            del primary_hdr[kwd]

        primary_hdr.insert('DATETIME',(kwd,value,comment),after=False) 
            
        # primary_hdr[kwd] = value

        # if comment is not None:
        #     primary_hdr.comments[kwd] = comment

    # Update the keyword DATETIME if requested (and it exists)
    now = astropy.time.Time(astropy.time.Time.now(), format='isot', scale='utc').fits.replace('(UTC)','')
    
    if (update_datetime is True) and ('DATETIME' in primary_hdr.keys()):
        primary_hdr['DATETIME'] = now

    # Create the primary HDU
    
    primary_hdu = _fits.PrimaryHDU(header=primary_hdr)
    
    # Create a HDU list and save it to a file
    
    hdulist = _fits.HDUList([primary_hdu, hdu])

    if checksum:
        # base it on the "now" time so they are connected (so checksum
        # in the writeto function call is always false
        datasums = [ext.add_datasum(when=now) for ext in hdulist]
        checksums = [ext.add_checksum(when=now,override_datasum=False) for ext in hdulist]
        
        
    hdulist.writeto(output_filename, checksum=False, overwrite=overwrite)

