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
import numpy as _np
from astropy.io import fits as _fits


def populate_fits_table_template(fits_template, data_dict, output_filename,
                                 primary_kwds={}, update_datetime=True,
                                 checksum=True, overwrite=False):
    """
    Populate a FITS table template with the provided data.

    Parameters
    ----------
    fits_template : list of str
        A FITS template with a primary HDU and a first extension with a table.
    data_dict : dict
        A dictionary with the data. Its keys should contain the name of the
        columns of the table in the first extension of the FITS template. Its
        values should be array-like with the data to populate the table.
    output_filename : str
        The name of the output file which will be created.
    primary_kwds : dict, optional
        A dictionary with a list of keywords and their corresponding values
        which will be written a in the primary header (updated or added).
    update_datetime : bool, optional
        Update DATETIME keyword in the output file.
    checksum : bool, optional
        Add CHECKSUM and DATASUM keywords in the output file.
    overwrite : bool, optional
        Overwrite the output FITS file.
    """
    
    # Read the FITS template
    
    template_hdulist = _fits.open(fits_template)
    template_primary_hdu = template_hdulist[0]
    template_hdu = template_hdulist[1]
    
    # Check that all the columns are available in the dictionary with the data
    
    template_column_names = [col.name for col in template_hdu.columns]

    for col_name in template_column_names:
        assert col_name in data_dict.keys(), "Template column {} is not present in the data".format(col_name)
        
    # Create list of columns as described in the template populated with the
    # provided data
    
    column_list = []
    
    for col in template_hdu.columns:
        _data = data_dict[col.name]
        column = _fits.Column(name=col.name, format=col.format,
                             disp=col.disp, unit=col.unit, null=col.null,
                              array=_data)
        
        column_list.append(column)
    
    # Create a HDU from the column list
    
    coldefs = _fits.ColDefs(column_list)
    hdu = _fits.BinTableHDU.from_columns(coldefs)
    
    # Copy the header from the template
    
    hdu.header = template_hdu.header
    
    # Create the primary extension and populate it with the provided information
    
    primary_hdu = template_primary_hdu
    
    for kwd in primary_kwds.keys():
        primary_hdu.header[kwd] = primary_kwds[kwd]
            
    # Update the keyword DATETIME if requested (and it exists)
    
    if (update_datetime is True) and ('DATETIME' in primary_hdu.header.keys()):
        datetime_str = _datetime.datetime.utcnow().strftime(
                           '%Y-%m-%dT%H:%M:%S.%f')
        primary_hdu.header['DATETIME'] = datetime_str
    
    # Create a HDU list and save it to a file
    
    hdulist = _fits.HDUList([primary_hdu, hdu])
    
    hdulist.writeto(output_filename, checksum=checksum, overwrite=overwrite)

