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
import os
import re
import numpy as np

from astropy.io import fits

from workflow.utils import check_equal_headers
from workflow.utils import populate_fits_table_template

def _filter_bundles(data, nfibres=37, central_spaxel_id='C04'):
    """
    Filter non-science mIFU bundles from fibre data.
    
    Output FITS catalogues should not contain mIFU bundles
    used for calibration or sky subtraction. This function
    first groups the fibres into bundles, and then determines
    if they are calib/sky bundles.

    The only instance where such bundles are retained within
    the target catalogue are when a non-standard processing is
    requested - for example, sky fibres with APS_FLAG values
    indicating additional APS analyis the the fibres. These 
    directives are communicated to the APS via the input FITS
    catalogues.

    Whilst users can in principle modify the calibration bundle,
    WASP may well reject the entry, particularly if the WD-sourced
    target data are changed. Modification of the sky targets 
    surrounding the WD might be permitted, subject to review.

    Parameters
    ----------
    data : astropy.io.fits.fitsrec.FITS_rec 
        Name of a FITS file containing a MOS catalogue.

    Returns
    -------
    _data : astropy.io.fits.fitsrec.FITS_rec
        Input data, with all bundles described above filtered out
    """

    selection = np.array([True for row in data])

    mifu_indices = np.where(np.array([int(row['PROGTEMP'][0]) > 6 for row in data]))[0]
    if not len(mifu_indices):
        return data
    bundles = []

    bundle_number = data['IFU_SPAXEL'][mifu_indices[0]][1:3]
    bundle = {'bundle':bundle_number,'indices':[],'central':None,'include':True,'spaxel_ids':[]}

    for index in mifu_indices:
        if (data[index]['IFU_SPAXEL'][1:3] == bundle_number) and (index != mifu_indices[-1]):
            bundle['indices'].append(index)
            bundle['spaxel_ids'].append(data[index]['IFU_SPAXEL'])
        else:
            # moved onto the next bundle - store the existing one and move on
            if index == mifu_indices[-1]:
                # this is the end of the list, add the final index
                bundle['indices'].append(index)
                bundle['spaxel_ids'].append(data[index]['IFU_SPAXEL'])
            bundle['indices'] = np.array(bundle['indices'])
            bundle['spaxel_ids'] = np.array(bundle['spaxel_ids'])
            assert len(bundle['indices']) == nfibres, "The number of fibres in bundle {} ({}) is not {}".format(bundle['bundle'],len(bundle['indices']),nfibres)
            assert len(np.unique(bundle['spaxel_ids'])) == nfibres, "There are duplicate spaxelIDs in bundle {}:\n{}".format(bundle['bundle'],str(bundle['spaxel_ids']))
            bundles.append(bundle)
            
            bundle_number = data['IFU_SPAXEL'][index][1:3]
            bundle = {'bundle':bundle_number,'indices':[index],'central':None,'include':True,'spaxel_ids':[data[index]['IFU_SPAXEL']]}

        if data[index]['IFU_SPAXEL'][3:6] == central_spaxel_id:
            bundle['central'] = index

    nrejected = 0
    for bundle in bundles:
        central_spaxel = data[bundle['central']]
        
        if central_spaxel['TARGUSE'] in ['C','R','S']:
            retain = False
            bundle_data = data[bundle['indices']]
            # This is a calib or sky bundle. Should only be included in FITS
            # catalogue if non-std analysis is specified for sky spaxels - this
            # is via:
            # - TARGCLASS != '' for TARGUSE=S
            # - APS_* inputs are non-null '' for TARGUSE=S
            # The TARGUSE='C','R' spaxels should inherit from the
            # WD catalogue and not be tampered with

            retention_reasons = []            
            if 'T' in bundle_data['TARGUSE']:
                # maybe this is an 'offset' target bundle with non-target in the centre
                # spaxel - do not filter this out!
                retention_reasons.append('(T)arget spaxel in bundle')
                retain = True
                
            sky_spaxels = bundle['indices'][np.where(bundle_data['TARGUSE'] == 'S')]
            sky_targclass = data[sky_spaxels]['TARGCLASS']


            

            if True in [targclass != 'SKY' for targclass in sky_targclass]:
                retain = True
                retention_reasons.append('TARGCLASS')
            
            for column in data.columns:
                if 'APS_' in column.name:
                    if column.null:
                        nullval = column.null
                    elif 'A' in column.format:
                        nullval = ''
                    elif column.format in ['I','J','K']:
                        nullval = -1
                    else:
                        nullval = np.nan
                        
                    if (nullval in ['',-1]) or not np.isnan(nullval):
                        test = np.array([spaxel_val != nullval for spaxel_val in data[sky_spaxels][column.name]])
                    else:
                        test = np.array([not np.isnan(spaxel_val) for spaxel_val in data[sky_spaxels][column.name]])
                    if True in test:
                        retention_reasons.append(column.name)
                        retain = True

            if len(retention_reasons):
                logging.debug("")
                logging.debug("Retaining bundle #{} [central spaxel @i={}]: custom entries: {} ".format(bundle['bundle'],bundle['central'],','.join(retention_reasons)))

            selection[bundle['indices']] = retain
            if not retain:
                nrejected += 1
                bundle['include'] = False
                bundle_type = 'sky bundle'
                if ('R' in  bundle_data['TARGUSE']) or ('C' in bundle_data['TARGUSE']):
                    bundle_type = 'calib bundle'

                logging.debug("\n")
                logging.debug("Removing: [{}]".format(bundle_type))
                for i in bundle['indices']:
                    central_str = ''
                    if i == bundle['central']:
                        central_str = ' [central spaxel]'
                    logging.debug("\t @i={}: {} {} {} {}{}".format(i,data['IFU_SPAXEL'][i],data['TARGUSE'][i],data['TARGCLASS'][i],data['APS_FLAG'][i],central_str))
            
    _data = data[selection]
    logging.debug("Filtered {}/{} mIFU bundles that are either calib or sky bundles".format(nrejected,len(bundles)))
    return _data
    
def _get_data_dict(fits_filename,filter=True):
    
    data_dict = {}
    
    with fits.open(fits_filename) as hdu_list:
        data = hdu_list[1].data
        if (filter) and len(data['CNAME']):
            data = _filter_bundles(data)

            
        for col in data.names:
            data_dict[col] = list(data[col])
    return data_dict


def _get_combo_data_dict(fits_filename1, fits_filename2, filter=True):
    assert (fits_filename1) or (fits_filename1)

    if fits_filename1:
        data_dict1 = _get_data_dict(fits_filename1,filter=filter)
        data_keys = data_dict1.keys()
    if fits_filename2:
        data_dict2 = _get_data_dict(fits_filename2,filter=filter)
        data_keys = data_dict2.keys()
        
    combo_data_dict = {}
    
    for key in data_keys:
        _data = []
        if fits_filename1:
            _data += data_dict1[key]
        if fits_filename2:
            _data += data_dict2[key]
        combo_data_dict[key] = _data
    
    return combo_data_dict


def create_filtered_ifu_fits_cat(ifu_cat, output_dir='output',
                                 output_filename=None, overwrite=False):
    """
    Create a filtered IFU catalogue (where no MOS counterpart exists)
    
    Parameters
    ----------
    ifu_cat : str
        Name of a FITS file containing a IFU catalogue.
    output_dir : str, optional
        The name of the directory for the output combo catalogue.
    output_filename : str, optional
        The name of the output combo catalogue (outdir value will be ignored if
        provided).
    overwrite : bool, optional
        Overwrite the output FITS file.

    Returns
    -------
    filtered_fits_cat_filename : str
        The name of the output combo FITS catalogue.
    """
    
    # Set the name of the output file

    if output_filename is None:
        
        ifu_cat_basename = os.path.basename(ifu_cat)
    
        if ifu_cat_basename.endswith('-ifu.fits'):
            
            prefix = ifu_cat_basename[:-9]
            
            filtered_fits_cat_filename = os.path.join(
                output_dir, '{}.fits'.format(prefix))
        
        else:
        
            filtered_fits_cat_filename = os.path.join(
                output_dir, 'filtered_cat.fits')
    
    else:
        
        filtered_fits_cat_filename = output_filename

    logging.debug("Output filename is: {}".format(filtered_fits_cat_filename))
        
    # Create the filtered FITS catalogue
    
    filtered_data_dict = _get_combo_data_dict(ifu_cat,None,filter=True)

    populate_fits_table_template(ifu_cat, filtered_data_dict,
                                 filtered_fits_cat_filename,
                                 update_datetime=True, overwrite=overwrite)
    
    return filtered_fits_cat_filename



def create_combo_fits_cat(mos_cat, ifu_cat, output_dir='output',
                          output_filename=None, overwrite=False):
    """
    Combine MOS and IFU catalogues to create a combo catalogue.
    
    Parameters
    ----------
    mos_cat : str
        Name of a FITS file containing a MOS catalogue.
    ifu_cat : str
        Name of a FITS file containing a IFU catalogue.
    output_dir : str, optional
        The name of the directory for the output combo catalogue.
    output_filename : str, optional
        The name of the output combo catalogue (outdir value will be ignored if
        provided).
    overwrite : bool, optional
        Overwrite the output FITS file.

    Returns
    -------
    combo_fits_cat_filename : str
        The name of the output combo FITS catalogue.
    """
    
    # Set the name of the output file
    
    if output_filename is None:
        
        mos_cat_basename = os.path.basename(mos_cat)
        ifu_cat_basename = os.path.basename(ifu_cat)
    
        if ((mos_cat_basename.endswith('-mos.fits') or
             mos_cat_basename.endswith('-ifu.fits')) and
            (ifu_cat_basename.endswith('-mos.fits') or
             ifu_cat_basename.endswith('-ifu.fits'))):
            
            prefix = mos_cat_basename[:-9]
            
            combo_fits_cat_filename = os.path.join(
                output_dir, '{}.fits'.format(prefix))
        
        else:
        
            combo_fits_cat_filename = os.path.join(
                output_dir, 'combo_cat.fits')
    
    else:
        
        combo_fits_cat_filename = output_filename
    
    # Check that the headers of both files are equal
    
    equal_headers = check_equal_headers(mos_cat, ifu_cat,
                                        ignore_values=['DATETIME'])
    
    if not equal_headers:
        logging.error('equal headers are expected to create a combo FITS file')
        
        assert equal_headers

        
    # Create the combo FITS catalogue
    
    combo_data_dict = _get_combo_data_dict(mos_cat, ifu_cat)

    populate_fits_table_template(mos_cat, combo_data_dict,
                                 combo_fits_cat_filename,
                                 update_datetime=True, overwrite=overwrite)
    
    return combo_fits_cat_filename


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create a combo catalogue')

    parser.add_argument('ifu_cat',
                        help='a FITS file containing an IFU catalogue')

    parser.add_argument('--mos_cat', default=None,
                        help='a FITS file containing a MOS catalogue')
    
    parser.add_argument('--outdir', dest='output_dir', default='output',
                        help="""name for the directory which will contain the
                        combo catalogue""")

    parser.add_argument('--out', dest='output_filename', default=None,
                        help="""name for the output file which will contain the
                        combo catalogue (outdir value will be ignored if
                        provided)""")

    parser.add_argument('--overwrite', dest='overwrite', action='store_true',
                        help='overwrite the output file')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    if args.output_filename is None:
        dirname = args.output_dir
    else:
        dirname = os.path.dirname(args.output_filename)
    
    if not os.path.exists(dirname):
        logging.info('Creating the output directory {}'.format(dirname))
        os.mkdir(dirname)

    if not args.mos_cat:
        out = create_filtered_ifu_fits_cat(args.ifu_cat,
                                           output_dir=args.output_dir,
                                           output_filename=args.output_filename,
                                           overwrite=args.overwrite)
    else:
        out = create_combo_fits_cat(args.mos_cat, args.ifu_cat,
                                    output_dir=args.output_dir,
                                    output_filename=args.output_filename,
                                    overwrite=args.overwrite)
        

    logging.info("Combined catalogue written to {}".format(out))
    logging.info("")
    logging.info("----- Please submit this ^^ to the WASP (with the correct filename) -----")
    
