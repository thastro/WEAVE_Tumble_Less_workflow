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

from workflow.utils import populate_fits_table_template
from workflow.utils.get_data_from_xmls import (
    get_datamver_from_xmls, get_halims_from_xmls,
    get_author_from_xmls, get_cc_report_from_xmls, get_trimester_from_xmls,
    get_tacalloc_from_xmls,
    get_spa_data_of_target_random_and_sky_fibres_from_xmls,
    get_spa_data_of_target_calib_random_and_sky_fibres_from_xmls,
    _get_attribute_in_simple_element_of_xml_file)


def _get_col_null_dict_of_template(fits_template):
    
    # Get an empty dictionary for the result
    
    col_null_dict = {}
    
    # Open the template and get list of keywords in extension 1
    
    with fits.open(fits_template) as hdu_list:
    
        hdu = hdu_list[1]
        hdr = hdu.header
        kwd_list = list(hdr.keys())
        
        # For each TTYPE keyword
        
        for kwd in kwd_list:
            
            if re.match('TTYPE([0-9]+)', kwd):
                
                # Get the column number
                
                num_str = kwd[5:]
                
                # Get the values of TTYPE and TFORM
                
                ttype_value = hdr['TTYPE{}'.format(num_str)]
                tform_value = hdr['TFORM{}'.format(num_str)]
                
                # Get the name of the keyword which COULD containe the TNULL
                
                tnull_kwd = 'TNULL{}'.format(num_str)
                
                # Get the TNULL depending on the TFORM value:
                #  - For integers, get it from the TNULL keyword
                #  - For reals, it will be NaN
                #  - For strings, it will be an empty string
                #  - If fails, assign None
                
                try:
                    if (re.match('([1]?)I', tform_value) or
                        re.match('([1]?)J', tform_value) or
                        re.match('([1]?)K', tform_value)):
                        
                        assert tnull_kwd in kwd_list
                        tnull_value = hdr[tnull_kwd]
                        
                    elif (re.match('([1]?)E', tform_value) or
                          re.match('([1]?)D', tform_value)):
                        
                        assert tnull_kwd not in kwd_list
                        tnull_value = np.nan
                        
                    elif re.match('([0-9]*)A', tform_value):
                    
                        assert tnull_kwd not in kwd_list
                        tnull_value = ''
                        
                    else:
                    
                        raise ValueError
                except:
                    tnull_value = None
                
                # Save the column name and its TNULL value
                
                col_null_dict[ttype_value] = tnull_value
    
    return col_null_dict


def _get_col_len(data_dict):

    # Get the lenght from one random column in the dictionary

    random_key = list(data_dict.keys())[0]
    col_len = len(data_dict[random_key])
    
    # Assert that all the columns in the dictionary have the same lenght
    
    for key in data_dict.keys():
        assert len(data_dict[key]) == col_len
    
    return col_len

def _add_missing_gaiadr_from_epoch(data_dict):
    # do this automatically based on the epoch. but *only* those epochs
    # we know about!
    epoch_lookup = {2015.5:'2',2016.0:'3'}
    gaia_dr = []
    dr_missed = 0
    for i in range(len(data_dict['GAIA_DR'])):
        if data_dict['GAIA_DR'][i] == '':
            try:
                gaia_dr.append(epoch_lookup[data_dict['GAIA_EPOCH'][i]])
            except KeyError:
                dr_missed += 1
                gaia_dr.append('')
        else:
            gaia_dr.append(data_dict['GAIA_DR'][i])

    data_dict['GAIA_DR'] = gaia_dr
    if dr_missed:
        logging.warning("Failed to automatically allocate {} GAIA_DR entries based on GAIA_EPOCH".
                        format(dr_missed))
    else:
        logging.info("Added GAIA_DR based on GAIA_EPOCH values")
    
def _add_missing_cols_of_template(data_dict, fits_template):
    
    # Get lenght of each column contained in the data dictionary
    
    col_len = _get_col_len(data_dict)
    
    # Get a dictionary with the names of the columns in the template and their
    # null values
    
    col_null_dict = _get_col_null_dict_of_template(fits_template)
    
    assert None not in col_null_dict.values()
    
    # For each missing column, add a list of nulls to the dictionary    

    out_data_dict = data_dict.copy()

    for col_name in col_null_dict.keys():
        if col_name not in data_dict.keys():
            out_data_dict[col_name] = [col_null_dict[col_name]] * col_len
    
    return out_data_dict


def _add_sim_extension(filename):

    # Open the FITS file in append mode

    with fits.open(filename, mode='append') as hdu_list:

        # Get the number of rows in the first extension

        num_rows = len(hdu_list[1].data)

        # Create the column defitions and its corresponding list of UCDs

        sim_col_list = []
        sim_ucd_list = []

        sim_col_list.append(
            fits.Column(name='SIM_TEMPLATE', format='30A',
                        array=['' for i in range(num_rows)]))
        sim_ucd_list.append('phot.flux.density;meta.file')

        sim_col_list.append(
            fits.Column(name='SIM_MAG', format='D', unit='mag',
                        array=[np.nan for i in range(num_rows)]))
        sim_ucd_list.append('phot.mag')

        sim_col_list.append(
            fits.Column(name='SIM_FILTERID', format='10A',
                        array=['' for i in range(num_rows)]))
        sim_ucd_list.append('instr.filter')

        sim_col_list.append(
            fits.Column(name='SIM_FWHM', format='D', unit='arcsec',
                        array=[np.nan for i in range(num_rows)]))
        sim_ucd_list.append('instr.obsty.seeing')

        sim_col_list.append(
            fits.Column(name='SIM_VELOCITY', format='D', unit='km/s',
                        array=[np.nan for i in range(num_rows)]))
        sim_ucd_list.append('phys.veloc')

        sim_col_list.append(
            fits.Column(name='SIM_REDSHIFT', format='D',
                        array=[np.nan for i in range(num_rows)]))
        sim_ucd_list.append('src.redshift')

        # Create an HDU from the column definitions

        sim_hdu = fits.BinTableHDU.from_columns(sim_col_list)

        # Add TUCD keywords to the HDU in nice positions

        for i, ucd in enumerate(sim_ucd_list):

            ucd_kwd = 'TUCD{}'.format(i + 1)

            unit_kwd = 'TUNIT{}'.format(i + 1)
            form_kwd = 'TFORM{}'.format(i + 1)

            if unit_kwd in sim_hdu.header:
                after_kwd = unit_kwd
            else:
                after_kwd = form_kwd

            sim_hdu.header.set(ucd_kwd, ucd, after=after_kwd)

        # Add an extension name to the HDU

        sim_hdu.name = 'SIM'

        # Append the sim HDU

        hdu_list.append(sim_hdu)


        
    
def create_ifu_fits_cat(fits_template, xml_files, ifu_driver_cat=None,
                        cat_nme1='', cat_nme2='', gaia_dr=None,
                        output_dir='output', output_filename=None,
                        sim_ext=False, add_calib=True, overwrite=False):
    """
    Create an IFU FITS catalogue.
    
    Parameters
    ----------
    fits_template : str
        A FITS template with a primary HDU and a first extension with a table.
    xml_files : str or list of str
        A string with a pattern of the input XML files or list of strings with
        the filenames of the XML files.
    ifu_driver_cat : str, optional
        the filename of the IFU driver for catalogue metadata. Overrides
        CAT_NME1, CAT_NME2 etc
    cat_nme1 : str, optional
        Value for populating CAT_NME1 keyword of the output file.
    cat_nme2 : str, optional
        Value for populating CAT_NME2 keyword of the output file.
    gaia_dr : str, optional
        Value for populating GAIA_DR column of the output file.
    output_dir : str
        Name of the directory which will contain the output IFU FITS catalogue.
    output_filename : str
        The name of the output IFU FITS catalogue (output_dir will be ignored if
        output_filename is not None).
    sim_ext : bool, optional
        Add an extra extension to the output file for the information needed to
        generate the simulations.
    add_calib : bool, optional
        Add fibres with targuse C to the FITS catalogue.
    overwrite : bool, optional
        Overwrite the output FITS file.

    Returns
    -------
    ifu_fits_cat_filename : str
        The name of the output IFU FITS catalogue.
    """
    

    assert os.path.isfile(fits_template), "FITS template {} provided by user does not exist: please download from the WASP at: http://wasp.ast.cam.ac.uk/mysurvey".format(fits_template)

    # Set the name of the output file

    assert '_CatalogueTemplate.fits' in fits_template, "Provided FITS template filename does not contain '_CatalogueTemplate.fits' - please ensure you are using an official WASP template"


    
    if output_filename is None:
    
        prefix = '-'.join(os.path.basename(xml_files[0]).split('-')[:-2])
    
        for xml_file in xml_files:
        
            prefix_of_this_xml_file = (
                '-'.join(os.path.basename(xml_file).split('-')[:-2]))
        
            if prefix_of_this_xml_file != prefix:
                prefix = ''
                break
        
        if prefix != '':
            ifu_fits_cat_filename = os.path.join(
                output_dir, '{}-ifu_from_xmls.fits'.format(prefix))

            # Check the prefix of the fits_template against the above inferred from XMLs
            template_prefix = os.path.basename(fits_template).split('_CatalogueTemplate.fits')[0]
            if not prefix == template_prefix:
                # check for the WEAVE Survey casu - ie prefix='WC_2020A1', template_prefix='WC'
                if prefix.split('_')[0] != template_prefix:
                    assert prefix == template_prefix, "Mismatch between XML prefix ('{}') and FITS template prefix ('{}') - are these for the same survey?".format(prefix,template_prefix)
        else:
            ifu_fits_cat_filename = os.path.join(
                output_dir, 'cat-ifu_from_xmls.fits')
    else:
        
        ifu_fits_cat_filename = output_filename
    
    # Get list of filenames depending on the type of input given for XML files
    
    if type(xml_files) == list:
        xml_filename_list = xml_files
    elif type(xml_files) == str:
        xml_filename_list = glob.glob(xml_files)
        xml_filename_list.sort()
    else:
        raise TypeError



    
    # Check that DATAMVER is consistent bewteen the XML files and the template

    xml_datamver = get_datamver_from_xmls(xml_filename_list)
    template_datamver = fits.getval(fits_template, 'DATAMVER')

    if template_datamver != xml_datamver:
        logging.critical(
            'DATAMVER mismatch ({} != {}) '.format(
                xml_datamver, template_datamver) +
            'between XML files and FITS template: ' +
            'Stop unless you are sure!')
        raise SystemExit(2)
    
    # Defer to the IFU driver cat, if provided (it really should be...) to get
    # the trimester, author and cc_report etc (which should be the
    # same for all them)
    if ifu_driver_cat:
        with fits.open(ifu_driver_cat) as hdu_list:
            logging.info("Extracting metadata from IFU driver {}".format(ifu_driver_cat))
            trimester = hdu_list[0].header['TRIMESTE']
            cat_nme1 = hdu_list[0].header['CAT_NME1']
            cat_nme2 = hdu_list[0].header['CAT_NME2']
            author = hdu_list[0].header['CAT_MAIL']
            cc_report = hdu_list[0].header['CAT_CC']
            tacalloc = hdu_list[0].header['TACALLOC']
            
            # check also the DATAMVER in the driver is consistent..
            driver_datamver = hdu_list[0].header['DATAMVER']     
            if driver_datamver != xml_datamver:
                logging.critical(
                    'DATAMVER mismatch ({} != {}) '.format(
                        xml_datamver, driver_datamver) +
                    'between driver file and FITS template: ' +
                    'Stop unless you are sure!')
                raise SystemExit(2)
            
    else:
        author = get_author_from_xmls(xml_filename_list)
        cc_report = get_cc_report_from_xmls(xml_filename_list)
        trimester = get_trimester_from_xmls(xml_filename_list)
        tacalloc = get_tacalloc_from_xmls(xml_filename_list)
        
        # ...and assume they have entered the cat_nme1, cat_nme2
        assert cat_nme1 != '', "Please provide a forename for the catalogue submitter"
        assert cat_nme2 != '', "Please provide a surname for the catalogue submitter"
    
    primary_kwds = {
        'CAT_NME1': cat_nme1,
        'CAT_NME2': cat_nme2,
        'CAT_MAIL': author,
        'CAT_CC': cc_report,
        'TRIMESTE': trimester,
        'TACALLOC': tacalloc
    }
    
    # Get a dictionary with the SPA data from targets and skies in the XML files
    
    if add_calib is True:
        spa_data_dict = get_spa_data_of_target_calib_random_and_sky_fibres_from_xmls(
                            xml_filename_list, sort_targets=True)
    else:
        logging.warning("Calibration fibres will not be added to FITS catalogue")
        logging.warning("This might cause problems with WASP validation")
        spa_data_dict = get_spa_data_of_target_random_and_sky_fibres_from_xmls(
                            xml_filename_list, sort_targets=True)

        
    # Add missing columns of the template to the dictionary
    
    data_dict = _add_missing_cols_of_template(spa_data_dict, fits_template)
    
    # Assign the requested value to the GAIA_DR column
    
    if gaia_dr is not None:
        data_dict['GAIA_DR'] = [gaia_dr
                                for i in range(len(data_dict['GAIA_DR']))]
    else:
        _add_missing_gaiadr_from_epoch(data_dict)
        
    # Populate the template

    populate_fits_table_template(fits_template, data_dict,
                                 ifu_fits_cat_filename,
                                 primary_kwds=primary_kwds,
                                 update_datetime=True, overwrite=overwrite)
    
    assert os.path.isfile(ifu_fits_cat_filename)
    logging.info("IFU FITS catalogue created: {}".format(ifu_fits_cat_filename))
    # Add an extra extension for the simulations if requested

    if sim_ext == True:
        logging.info("Adding simulation extension for FITS catalogue")
        _add_sim_extension(ifu_fits_cat_filename)
    
    return ifu_fits_cat_filename


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create an IFU FITS catalogue')

    parser.add_argument('template',
                        help="""a FITS template with a primary HDU and a first
                        extension with a table""")

    parser.add_argument('xml_file', nargs='+', help='name of a XML file')

    parser.add_argument('--cat_nme1', default='',
                        help="""value for populating CAT_NME1 keyword of the
                        output file""")

    parser.add_argument('--cat_nme2', default='',
                        help="""value for populating CAT_NME1 keyword of the
                        output file""")

    parser.add_argument('--gaia_dr', default=None,
                        help="""value for populating GAIA_DR column of the
                        output file""")

    parser.add_argument('--driver', dest='ifu_driver_cat', default=None,
                        help="""name of the IFU driver catalogue used in
                        stages1/2 with catalogue metadata (cat_nme1 and
                        cat_nme2 will be ignored if provided)""")

    parser.add_argument('--outdir', dest='output_dir', default='output',
                        help="""name of the directory which will contain the
                        output IFU FITS catalogue""")

    parser.add_argument('--out', dest='output_filename', default=None,
                        help="""name of the file which will contain the output
                        IFU FITS catalogue (outdir value will be ignored if
                        provided)""")

    parser.add_argument('--sim_ext',  action='store_true',
                        help="""add an extra extension to the output file for
                        the information needed to generate the simulations""")

    parser.add_argument('--add_calib', dest='add_calib', action='store_true',
                        help='add fibres with targuse C to the FITS catalogue')

    parser.add_argument('--not_add_calib', dest='add_calib',
                        action='store_false',
                        help="""not add fibres with targuse C to the FITS
                        catalogue""")
    parser.set_defaults(add_calib=True)

    parser.add_argument('--overwrite', action='store_true',
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
        logging.info('Creating the output directory')
        os.mkdir(dirname)

    create_ifu_fits_cat(args.template, args.xml_file,
                        ifu_driver_cat=args.ifu_driver_cat,
                        cat_nme1=args.cat_nme1, cat_nme2=args.cat_nme2,
                        gaia_dr=args.gaia_dr, output_dir=args.output_dir,
                        output_filename=args.output_filename,
                        sim_ext=args.sim_ext, add_calib=args.add_calib,
                        overwrite=args.overwrite)

