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
import sys
from collections import OrderedDict, Counter

import numpy as np
from astropy.io import fits

from workflow.utils import check_equal_headers
from workflow.utils.get_obstemp_info import (get_obstemp_info,
                                             get_obstemp_dict)
from workflow.utils.get_progtemp_info import (get_progtemp_dict,
                                              get_progtemp_info,
                                              get_obsmode_from_progtemp)
from workflow.utils.get_resources import (get_progtemp_file,
                                          get_obstemp_file)
from workflow.utils.get_obslimit_assets import (get_obslimits_tool,
                                               get_configure_cfg)


def _check_versus_template(cat_filename, template):
    
    """
    A function to check the header values between the catalogue file and the template file are consistent.

    :param cat_filename: (Str) the name of the catalogue file.
    :param template: (Str) the name of the template file.
    :return: (bool) result.
    """

    ignore_values = ['TRIMESTE', 'DATETIME', 'VERBOSE', 'TACALLOC',
                     'CAT_NME1', 'CAT_NME2', 'CAT_MAIL', 'CAT_CC']

    result = check_equal_headers(cat_filename, template,
                                 ignore_values=ignore_values)
    
    return result


def _check_trimester(trimester):
    
    """
    A function to check the trimester header value formatting.

    :param trimester: (Str) the trimester value to be validated.
    :return: (bool) result.
    """

    if re.match('^[0-9]{4}[AB][12]$', trimester):
        result = True
    else:
        result = False

    return result


def _check_verbose(verbose):
    
    """
    A function to check the verbose mode value.

    :param verbose: (int) the verbose value to be validated.
    :return: (bool) result.
    """

    result = (verbose in [0, 1])

    return result


def _check_email(string):
    
    """
    A function to check the formatting of a supplied email address.

    :param string: (Str) the email address to be validated.
    :return: (bool) result.
    """

    arroba_index = string.find('@')

    if arroba_index != -1:

        substring = string[arroba_index + 1:]
        dot_index = substring.find('.')

        if dot_index != -1:
            result = True
        else:
            result = False

    else:
        result = False

    return result


def _check_author(author):
    
    """
    A function to cross-check the supplied author name with the email address.

    :param author: (Str) the author name to be validated.
    :return: (bool) result.
    """

    result = _check_email(author)

    return result


def _check_ccreport(ccreport):
    
    """
    A function to cross-check the supplied ccreport with the email address.

    :param ccreport: (Str) the ccreport name to be validated.
    :return: (bool) result.
    """

    result = True

    for string in ccreport.split(','):
        if not _check_email(string):
           result = False

    return result


def _check_targsrvy(array, trimester, cat_filename):

    icd_030_targsrvy_set = [
        'Ga-LRDISC', 'GA-LRHIGHLAT', 'GA-HR', 'GA-OC', 'GA-CALIB', 'STEPS',
        'SCIP-AC', 'SCIP-CYG', 'SCIP-LR', 'WA', 'WC', 'WL-WIDE', 'WL-MID',
        'WL-DEEP', 'WQ', 'ASTRO-CALIB', 'WD', 'GS', 'ING-SYSCAT'
    ]
    no_trimester_targsrvy_set = ['ASTRO-CALIB', 'WD', 'GS', 'ING-SYSCAT']

    sv_ot_regex = '^W[SV]([0-9]{4}[AB][12])-[0-9]{3}$'

    result = True

    targsrvy_set = set(array)

    for targsrvy_value in targsrvy_set:
        if targsrvy_value not in icd_030_targsrvy_set:

            if not re.match(sv_ot_regex, targsrvy_value):
                logging.error('unexpected TARGSRVY value: {}'.format(
                    targsrvy_value))
                result = False
            else:
                match = re.match(sv_ot_regex, targsrvy_value)

                trimester_in_sv_ot, = match.groups()

                if trimester_in_sv_ot != trimester:
                    # The trimester inferred from the survey name does not match the submission trimester
                    logging.warning(
                        'unexpected TARGSRVY value for TRIMESTER {}: {} (this could be a carryover survey)'.format(
                            trimester, targsrvy_value))

    if len(targsrvy_set) > 1:
        logging.error(
            'there is more than one TARGSRV value in the catalogue')
        result = False

    if result == True:

        targsrvy = list(targsrvy_set)[0]
        
        if re.match(sv_ot_regex, targsrvy_value) or targsrvy_value in no_trimester_targsrvy_set:

            expected_filename = '{}-ifu_driver_cat.fits'.format(targsrvy)

        else:

            expected_filename = '{}_{}-ifu_driver_cat.fits'.format(
                targsrvy, trimester)

        if expected_filename != os.path.basename(cat_filename):
            logging.error(
                'unexpected filename for TARGSRVY/TRIMESTE {}/{}: {}'.format(
                    targsrvy, trimester, expected_filename))
            result=False

    return result


def _check_non_empty_string(array):

    result = True

    set_of_values = set(array)

    if '' in set_of_values:
        result = False

    return result


def _check_targid(array):

    result = _check_non_empty_string(array)

    if result is False:
        logging.error('TARGID value should not be populated by empty strings')

    return result


def _check_targname(array):

    result = _check_non_empty_string(array)

    if result is False:
        logging.error('TARGNAME value should not be populated by empty strings')

    return result


def _check_non_null_value(array):

    result = True

    for value in array:
        if np.isnan(value):
            result = False
            break

    return result


def _get_tl(hdu, colname, minmax):

    assert minmax in ['min', 'max']

    if colname in hdu.data.columns.names:

        index = hdu.data.columns.names.index(colname)

        kwd = 'TL{}{}'.format(minmax.upper(), index + 1)

        if kwd in hdu.header:
            value = hdu.header[kwd]
        else:
            value = None

    else:
        value = None

    return value


def _check_in_tl_range(hdu, colname):

    result = True

    min_value = _get_tl(hdu, colname, 'min')
    max_value = _get_tl(hdu, colname, 'max')
    
    if min_value is None:
        logging.warning(
            'TLMIN value is not available for column {}'.format(colname))
    
    if max_value is None:
        logging.warning(
            'TLMAX value is not available for column {}'.format(colname))

    for value in hdu.data[colname]:

        if not np.isnan(value):

            if min_value is not None:
                if not (value >= min_value):
                    result = False
                    break

            if max_value is not None:
                if not (value <= max_value):
                    result = False
                    break

    return result


def _check_targprio(hdu, progtemp_dict):

    result = True

    null_result = _check_non_null_value(hdu.data['TARGPRIO'])

    if null_result is False:
        logging.error('TARGPRIO values contain NaN')
        result = False

    range_result = _check_in_tl_range(hdu, 'TARGPRIO')

    if range_result is False:
        logging.error('TARGPRIO values out of range')
        result = False

    for targprio, progtemp in zip(hdu.data['TARGPRIO'], hdu.data['PROGTEMP']):

        obsmode = get_obsmode_from_progtemp(progtemp,
                                            progtemp_dict=progtemp_dict)

        if obsmode == 'LIFU':
            if targprio != 10.0:
                logging.error('TARGPRIO value must be 10.0 for LIFU mode')
                result = False

    return result


def _check_progtemp(array, progtemp_dict, forbidden_dict):

    result = True

    progtemp_set = set(array)

    for progtemp_value in progtemp_set:

        try:
            spectrograph_dict = get_progtemp_info(progtemp_value,
                                                  progtemp_dict=progtemp_dict)
        except:
            logging.error('unexpected PROGTEMP value: {}'.format(value))
            result = False

        for survey_type in forbidden_dict.keys():
            for regex in forbidden_dict[survey_type]:
                if re.match(regex, progtemp_value):
                    logging.warning(
                        'PROGTEMP value forbidden for {} surveys: {}'.format(
                            survey_type.upper(), progtemp_value))

    return result


def _check_obstemp(array, obstemp_dict):

    result = True

    obstemp_set = set(array)

    for obstemp_value in obstemp_set:

        try:
            obsconstraints_dict = get_obstemp_info(obstemp_value,
                                                   obstemp_dict=obstemp_dict)
        except:
            logging.error('unexpected OBSTEMP value: {}'.format(value))
            result = False

    return result


def _check_gaia_ra(hdu):

    result = True

    null_result = _check_non_null_value(hdu.data['GAIA_RA'])

    if null_result is False:
        logging.error('GAIA_RA values contain NaN')
        result = False

    range_result = _check_in_tl_range(hdu, 'GAIA_RA')

    if range_result is False:
        logging.error('GAIA_RA values out of range')
        result = False

    return result


def _check_gaia_dec(hdu):

    result = True

    null_result = _check_non_null_value(hdu.data['GAIA_DEC'])

    if null_result is False:
        logging.error('GAIA_DEC values contain NaN')
        result = False

    range_result = _check_in_tl_range(hdu, 'GAIA_DEC')

    if range_result is False:
        logging.error('GAIA_DEC values out of range')
        result = False

    return result


def _check_gaia_epoch(hdu):

    result = True

    null_result = _check_non_null_value(hdu.data['GAIA_EPOCH'])

    if null_result is False:
        logging.error('GAIA_EPOCH values contain NaN')
        result = False

    range_result = _check_in_tl_range(hdu, 'GAIA_EPOCH')

    if range_result is False:
        logging.error('GAIA_EPOCH values out of range')
        result = False

    return result


def _check_gaia_pm(hdu):

    result = True

    pmra_null_result = _check_non_null_value(hdu.data['GAIA_PMRA'])

    if pmra_null_result is False:
        logging.error('GAIA_PMRA values contain NaN')
        result = False

    pmdec_null_result = _check_non_null_value(hdu.data['GAIA_PMDEC'])

    if pmdec_null_result is False:
        logging.error('GAIA_PMDEC values contain NaN')
        result = False
    
    return result


def _check_gaia_paral(hdu):

    result = True

    null_result = _check_non_null_value(hdu.data['GAIA_PARAL'])

    if null_result is False:
        logging.error('GAIA_PARAL values contain NaN')
        result = False
    
    return result


def _check_ifu_pa_request(hdu, progtemp_dict):

    result = True

    range_result = _check_in_tl_range(hdu, 'IFU_PA_REQUEST')

    if range_result is False:
        logging.error('IFU_PA_REQUEST values out of range')
        result = False

    ifu_pa_request_array = hdu.data['IFU_PA_REQUEST']
    progtemp_array = hdu.data['PROGTEMP']

    for i, (ifu_pa_request, progtemp) in enumerate(
            zip(ifu_pa_request_array, progtemp_array)):

        if (not np.isnan(ifu_pa_request)):

            obsmode = get_obsmode_from_progtemp(progtemp, progtemp_dict)

            if obsmode != 'LIFU':

                logging.error(
                    'unexpected non-null IFU_PA_REQUEST for non-LIFU PROGTEMP' +
                    ' in row {}'.format(i + 1))
                result = False

    return result


def _check_ifu_dither(hdu, progtemp_dict):

    result = True

    range_result = _check_in_tl_range(hdu, 'IFU_DITHER')

    if range_result is False:
        logging.error('IFU_DITHER values out of range')
        result = False

    ifu_dither_array = hdu.data['IFU_DITHER']
    progtemp_array = hdu.data['PROGTEMP']

    for i, (ifu_dither, progtemp) in enumerate(
            zip(ifu_dither_array, progtemp_array)):


        obsmode = get_obsmode_from_progtemp(progtemp, progtemp_dict)

        if (obsmode == 'LIFU'):
            if (ifu_dither not in [-1, 0, -3, 3, 4, 5, 6]):
                logging.error(
                    'unexpected IFU_DITHER for LIFU PROGTEMP in row {}'.format(
                        i + 1))
                result = False
        elif (obsmode == 'mIFU'):
            if (ifu_dither not in [0, -3, 3, 4, 5, 6]):
                logging.error(
                    'unexpected IFU_DITHER for mIFU PROGTEMP in row {}'.format(
                        i + 1))
                result = False
        elif (ifu_dither != 0):
            logging.error(
                'unexpected IFU_DITHER and PROGTEMP in row {}'.format(i + 1))
            result = False

    return result

def _nan_to_none(value):

    if np.isnan(value):
        result = None
    else:
        result = value

    return result


def _check_repeated_objects(hdu, offset_tol_arcsec=10.0, include_progtemp=True):

    result = True

    # Create a dictionary with TARGNAME/TARGID as keys and empty lists as values

    targname_array = hdu.data['TARGNAME']
    targid_array = hdu.data['TARGID']
    if include_progtemp:
        progtemp_array = hdu.data['PROGTEMP']

    if include_progtemp:
        logging.warning("Including PROGTEMP in uniqueness key")

        object_dict = {(targname, targid, progtemp): []
                       for targname, targid, progtemp in zip(targname_array, targid_array, progtemp_array)}

    else:
        object_dict = {(targname, targid): []
                       for targname, targid in zip(targname_array, targid_array)}

    # Populate the lists with the information of the dithering pattern and
    # coordinates
    
    for i in range(len(hdu.data)):

        targname = hdu.data['TARGNAME'][i]
        targid = hdu.data['TARGID'][i]

        ifu_dither = hdu.data['IFU_DITHER'][i]

        gaia_ra = hdu.data['GAIA_RA'][i]
        gaia_dec = hdu.data['GAIA_DEC'][i]
        gaia_epoch = hdu.data['GAIA_EPOCH'][i]
        gaia_pmra = hdu.data['GAIA_PMRA'][i]
        gaia_pmdec = hdu.data['GAIA_PMDEC'][i]
        gaia_paral = hdu.data['GAIA_PARAL'][i]
        if include_progtemp:
            progtemp = hdu.data['PROGTEMP'][i]

            object_dict[(targname, targid, progtemp)].append(
                (ifu_dither,
                 (gaia_ra, gaia_dec, _nan_to_none(gaia_epoch),
                  _nan_to_none(gaia_pmra), _nan_to_none(gaia_pmdec),
                  _nan_to_none(gaia_paral),
                  progtemp))
            )
        else:
            object_dict[(targname, targid)].append(
                (ifu_dither,
                 (gaia_ra, gaia_dec, _nan_to_none(gaia_epoch),
                  _nan_to_none(gaia_pmra), _nan_to_none(gaia_pmdec),
                  _nan_to_none(gaia_paral)))
            )

            
    # Check the uniqueness of the coordinates
    for key in object_dict.keys():

        if include_progtemp:
            targname, targid, progtemp = key        
        else:
            targname, targid = key        
        list_in_value = object_dict[key]

        ifu_dither_list = [elem[0] for elem in list_in_value]
        coord_list = [elem[1] for elem in list_in_value]

        # Case 1: We do not have custom dithers
        if (-1 not in ifu_dither_list):
            # all the coordinates shoule be the same:
            if len(set(coord_list)) != 1:
                result = False
                logging.error(
                    'Multiple sets of coordinates in non-custom dither(s) of (TARGNAME:TARGID) ({}:{}) :'.format(
                        targname, targid))
                for _coord in set(coord_list):
                    logging.error(str(_coord))
                logging.error(
                    "Coordinates must be the same for a provided TARGID and TARGNAME combination")
                logging.error(
                    "If you wish to specify a new target to observe, please use a different TARGID or TARGNAME")

            if len(set(ifu_dither_list)) != len(ifu_dither_list):
                result = False
                logging.error(
                    'Target defined by (TARGNAME:TARGID) ({}:{}) has been specified multiple times for a given dither pattern :'.format(
                        targname, targid))
                _counts = Counter(ifu_dither_list)

                for _dither in _counts:
                    if _counts[_dither] > 1:
                        logging.error("Dither pattern {} - {} times".format(_dither,_counts[_dither]))
                logging.error("If you wish to observe an object across multiple fields, either re-group that target at input level")
                logging.error("with a new TARGID+TARGNAME (alonside other similar targets), or input the same field multiple times to configure")
                
        else:

            # If we have a mix of non-custom and custom dithers, take the non-
            # custom as reference. If we have only custom dithers, take the
            # first position as reference
            
            non_custom_coord_list = [elem[1] for elem in list_in_value
                                     if elem[0] != -1]

            non_custom_coord_set = set(non_custom_coord_list)

            # All the non-custom dithers should have the same coordinates
            
            if len(non_custom_coord_set) > 1:
                logging.error(
                    'unexpected change of coordinates for {}:{}'.format(
                        targname, targid))
                result = False
                continue
            elif len(non_custom_coord_set) == 1:
                ref_coord = non_custom_coord_list[0]
            else:
                ref_coord = coord_list[0]

            # Check that the offsets are not far from the reference
            
            for ifu_dither, coord in list_in_value:
                
                if ifu_dither == -1:

                    ra_offset = np.abs((coord[0] - ref_coord[0]) *
                                       np.cos(np.deg2rad(ref_coord[1]))) * 3600
                    dec_offset = np.abs(coord[1] - ref_coord[1]) * 3600

                    if ((ra_offset > offset_tol_arcsec) or
                        (dec_offset > offset_tol_arcsec)):
                        logging.error(
                            'too big offset for custom dither in {}:{}'.format(
                                targname, targid))
                        result = False
                        continue

                    for i in range(2, 6):
                        if ((coord[i] != ref_coord[i]) and not
                            (np.isnan(coord[i]) and np.isnan(ref_coord[i]))):
                            logging.error(
                                'unexpected change of GAIA_EPOCH/GAIA_PMRA/'
                                'GAIA_PMDEC/GAIA_PARAL for {}:{}'.format(
                                    targname, targid))
                            result = False
                            continue

    return result


def _check_progtemp_ifu_dither_compatibility(hdu, progtemp_dict):

    result = True

    progtemp_array= hdu.data['PROGTEMP']
    ifu_dither_array = hdu.data['IFU_DITHER']

    # Check the non-custom dithers

    for i, (progtemp, ifu_dither) in enumerate(
            zip(progtemp_array, ifu_dither_array)):

        # doesn't correctly deal with ifu_dither==0
        #if np.abs(ifu_dither) != 1:
        if np.abs(ifu_dither) >= 3:

            num_dither = np.abs(ifu_dither)

            spectrograph_dict = get_progtemp_info(progtemp,
                                                  progtemp_dict=progtemp_dict)

            num_exp_r = spectrograph_dict['red_num_exposures']
            num_exp_b = spectrograph_dict['blue_num_exposures']
 
            if ((num_exp_r != num_dither) or 
                (num_exp_b != num_dither)):
                logging.error(
                    'PROGTEMP ({}/{} exposures in red/blue) and IFU_DITHER ({} exposures) are not compatible in row {}'.format(
                        num_exp_r,num_exp_b,num_dither,i + 1))
                result = False
 
            if not (num_exp_r == num_exp_b):
                logging.error(
                    'PROGTEMP-defined number of red and blue exposures are not equal in row {}'.format(
                        i + 1))
                result = False

    # Check the custom dithers

    targname_array = hdu.data['TARGNAME']
    targid_array = hdu.data['TARGID']
    obstemp_array = hdu.data['OBSTEMP']

    custom_dict = OrderedDict()

    for i, (ifu_dither, targname, targid, progtemp, obstemp) in enumerate(
            zip(ifu_dither_array,
                targname_array, targid_array, progtemp_array, obstemp_array)):

        if ifu_dither == -1:

            key = (targname, targid, progtemp, obstemp)

            if key in custom_dict:
                custom_dict[key] += 1
            else:
                custom_dict[key] = 1

    for key in custom_dict.keys():

        targname, targid, progtemp, obstemp = key

        num_dither = custom_dict[key]

        spectrograph_dict = get_progtemp_info(progtemp,
                                              progtemp_dict=progtemp_dict)

        num_exp_r = spectrograph_dict['red_num_exposures']
        num_exp_b = spectrograph_dict['blue_num_exposures']

        if ((num_exp_r % num_dither != 0) or 
            (num_exp_b % num_dither != 0)):
            logging.error(
                'PROGTEMP and custom dither {}:{}:{}:{} are not compatible'.format(
                    targname, targid, progtemp, obstemp)
            )
            result = False
 
        if not ((num_exp_r % num_exp_b == 0) or
                (num_exp_b % num_exp_r == 0)):
            logging.error(
                'PROGTEMP and custom dither {}:{}:{}:{} are not compatible'.format(
                    targname, targid, progtemp, obstemp)
            )
            result = False

    return result


def _check_lifu_rotlim_visibility(hdu, ha_range=2.0, cache=False):
    """
    Check the observability of the pointing given a canonical hour angle range
    """

    # Get the obslimit assets first....
    # 1. configure.xml:
    xml_configure = get_configure_cfg(cache=cache)

    # 2. point to where obs_limits.py is:
    obslim_py = get_obslimits_tool(cache=cache)

    # 3. Add the parent directory of the py file to system path 
    # sys.path.append(os.path.dirname(obslim_py))
    sys.path = [os.path.dirname(obslim_py)] + sys.path
    

    # 4. Now import the dir containing the .py file
    from obs_limits import ObsLimits

    result = True

    progtemp_array = hdu.data['PROGTEMP']
    dec_array = hdu.data['GAIA_DEC']
    pa_array = hdu.data['IFU_PA_REQUEST']
    ha_array = hdu.data['HA_RANGE']

    # Filter for only LIFU pointings
    mask = np.array([progtemp[0] in ['4','5','6'] for progtemp in progtemp_array])
    sel_index = np.arange(np.shape(mask)[0])
    
    dec_lifu = dec_array[mask]
    pa_lifu = pa_array[mask]
    ha_lifu = ha_array[mask]
    sel_index = sel_index[mask]

    pa_check = ObsLimits(cfg=xml_configure)
    

    for  _dec, _pa, _ha, _index in zip(dec_lifu,pa_lifu,ha_lifu,sel_index):
        # pa_valid = pa_check.pa_check(_ha*15.,_dec,_pa,verbose=False)

        pa_valid, lim_data = pa_check.limit_check(_ha*15.,_dec,_pa,detail=True,verbose=False)

        if not pa_valid:
            result = False
            row = hdu.data[_index]
            pa_ranges = pa_check.pa_range(ha_range*15.,_dec,verbose=False,explicit=True)
            logging.error("{} LIFU entry at index {} (TARGID: {}) has Dec ({:.3f}) and PA ({:.1f}) that exceeds rotator or horizon limits for requested HA range (+/-{:.1f}). ".format(
                row['TARGPROG'],_index,row['TARGID'],_dec,_pa,_ha))

            dha = 0.2
            logging.error("Viable PA range across the HA range requested:")
            
            for _ha in np.arange(-1.0*abs(_ha),abs(_ha)+(0.5*dha),dha):
                if abs(_ha-0.0) < 1e-6:
                    # skip
                    continue
                _minmax = pa_check.pa_range(_ha*15.,_dec,verbose=False,explicit=True)
                range_plural = ''
                if len(_minmax) > 1:
                    range_plural = 's'
                logging.error("@HA={:1.3f}: PA range{} availble: {}".format(_ha,range_plural,_minmax))

    if not result:
        logging.error("")
        logging.error("One or more requested LIFU pointings cannot be observed due to WEAVE rotator limits")
        logging.error("This arises for pointings at a particular Declination and PA under a {}hr HA angle range".format(ha_range))
        logging.error("Please use the obs_limits.py tool (local copy: {}) to check the legal PA for your pointing".format(obslim_py))
        logging.error("Note that this might impact selection of guide stars for the chosen PA")
        logging.error("We advise users to identify the viable PA range and select a PA at the mid point of this")
        
    return result
        

def check_ifu_driver_cat(cat_filename, template=None, check_vs_template=True,
                         progtemp_file=None, obstemp_file=None, ha_range=2.0,
                         remove_on_fail=False, cache=False):
    """
    Check the contents of an IFU driver catalogue.

    Parameters
    ----------
    cat_filename : str
        A FITS file with an IFU driver catalogue.
    template : str, optional
        A FITS file containing an IFU driver template.
    check_vs_template: bool, optional
        An option to indicate whether the file should be checked against a
        template or not.
    progtemp_file : str, optional
        A progtemp.dat file with the definition of PROGTEMP.
    obstemp_file : str, optional
        A obstemp.dat file with the definition of OBSTEMP.
    ha_range : float, optional
        Set the hour angle limits for determining LIFU observability
    remove_on_fail : bool, optional
        If the catalogue fails, delete it
    cache : bool, optional
        Do not download catalogue files if they aready exist

    Returns
    -------
    result : bool
        True if the file passes all the checks, otherwise False.
    """
    
    result = True

    # Get the dictionaries to interpret PROGTEMP and OBSTEMP, and their
    # DATAMVER values

    progtemp_datamver, progtemp_dict, forbidden_dict = get_progtemp_dict(
        filename=progtemp_file, assert_orb=True, cache=cache)

    obstemp_datamver, obstemp_dict = get_obstemp_dict(
        filename=obstemp_file, cache=cache)

    # Check the file versus a template (if requested)
    
    if check_vs_template is True:

        if template is not None:

            if not _check_versus_template(cat_filename, template):
                logging.error(
                    '{} does not match the template {}'.format(
                        cat_filename, template))
                result = False

        else:
            logging.error(
                'a template for making the checks has not been provided')
            result = False

    # Open the catalogue file

    with fits.open(cat_filename) as hdu_list:
        
        # Check DATAMVER of the IFU driver cat, PROGTEMP file and OBSTEMP file
        # are consistent

        datamver = hdu_list[0].header['DATAMVER']

        if datamver != progtemp_datamver:
            logging.critical(
                'DATAMVER mismatch ({} != {}) for PROGTEMP file: '.format(
                    datamver, progtemp_datamver) +
                'Stop unless you are sure!')
            raise SystemExit(2)

        if datamver != obstemp_datamver:
            logging.critical(
                'DATAMVER mismatch ({} != {}) for OBSTEMP file: '.format(
                    datamver, obstemp_datamver) +
                'Stop unless you are sure!')
            raise SystemExit(2)

        # Check the values of its primary header

        trimester = hdu_list[0].header['TRIMESTE']
        if not _check_trimester(trimester):
            logging.error(
                'invalid value in the TRIMESTE keyword: {}'.format(trimester))
            result = False

        verbose = hdu_list[0].header['VERBOSE']
        if not _check_verbose(verbose):
            logging.error(
                'invalid value in the VERBOSE keyword: {}'.format(verbose))
            result = False

        author = hdu_list[0].header['CAT_MAIL']
        if not _check_author(author):
            logging.error(
                'invalid email address in the CAT_MAIL keyword: {}'.format(
                    author))
            result = False

        ccreport = hdu_list[0].header['CAT_CC']
        if (len(ccreport)) and (not _check_ccreport(ccreport)):
            logging.error(
                'invalid email addresses in the CAT_CC keyword: {}'.format(
                    ccreport))
            result = False

        # Check the data in the table of its first extension

        hdu1 = hdu_list[1]
        data = hdu_list[1].data

        if not _check_targsrvy(data['TARGSRVY'], trimester, cat_filename):
            result = False

        # Nothing to check in TARGPROG

        if not _check_targid(data['TARGID']):
            result = False

        if not _check_targname(data['TARGNAME']):
            result = False

        if not _check_targprio(hdu1, progtemp_dict):
            result = False

        if not _check_progtemp(data['PROGTEMP'], progtemp_dict, forbidden_dict):
            result = False

        if not _check_obstemp(data['OBSTEMP'], obstemp_dict):
            result = False

        if not _check_gaia_ra(hdu1):
            result = False

        if not _check_gaia_dec(hdu1):
            result = False

        if not _check_gaia_epoch(hdu1):
            result = False

        if not _check_gaia_pm(hdu1):
            result = False

        if not _check_gaia_paral(hdu1):
            result = False

        if not _check_ifu_pa_request(hdu1, progtemp_dict):
            result = False

        if not _check_ifu_dither(hdu1, progtemp_dict):
            result = False

        if not _check_lifu_rotlim_visibility(hdu1, ha_range=ha_range, cache=cache):
            result = False
            
        if not _check_repeated_objects(hdu1):
            result = False

        if not _check_progtemp_ifu_dither_compatibility(hdu1, progtemp_dict):
            result = False

    # Write a message

    if not cat_filename.endswith('-ifu_driver_cat.fits'):
        logging.warning(
            '{} does not follow the filename convention'.format(cat_filename) +
            ': it should be the name of the catalogue FITS file submitted to ' +
            'WASP with the suffix "-ifu_driver_cat.fits". Please, follow this '+
            'convention unless you know what you are doing and its impact.' +
            'Do not request for help for not following the convention, ' +
            'unless you have read the whole documentation from top to bottom ' +
            'and you have strong reasons for doing it.')

    # Write a message

    if result == True:
        logging.info('{} has passed the checks'.format(cat_filename))
    else:
        logging.error('{} has NOT passed the checks'.format(cat_filename))
        if remove_on_fail:
            logging.warning("Removing file {}".format(cat_filename))
            os.unlink(cat_filename)
    return result

    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(
        description='Check the contents of an IFU driver catalogue')
    
    parser.add_argument('catalogue',
                        help='a FITS file with an IFU driver catalogue')
    
    parser.add_argument('--template',
                        default=os.path.join('aux', 'ifu_driver_template.fits'),
                        help='a FITS file containing an IFU driver template')

    parser.add_argument('--no_check_vs_template', dest='check_vs_template',
                        action='store_false',
                        help=
                        'skip the check of the catalogue versus the template')

    parser.add_argument('--progtemp_file',
                        default=os.path.join('aux', 'progtemp.dat'),
                        help="""a progtemp.dat file with the definition of
                        PROGTEMP""")

    parser.add_argument('--obstemp_file',
                        default=os.path.join('aux', 'obstemp.dat'),
                        help="""a obstemp.dat file with the definition of
                        OBSTEMP""")

    parser.add_argument('--ha_range', default=2.0, type=float,
                        help="""The hour angle limits for determining 
                        LIFU observability based on provided Dec and PA
                        """)

    parser.add_argument('--remove_on_fail', action='store_true',
                        help="""Delete the input catalogue if it does not
                        pass the checks in this module""")
    
    parser.add_argument('--cache', action='store_true',
                        help="""Do not download observing limit assets if they
                        already exist""")
    
    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')
    
    args = parser.parse_args()
    
    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    if not os.path.exists(args.progtemp_file):
        logging.info('Downloading the progtemp file')
        get_progtemp_file(file_path=args.progtemp_file)

    if not os.path.exists(args.obstemp_file):
        logging.info('Downloading the obstemp file')
        get_obstemp_file(file_path=args.obstemp_file)
    
    check_ifu_driver_cat(args.catalogue, template=args.template,
                         check_vs_template=args.check_vs_template,
                         progtemp_file=args.progtemp_file,
                         obstemp_file=args.obstemp_file,
                         remove_on_fail=args.remove_on_fail,
                         cache=args.cache)

