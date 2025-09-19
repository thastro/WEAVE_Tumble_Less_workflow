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


from collections import OrderedDict
import xml.dom.minidom as _minidom

import numpy as _np


def _get_lookup(namespace='fits', post_configure=True):

    assert namespace in ['fits', 'xml']

    lookup = {}

    # Set the dictionary for the 'fits' namespace

    lookup['target:cname'] = 'CNAME'
    lookup['target:targsrvy'] = 'TARGSRVY'
    lookup['target:targprog'] = 'TARGPROG'
    lookup['target:targcat'] = 'TARGCAT'
    lookup['target:targid'] = 'TARGID'
    lookup['target:targname'] = 'TARGNAME'
    lookup['target:targprio'] = 'TARGPRIO'
    lookup['target:targuse'] = 'TARGUSE'
    lookup['target:targclass'] = 'TARGCLASS'
    lookup['observation:progtemp'] = 'PROGTEMP'
    lookup['observation:obstemp'] = 'OBSTEMP'
    lookup[''] = 'GAIA_ID'
    lookup[''] = 'GAIA_DR'
    lookup['target:targra'] = 'GAIA_RA'
    lookup['target:targdec'] = 'GAIA_DEC'
    lookup['target:targepoch'] = 'GAIA_EPOCH'
    lookup['target:targpmra'] = 'GAIA_PMRA'
    lookup[''] = 'GAIA_PMRA_ERR'
    lookup['target:targpmdec'] = 'GAIA_PMDEC'
    lookup[''] = 'GAIA_PMDEC_ERR'
    lookup['target:targparal'] = 'GAIA_PARAL'
    lookup[''] = 'GAIA_PARAL_ERR'
    lookup[''] = 'HEALPIX'
    lookup['target:ifu_spaxel'] = 'IFU_SPAXEL'
    if post_configure is True:
        lookup['target:ifu_pa'] = 'IFU_PA'
    else:
        lookup['observation:pa'] = 'IFU_PA'
    lookup['dithering:apply_dither'] = 'IFU_DITHER'
    lookup['photometry:mag_g'] = 'MAG_G'
    lookup['photometry:emag_g'] = 'MAG_G_ERR'
    lookup['photometry:mag_r'] = 'MAG_R'
    lookup['photometry:emag_r'] = 'MAG_R_ERR'
    lookup['photometry:mag_i'] = 'MAG_I'
    lookup['photometry:emag_i'] = 'MAG_I_ERR'
    lookup['photometry:mag_gg'] = 'GAIA_MAG_G'
    lookup['photometry:emag_gg'] = 'GAIA_MAG_G_ERR'
    lookup['photometry:mag_bp'] = 'GAIA_MAG_BP'
    lookup['photometry:emag_bp'] = 'GAIA_MAG_BP_ERR'
    lookup['photometry:mag_rp'] = 'GAIA_MAG_RP'
    lookup['photometry:emag_rp'] = 'GAIA_MAG_RP_ERR'

    lookup.pop('')
    
    # Use the dictionary for 'fits' namespace to set dictionary for 'xml'
    
    if namespace == 'xml':
        
        # Remove some keys from the dictionary
        
        keys_to_remove_for_xml = [
            'observation:progtemp', 'observation:obstemp',
            'dithering:apply_dither'
        ]
        
        for key in keys_to_remove_for_xml:
            lookup.pop(key)
        
        # Set the values of the dictionary from the attrib name
        
        for key in lookup.keys():
            lookup[key] = key.split(':')[1]
        
        # Set additional keys
        
        lookup['target:automatic'] = 'automatic'
        lookup['target:configid'] = 'configid'
        lookup['target:fibreid'] = 'fibreid'
        lookup['target:targx'] = 'targx'
        lookup['target:targy'] = 'targy'
        
        lookup['simulation:filterid'] = 'sim_filterid'
        lookup['simulation:fwhm'] = 'sim_fwhm'
        lookup['simulation:mag'] = 'sim_mag'
        lookup['simulation:redshift'] = 'sim_redshift'
        lookup['simulation:template'] = 'sim_template'
        lookup['simulation:velocity'] = 'sim_velocity'

    return lookup


def _get_formats(post_configure=True):

    formats = {}

    formats['target:targprio'] = float
    formats['target:targra'] = float
    formats['target:targdec'] = float
    formats['target:targepoch'] = float
    formats['target:targpmra'] = float
    formats[''] = float # 'GAIA_PMRA_ERR'
    formats['target:targpmdec'] = float
    formats[''] = float # 'GAIA_PMDEC_ERR'
    formats['target:targparal'] = float
    formats[''] = float # 'GAIA_PARAL_ERR'
    formats[''] = int # 'HEALPIX'
    if post_configure is True:
        formats['target:ifu_pa'] = float
    else:
        formats['observation:pa'] = float
    formats['dithering:apply_dither'] = int
    formats['photometry:mag_g'] = float
    formats['photometry:emag_g'] = float
    formats['photometry:mag_r'] = float
    formats['photometry:emag_r'] = float
    formats['photometry:mag_i'] = float
    formats['photometry:emag_i'] = float
    formats['photometry:mag_gg'] = float
    formats['photometry:emag_gg'] = float
    formats['photometry:mag_bp'] = float
    formats['photometry:emag_bp'] = float
    formats['photometry:mag_rp'] = float
    formats['photometry:emag_rp'] = float
    
    formats['target:automatic'] = int
    formats['target:configid'] = int
    formats['target:fibreid'] = int
    formats['target:targx'] = float
    formats['target:targy'] = float
    formats['simulation:fwhm'] = float
    formats['simulation:mag'] = float
    formats['simulation:redshift'] = float
    formats['simulation:velocity'] = float

    formats.pop('')

    return formats


def _get_xml_data(xml_filename):

    xml_data = {}

    dom = _minidom.parse(xml_filename)

    observation = dom.getElementsByTagName('observation')[0]
    dithering = dom.getElementsByTagName('dithering')[0]
    field = dom.getElementsByTagName('field')
    target = dom.getElementsByTagName('target')

    xml_data['observation'] = observation
    xml_data['dithering'] = dithering
    xml_data['field'] = field
    xml_data['target'] = target

    return xml_data


def _get_value_from_xml_data(xml_data, target, key, formats):

    # Get the column, element and attribute name
    
    element_name, attribute_name = key.split(':')
    
    # It the name of the element is 'target', choose as the element
    # the current target; otherwise, get the element from the
    # xml_data
    
    if element_name == 'target':
        element = target
    elif element_name == 'photometry':
        
        photometry_element_list = target.getElementsByTagName('photometry')
        
        # If the photometry element exist, select it and continue
        
        if len(photometry_element_list) > 0:
            element = photometry_element_list[0]
        
        # If there is not a photometry element (like for the sky fibres), the
        # key should correspond to a float in the formats dictionary, so a NaN
        # will be returned
        
        else:
            
            assert (formats[key] == float)
            
            formatted_value = _np.nan
            
            return formatted_value
        
    elif element_name == 'simulation':
        
        simulation_element_list = target.getElementsByTagName('simulation')
        
        # If the simulation element exist, select it and continue
        
        if len(simulation_element_list) > 0:
            element = simulation_element_list[0]
        
        # If there is not a simulation element (like for the sky fibres), the
        # key should correspond to a float in the formats dictionary, so a NaN
        # will be returned
        
        else:
            
            formatted_value = None
            
            return formatted_value
        
    else:
        element = xml_data[element_name]
    
    # Get the raw value from the XML

    raw_value = str(element.getAttribute(attribute_name))

    if key in formats.keys():
        
        format_func = formats[key]
        
        if raw_value not in ['', '%%%']:
            formatted_value = format_func(raw_value)
        elif ((format_func == float) and
              ((raw_value == '') or (raw_value == '%%%'))):
            formatted_value = _np.nan
        elif ((format_func == int) and
              ((raw_value == '') or (raw_value == '%%%'))):
            formatted_value = _np.nan
        else:
            raise ValueError(
                      'raw value {}, format {}'.format(raw_value, format_func))
    else:
        formatted_value = raw_value
    
    return formatted_value


def _get_data_from_targuse_per_field_of_one_xml(
        xml_filename, targuse_list=['T', 'S', 'R'], replace_triple_percent=True,
        post_configure=True, only_allocated=False, sort_targets=False,
        namespace='fits'):

    # Get dictionaries with the lookup information and the formats

    lookup = _get_lookup(namespace=namespace, post_configure=post_configure)
    formats = _get_formats(post_configure=post_configure)

    # Get the data from the XML file

    xml_data = _get_xml_data(xml_filename)

    # Create a dictionary with the desired keywords and empty lists in order to
    # get ready to store them the data

    data_dict_list = [{lookup[key]: [] for key in lookup.keys()}
                      for i in range(len(xml_data['field']))]
    
    # For each field in the XML
    
    for i, field in enumerate(xml_data['field']):
        
        # Get an ordered dictionary with the targets
        # It will be useful if key decide to sort the result or not
        
        target_ord_dict = OrderedDict()
        
        for j, target in enumerate(field.getElementsByTagName('target')):
        
            # The key will have 4 components:
            #  - ifu_bundle: 'L' (for LIFU) or 'm01' ... 'm20' (for mIFU)
            #  - automatic: 0 or 1
            #  - ifu_spaxel
            #  - j (the number of target, to prevent duplicities in the key)
            
            # This will allow us to short the fibres by bundles, having the
            # first one the central fibre
        
            if 'ifu_spaxel' in target.attributes.keys():
                ifu_spaxel = str(target.getAttribute('ifu_spaxel'))
                
                if len(ifu_spaxel) == 3:
                    ifu_bundle = 'L'
                elif len(ifu_spaxel) == 6:
                    ifu_bundle = ifu_spaxel[:3]
                else:
                    ifu_bundle = ifu_spaxel
            else:
                ifu_spaxel = ''
                ifu_bundle = ''
            
            if 'automatic' in target.attributes.keys():
                automatic = int(target.getAttribute('automatic'))
            else:
                automatic = 0
            
            target_key = (ifu_bundle, automatic, ifu_spaxel, j)
            
            target_ord_dict[target_key] = target
        
        # Get the list of keywords of the ordered dict, and sort them if
        # requested
        
        target_key_list = list(target_ord_dict.keys())
        
        if sort_targets == True:
            target_key_list.sort()

        # For each target in the XML

        for target_key in target_key_list:
        
            target = target_ord_dict[target_key]

            # Skip the target if it is not in the targuse list

            if targuse_list is not None:
                if str(target.getAttribute('targuse')) not in targuse_list:
                    continue
            
            # If we want to skip the non-allocated targets, skip the target
            # if it is the case
            
            if only_allocated == True:
                if 'fibreid' not in target.attributes.keys():
                    continue

            # For each key in the lookup dictionary
            
            for key in lookup.keys():
                
                # Get its column name and its value and save them
                
                col_name = lookup[key]
                
                value = _get_value_from_xml_data(
                            xml_data, target, key, formats)
                
                # Replace the triple percent if requested
                
                if (replace_triple_percent is True) and (value == '%%%'):
                    value = ''

                # Save the value in the corresponding column
                
                data_dict_list[i][col_name].append(value)

    return data_dict_list


def _get_spa_data_from_targuse(xml_filename_list, targuse_list=['T', 'S', 'R'],
                               replace_triple_percent=True,
                               post_configure=True, only_allocated=False,
                               sort_targets=False):
                               
    # Set the value of the namespace to be used in the calls of _get_lookup and
    # _get_data_from_targuse_per_field_of_one_xml
    
    namespace = 'fits'

    # Get dictionaries with the lookup

    lookup = _get_lookup(namespace=namespace, post_configure=post_configure)

    # Create a dictionary with the desired keywords and empty lists in order to
    # get ready to store them the data

    data_dict = {lookup[key]: [] for key in lookup.keys()}

    # For each XML file
    
    for xml_filename in xml_filename_list:

        data_dict_list = _get_data_from_targuse_per_field_of_one_xml(
            xml_filename, targuse_list=targuse_list,
            replace_triple_percent=replace_triple_percent,
            post_configure=post_configure, only_allocated=only_allocated,
            sort_targets=sort_targets, namespace=namespace)
        
        for field_data_dict in data_dict_list:
            for col_name in data_dict.keys():
                data_dict[col_name].extend(field_data_dict[col_name])

    return data_dict


def get_spa_data_of_target_fibres_from_xmls(
        xml_filename_list, replace_triple_percent=True, post_configure=True,
        only_allocated=False, sort_targets=False):
    """
    Get SPA data of the target fibres contained in a list of XML files.

    Parameters
    ----------
    xml_filename_list : list of str
        A list of configure XML files.
    replace_triple_percent : bool, optional
        It will replace the values of '%%%' by '' also for data of type string.
    post_configure : bool, optional
        An option to indicate whether the input XMLs has been processed by
        configure or not.
    only_allocated : bool, optional
        It will include only information from allocated fibres.
    sort_targets : bool, optional
        Sort the data inside each field grouping by IFU bundle.

    Returns
    -------
    data_dict : dict
        A dictionary with the data. Its keys are the name of the SPA columns
        which are expected to be potencially in the configure XML files, while
        its values are lists containing these data.
    """
    data_dict = \
        _get_spa_data_from_targuse(
            xml_filename_list, targuse_list=['T'],
            replace_triple_percent=replace_triple_percent,
            post_configure=post_configure, only_allocated=only_allocated,
            sort_targets=sort_targets)

    return data_dict


def get_spa_data_of_target_random_and_sky_fibres_from_xmls(
        xml_filename_list, replace_triple_percent=True, post_configure=True,
        only_allocated=True, sort_targets=False):
    """
    Get SPA data of the T/R/S fibres contained in a list of XML files.

    Parameters
    ----------
    xml_filename_list : list of str
        A list of configure XML files.
    replace_triple_percent : bool, optional
        It will replace the values of '%%%' by '' also for data of type string.
    post_configure : bool, optional
        An option to indicate whether the input XMLs has been processed by
        configure or not.
    only_allocated : bool, optional
        It will include only information from allocated fibres.
    sort_targets : bool, optional
        Sort the data inside each field grouping by IFU bundle.

    Returns
    -------
    data_dict : dict
        A dictionary with the data. Its keys are the name of the SPA columns
        which are expected to be potencially in the configure XML files, while
        its values are lists containing these data.
    """
    data_dict = \
        _get_spa_data_from_targuse(
            xml_filename_list, targuse_list=['T', 'S', 'R'],
            replace_triple_percent=replace_triple_percent,
            post_configure=post_configure, only_allocated=only_allocated,
            sort_targets=sort_targets)

    return data_dict


def get_spa_data_of_target_calib_random_and_sky_fibres_from_xmls(
        xml_filename_list, replace_triple_percent=True, post_configure=True,
        only_allocated=True, sort_targets=False):
    """
    Get SPA data of the T/C/R/S fibres contained in a list of XML files.

    Parameters
    ----------
    xml_filename_list : list of str
        A list of configure XML files.
    replace_triple_percent : bool, optional
        It will replace the values of '%%%' by '' also for data of type string.
    post_configure : bool, optional
        An option to indicate whether the input XMLs has been processed by
        configure or not.
    only_allocated : bool, optional
        It will include only information from allocated fibres.
    sort_targets : bool, optional
        Sort the data inside each field grouping by IFU bundle.

    Returns
    -------
    data_dict : dict
        A dictionary with the data. Its keys are the name of the SPA columns
        which are expected to be potencially in the configure XML files, while
        its values are lists containing these data.
    """
    data_dict = \
        _get_spa_data_from_targuse(
            xml_filename_list, targuse_list=['C', 'T', 'S', 'R'],
            replace_triple_percent=replace_triple_percent,
            post_configure=post_configure, only_allocated=only_allocated,
            sort_targets=sort_targets)

    return data_dict


def get_data_per_field_of_one_xml(xml_filename):
    """
    Get SPA data of a XML file divided per field.

    Parameters
    ----------
    xml_filename : str
        A configure XML files.

    Returns
    -------
    data_dict_list : list of dict
        A list of dictionaries with the data. Each element of the list
        corresponds to a different field. Each key of each dictionary is the
        name of a target attribute which are expected to be potencially in the
        configure XML file, while its values are lists containing these data.
    """
    
    data_dict_list = _get_data_from_targuse_per_field_of_one_xml(
        xml_filename, targuse_list=['C', 'T', 'S', 'R'],
        replace_triple_percent=True, post_configure=True, only_allocated=True,
        sort_targets=False, namespace='xml')
    
    return data_dict_list


def _get_attribute_in_simple_element_of_xml_file(xml_file, element_name,
                                                 attribute_name):
    
    dom = _minidom.parse(xml_file)
    
    element = dom.getElementsByTagName(element_name)[0]
    
    attribute = element.getAttribute(attribute_name)
    
    return attribute


def _get_single_value_from_list(value_list):

    value_set = set(value_list)
    
    assert len(value_set) == 1
    
    value = value_list[0]
    
    return value


def get_halims_from_xmls(input_xmls):
    """
    Get the ha limit ranges from a list of XML files or a single file.

    Parameters
    ----------
    input_xmls : list of str or str
        A list of XML filenames or a single filename.

    Returns
    -------
    ha_range : tuple
        The min/max ha_range defined in the XML
    """
    
    if type(input_xmls) is list:
        xml_filename_list = input_xmls
    else:
        xml_filename_list = [input_xmls]
    
    
    min_list = [_get_attribute_in_simple_element_of_xml_file(
                          xml_file, 'hour_angle_limits', 'earliest')
                      for xml_file in xml_filename_list]

    max_list = [_get_attribute_in_simple_element_of_xml_file(
                          xml_file, 'hour_angle_limits', 'latest')
                      for xml_file in xml_filename_list]

    ha_range = [(ha_min,ha_max) for ha_min, ha_max in zip(min_list,max_list)]
    
    # trimester = _get_single_value_from_list(trimester_list)
    
    return ha_range



def get_trimester_from_xmls(input_xmls):
    """
    Get the trimester from a list of XML files or a single file.

    Parameters
    ----------
    input_xmls : list of str or str
        A list of XML filenames or a single filename.

    Returns
    -------
    trimester : str
        The trimester present in the XML file.
    """
    
    if type(input_xmls) is list:
        xml_filename_list = input_xmls
    else:
        xml_filename_list = [input_xmls]
    
    
    trimester_list = [_get_attribute_in_simple_element_of_xml_file(
                          xml_file, 'observation', 'trimester')
                      for xml_file in xml_filename_list]
    
    trimester = _get_single_value_from_list(trimester_list)
    
    return trimester


def get_tacalloc_from_xmls(input_xmls):
    """
    Get the TACALLOC code from a list of XML files or a single file.

    Parameters
    ----------
    input_xmls : list of str or str
        A list of XML filenames or a single filename.

    Returns
    -------
    tacalloc : str
        The tac_alloc attribute present in the XML file.
    """
    
    if type(input_xmls) is list:
        xml_filename_list = input_xmls
    else:
        xml_filename_list = [input_xmls]
    
    
    tacalloc_list = [_get_attribute_in_simple_element_of_xml_file(
        xml_file, 'observation', 'tac_alloc')
                    for xml_file in xml_filename_list]
    
    tacalloc = _get_single_value_from_list(tacalloc_list)
    
    return tacalloc



def get_obsmode_from_xmls(input_xmls):
    """
    Get the obsmode from a list of XML files or a single file.

    Parameters
    ----------
    input_xmls : list of str or str
        A list of XML filenames or a single filename.

    Returns
    -------
    obsmode : str
        The obsmode present in the XML file.
    """
    
    if type(input_xmls) is list:
        xml_filename_list = input_xmls
    else:
        xml_filename_list = [input_xmls]
    
    
    obsmode_list = [_get_attribute_in_simple_element_of_xml_file(
                          xml_file, 'observation', 'obs_mode')
                      for xml_file in xml_filename_list]
    
    obsmode = _get_single_value_from_list(obsmode_list)
    
    return obsmode


def get_author_from_xmls(input_xmls):
    """
    Get the author from a list of XML files or a single file.

    Parameters
    ----------
    input_xmls : list of str or str
        A list of XML filenames or a single filename.

    Returns
    -------
    author : str
        The author present in the XML file.
    """
    
    if type(input_xmls) is list:
        xml_filename_list = input_xmls
    else:
        xml_filename_list = [input_xmls]
    
    author_list = [_get_attribute_in_simple_element_of_xml_file(
                       xml_file, 'weave', 'author')
                   for xml_file in xml_filename_list]
    
    author = _get_single_value_from_list(author_list)
    
    return author


def get_cc_report_from_xmls(input_xmls):
    """
    Get the cc_report from a list of XML files or a single file.

    Parameters
    ----------
    input_xmls : list of str or str
        A list of XML filenames or a single filename.

    Returns
    -------
    cc_report : str
        The cc_report present in the XML file.
    """
    
    if type(input_xmls) is list:
        xml_filename_list = input_xmls
    else:
        xml_filename_list = [input_xmls]
    
    cc_report_list = [_get_attribute_in_simple_element_of_xml_file(xml_file,
                          'weave', 'cc_report')
                      for xml_file in xml_filename_list]
    
    cc_report = _get_single_value_from_list(cc_report_list)
    
    return cc_report


def get_report_verbosity_from_xmls(input_xmls):
    """
    Get the report_verbosity from a list of XML files or a single file.

    Parameters
    ----------
    input_xmls : list of str or str
        A list of XML filenames or a single filename.

    Returns
    -------
    report_verbosity : str
        The report_verbosity present in the XML file.
    """
    
    if type(input_xmls) is list:
        xml_filename_list = input_xmls
    else:
        xml_filename_list = [input_xmls]
    
    report_verbosity_list = [_get_attribute_in_simple_element_of_xml_file(
                                 xml_file, 'weave', 'report_verbosity')
                      for xml_file in xml_filename_list]
    
    report_verbosity = _get_single_value_from_list(report_verbosity_list)
    
    report_verbosity = int(report_verbosity)
    
    return report_verbosity


def get_datamver_from_xmls(input_xmls):
    """
    Get the datamver from a list of XML files or a single file.

    Parameters
    ----------
    input_xmls : list of str or str
        A list of XML filenames or a single filename.

    Returns
    -------
    datamver : str
        The datamver present in the XML file.
    """
    
    if type(input_xmls) is list:
        xml_filename_list = input_xmls
    else:
        xml_filename_list = [input_xmls]
    
    datamver_list = [_get_attribute_in_simple_element_of_xml_file(xml_file,
                          'weave', 'datamver')
                      for xml_file in xml_filename_list]
    
    datamver = _get_single_value_from_list(datamver_list)
    
    return datamver


def get_obs_mode(input_xml):
    """
    Get the obs_mode from a XML file.

    Parameters
    ----------
    input_xml : str
        A XML filename.

    Returns
    -------
    obs_mode : str
        The obs_mode present in the XML file.
    """
    
    obs_mode = _get_attribute_in_simple_element_of_xml_file(
                   input_xml, 'observation', 'obs_mode')
    
    return obs_mode


def get_coord_of_first_field(input_xml):
    """
    Get the coordinates of the first field from a XML file.

    Parameters
    ----------
    input_xml : str
        A XML filename.

    Returns
    -------
    ra : float
        The right ascension of the first field.
    dec : float
        The declination of the first field.
    """
    
    ra = float(_get_attribute_in_simple_element_of_xml_file(
                   input_xml, 'field', 'RA_d'))
    dec = float(_get_attribute_in_simple_element_of_xml_file(
                    input_xml, 'field', 'Dec_d'))
    
    return ra, dec

