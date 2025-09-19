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
import glob
import logging
import os
import re
import xml.dom.minidom
from collections import OrderedDict
from astropy.time import Time
import numpy as np
from astropy.io import fits

from workflow.utils import get_version
from workflow.utils.classes import OBXML
from workflow.utils.get_obstemp_info import (get_obstemp_info,
                                             get_obstemp_dict)
from workflow.utils.get_progtemp_info import (get_progtemp_dict,
                                              get_progtemp_info,
                                              get_obsmode_from_progtemp)
from workflow.utils.get_resources import (get_blank_xml_template,
                                          get_progtemp_file,
                                          get_obstemp_file)


class _IFUDriverCat:
    """
    Convert the target-level data from an IFU driver catalogue to a set of XMLs.
    
    This class provides code to take an IFU driver catalogue containing
    information of the central fibres and convert it into XML targets that are
    written into a supplied XML.
    
    Parameters
    ----------
    filename : str
        A FITS file containing an IFU driver cat.
    targcat : str, optional
        Filename of the FITS catalogue which will be submitted to WASP (if None,
        it will be derived from the input filename).
    """

    
    def __init__(self, filename, targcat=None):

        # Save the filename

        self.filename = filename

        # Set the default prefix to be used for output XMLs

        self.default_prefix = os.path.splitext(
            os.path.basename(self.filename))[0].replace('-ifu_driver_cat', '')

        if self.default_prefix[-1] != '-':
            self.default_prefix = self.default_prefix + '-'

        # Set the targcat value

        if targcat is None:
            self.targcat = os.path.basename(self.filename).replace(
                '-ifu_driver_cat', '')
        else:
            self.targcat = targcat

        # Create a blank OB name dictionary (this is written to the
        # <observation name> attribute)
        self.obs_names = {}
            
        # Read and save the information into the IFU driver catalogue
        
        logging.info('Reading {}'.format(self.filename))

        with fits.open(self.filename) as hdu_list:

            self.datamver = hdu_list[0].header['DATAMVER']
            self.trimester = hdu_list[0].header['TRIMESTE']
            self.report_verbosity = hdu_list[0].header['VERBOSE']
            self.tacalloc = hdu_list[0].header['TACALLOC']
            self.author = hdu_list[0].header['CAT_MAIL']
            self.cc_report = hdu_list[0].header['CAT_CC']

            self.data = hdu_list[1].data

    def _get_trimester_midpoint(self,trimester):
        mpoint_month = {'A1':3,'A2':6,'B1':9,'B2':12}
        mp = mpoint_month[trimester[-2:]]
        yr = int(trimester[:4])
        tstr = '{}-{}-15T00:00:00'.format(yr,str(mp).zfill(2))
        tmid = Time(tstr, format='isot', scale='utc')
        return tmid
        
    def _get_previous_files(self, output_dir='', prefix=None, suffix='-t'):

        if prefix is None:
            prefix = self.default_prefix

        glob_pattern = os.path.join(output_dir, '{}*{}.xml'.format(prefix,
                                                                   suffix))
        regex = '^{}[lm]ifu_[0-9]+{}.xml$'.format(prefix, suffix)

        candidate_list = glob.glob(glob_pattern)
        candidate_list.sort()

        filename_list = [filename for filename in candidate_list
                        if re.match(regex, os.path.basename(filename))]

        return filename_list

    
    def remove_xmls(self, output_dir='', prefix=None, suffix='-t'):

        prev_filename_list = self._get_previous_files(output_dir=output_dir,
                                                      prefix=prefix,
                                                      suffix=suffix)

        if len(prev_filename_list) > 0:

            logging.info('Removing previous files: {}'.format(
                prev_filename_list))

            for filename in prev_filename_list:
                os.remove(filename)


    def _get_obs_name(self,data,delimiter='_',illegal=['-', '\\', '/', ':', '*', '?', '<', '>', '|', ' ']):

        assert delimiter not in illegal, "The chosen delimiter ({}) is an OCS illegal character"

        # is this lifu or mifu?
        is_lifu = False
        ins_mode = 'm'
        if str(data['PROGTEMP'][0]) in ['4','5','6']:
            is_lifu = True
            # do LIFU rules (include the targid)
            ins_mode = 'L'
            start_len = 5
            # mode/targsrvy/targname/targid/dither/index
        else:            
            ins_mode = 'm'
            start_len = 4
            # mode/targsrvy/targname/dither/index
            
        targsrvy = data['TARGSRVY']
        if (len(targsrvy) == len('WS2023A1-022')) and (targsrvy[-4] == '-') and (targsrvy[:2] in ['WV','WS']):
            # this is an SV/opentime survey. Just want the numeric code
            targsrvy = data['TARGSRVY'].split('-')[1]
            
        _name = [ins_mode,targsrvy,data['TARGNAME']]
        if is_lifu:
            _name.append(data['TARGID'])

        dither = "d{}".format(str(data['IFU_DITHER']))
        if data['IFU_DITHER'] < -1:
            dither = "dm{}".format(str(abs(data['IFU_DITHER'])))
        elif data['IFU_DITHER'] == -1:
            dither = "dc"
        _name.append(dither)
            
        obs_name = delimiter.join(_name)
        _on = obs_name
        while _on:
            try:
                _on = self.obs_names[obs_name]
                if len(_name) == start_len:
                    # need to add the index, set to 1 so it will +1 to 2
                    _name.append('1')
                
                _name[-1] = str(int(_name[-1])+1)
                obs_name = delimiter.join(_name)
            except KeyError:
                _on = None
                self.obs_names[obs_name] = 1


        # strip out any illegal characters - they could have come from targsrvy, targid or targname
        obs_name = ''.join(delimiter if c in illegal else c for c in obs_name)
        return obs_name
                
    def _get_output_path(self, obsmode, output_dir='', prefix=None,
                         suffix='-t'):

        # Set the prefix of the output file if it has not been provided

        if prefix is None:
            prefix = self.default_prefix

        # Choose the first filename which does not exist

        index = 0
        output_path = ''
        
        while (index < 1) or os.path.exists(output_path):

            index += 1

            output_basename = '{}{}_{:02d}{}.xml'.format(
                prefix, obsmode.lower(), index, suffix)

            output_path = os.path.join(output_dir, output_basename)

        return output_path

    
    def _process_ob(self, entry_group, xml_template, progtemp_dict,
                    obstemp_dict, spatial_binning, output_dir='', prefix=None,
                    suffix='-t'):

        # Get some information from the first entry of the group

        first_entry = entry_group[0]

        progtemp = first_entry['PROGTEMP']
        obstemp = first_entry['OBSTEMP']
        ifu_dither = first_entry['IFU_DITHER']
        ifu_pa_request = first_entry['IFU_PA_REQUEST']
        ha_range = first_entry['HA_RANGE']        
        ha_req = first_entry['HA_REQUEST']
        tac_id = first_entry['TAC_ID'] 
        
        # Get some information from all the entries in the group

        targsrvy_list = list(set([entry['TARGSRVY'] for entry in entry_group]))
        targsrvy_list.sort()

        # Guess some parameters which depends on PROGTEMP
        
        chained = ('+' in progtemp)

        spectrograph_dict = get_progtemp_info(progtemp,
                                              progtemp_dict=progtemp_dict)

        obsmode = spectrograph_dict['obsmode']

        assert obsmode in ['LIFU', 'mIFU']

        assert (spectrograph_dict['red_resolution'] ==
                spectrograph_dict['blue_resolution'])
        resolution = spectrograph_dict['red_resolution']
                
        red_vph = spectrograph_dict['red_vph']
        blue_vph = spectrograph_dict['blue_vph']

        assert (spectrograph_dict['red_num_exposures'] ==
                spectrograph_dict['blue_num_exposures'])
        num_science_exposures = spectrograph_dict['red_num_exposures']
        
        assert (spectrograph_dict['red_exp_time'] ==
                spectrograph_dict['blue_exp_time'])
        science_exptime = spectrograph_dict['red_exp_time']

        assert (spectrograph_dict['red_binning_x'] ==
                spectrograph_dict['blue_binning_x'])
        binning_x = spectrograph_dict['red_binning_x']

        # Set the spatial binning from the input of this method
        
        binning_y = spatial_binning

        # Set the some paremeters which depends on OBSMODE:
        #  - Name of the observation
        #  - Plate
        #  - max_fibres

        if obsmode == 'LIFU':
            observation_name = '{}-{}'.format(first_entry['TARGNAME'],
                                              first_entry['TARGID'])
            observation_name = self._get_obs_name(first_entry,delimiter='_')
            plate = 'LIFU'
            max_calibration = 603
            max_guide = 1
            max_sky = 603
            num_sky_fibres = 603
            max_fibres = 603
        elif obsmode == 'mIFU':
            # observation_name = first_entry['TARGNAME']
            observation_name = self._get_obs_name(first_entry,delimiter='_')

            plate = 'PLATE_B'
            max_calibration = 37*20
            max_guide = 8
            max_sky = 37*20
            num_sky_fibres = 740
            max_fibres = 740

        # Guess the number of dither positions

        if (obsmode == 'LIFU') and (len(entry_group) > 1):
            num_dither_positions = len(entry_group)
        else:
            num_dither_positions = ifu_dither

        if abs(ifu_dither) >= 3:
            # this is not a custom dither, not an undithered obs
            assert abs(ifu_dither) == num_science_exposures
        elif ifu_dither < 0:
            # custom dither - the number of entries == number of sci exp?
            assert len(entry_group) == num_science_exposures
            
        # Guess some parameters from obstemp

        obsconstraints_dict = get_obstemp_info(obstemp,
                                               obstemp_dict=obstemp_dict)

        elevation_min = obsconstraints_dict['elevation_min']
        moondist_min = obsconstraints_dict['moondist_min']
        seeing_max = obsconstraints_dict['seeing_max']
        skybright_max = obsconstraints_dict['skybright_max']
        transparency_min = obsconstraints_dict['transparency_min']

        # Set the position angle (if possible)

        if obsmode == 'LIFU':

            if not np.isnan(ifu_pa_request):
                pa = ifu_pa_request
            else:
                pa = None

        elif obsmode == 'mIFU':

            pa = None

        # Create an OB from the XML template

        ob_xml = OBXML(xml_template)

        # Remove the elements and attributes which are configure outputs
                        
        ob_xml.remove_configure_outputs()

        # Remove the non-used elements
                        
        ob_xml.remove_non_used_elements()

        # Set the attributes of the root element
        
        comment = 'WIW{}'.format(get_version())
        
        root_attrib_dict = {
            'author': self.author,
            'cc_report': self.cc_report,
            'comment': comment,
            'report_verbosity': self.report_verbosity
        }

        ob_xml.set_root_attrib(root_attrib_dict)

        # Set the contents of the spectrograph element

        ob_xml.set_spectrograph(binning_x=binning_x, binning_y=binning_y,
                                resolution=resolution, red_vph=red_vph,
                                blue_vph=blue_vph)

        # Set the contents of the exposures element

        ob_xml.set_exposures(num_science_exposures, science_exptime)

        # Set the attributes of the observation element

        observation_attrib_dict = {
            'chained': chained,
            'name': observation_name,
            'obstemp': obstemp,
            'obs_mode': obsmode,
            'progtemp': progtemp,
            'trimester': self.trimester,
            'tac_alloc': self.tacalloc,
            'tac_id': tac_id
        }

        if pa is not None:
            observation_attrib_dict['pa'] = pa

        ob_xml.set_observation(observation_attrib_dict)
            
        # Set the attributes of the dithering element

        configure_attrib_dict = {
            'max_calibration': max_calibration,
            'max_guide': max_guide,
            'max_sky': max_sky,
            'num_sky_fibres': num_sky_fibres,
            'plate': plate,
            'seed': ''
        }

        ob_xml.set_configure(configure_attrib_dict)

        # Set the attributes of the obsconstraints element

        obsconstraints_attrib_dict = {
            'elevation_min': elevation_min,
            'moondist_min': moondist_min,
            'seeing_max': seeing_max,
            'skybright_max': skybright_max,
            'transparency_min': transparency_min
        }

        ob_xml.set_obsconstraints(obsconstraints_attrib_dict)

        # Set the attributes of the dithering element

        dithering_attrib_dict = {
            'apply_dither': ifu_dither
        }

        ob_xml.set_dithering(dithering_attrib_dict)

        # Set the contents of the survey element

        ob_xml.set_surveys(targsrvy_list, max_fibres)

        # Set the contents of the fields element

        ob_xml.set_fields(obsmode, entry_group, targcat=self.targcat,
                          num_science_exposures=num_science_exposures)

        # Write the OB XML to a file

        output_path = self._get_output_path(obsmode, output_dir=output_dir,
                                            prefix=prefix, suffix=suffix)

        ob_xml.write_xml(output_path)

        # create the epoch file
        with open(output_path.replace('-t.xml','.epoch'),'w') as file:
            tmid = self._get_trimester_midpoint(self.trimester)
            file.write("{:.9f}".format(tmid.decimalyear))

        # create the ha_limits file
        with open(output_path.replace('-t.xml','.ha'),'w') as file:
            file.write("{:.2f}".format(ha_req))

        
        return output_path
        

        
    def _generate_lifu_xmls(self, lifu_entry_list, xml_template, progtemp_dict,
                            obstemp_dict, spatial_binning, output_dir='',
                            prefix='', suffix='-t'):
        
        output_file_list = []

        # How do you group entries belonging to the same OB?
        # (for custom dithers)

        group_id = ('TARGID', 'TARGNAME', 'PROGTEMP', 'OBSTEMP', 'TAC_ID')

        # Detect the non-custom dithers and save them into a list

        non_custom_dithers_list = []
        
        for lifu_entry in lifu_entry_list:

            if lifu_entry['IFU_DITHER'] != -1:

                non_custom_dithers_list.append(lifu_entry)

        # Detect the custom dithers and save them into a dict

        custom_dithers_dict = OrderedDict()

        for lifu_entry in lifu_entry_list:

            if lifu_entry['IFU_DITHER'] == -1:

                key = tuple(lifu_entry[col] for col in group_id)

                if key not in custom_dithers_dict.keys():
                    custom_dithers_dict[key] = []
                
                custom_dithers_dict[key].append(lifu_entry)

        # Generate the IFU XMLs of the non-custom dithers

        logging.info(
            'Processing {} non-custom dither OBs for LIFU'.format(
                len(non_custom_dithers_list)))


        for lifu_entry in non_custom_dithers_list:

            entry_group = [lifu_entry]

            output_file = self._process_ob(entry_group, xml_template,
                                           progtemp_dict, obstemp_dict,
                                           spatial_binning,
                                           output_dir=output_dir, prefix=prefix,
                                           suffix=suffix)
            
            output_file_list.append(output_file)


                
        # Process the custom dither pointing groupings
        
        logging.info(
            'Processing {} custom dither OBs for LIFU'.format(
                len(custom_dithers_dict)))

        for entry_group in custom_dithers_dict.values():
        
            output_file = self._process_ob(entry_group, xml_template,
                                           progtemp_dict, obstemp_dict,
                                           spatial_binning,
                                           output_dir=output_dir, prefix=prefix,
                                           suffix=suffix)
            
            output_file_list.append(output_file)

        return output_file_list

                
    def _generate_mifu_xmls(self, mifu_entry_list, xml_template, progtemp_dict,
                            obstemp_dict, spatial_binning, mifu_mode='all',
                            mifu_num_calibs=2, mifu_num_extra=0, output_dir='',
                            prefix='', suffix='-t'):

        # What happens if there are (e.g.) 100 bundles in a given key and
        # mifu_num_calibs is 2?
        #
        # User needs to decide 3 options:
        # - all: Include all
        # - aufbau: Aufbau principle, i.e. fill up to 18, leaving 2 for
        #           calibration bundles) then make a new XML
        # - equipartition: divide up the bundles equally:
        #                  100 / 18 = 5.5 --> 6 OBs
        #                  100 / 6 = 17 bundles per OB
        #                            (with remainder added to a final OB)
        #                  i.e. 5 OBs x 17 bundles + 1 OB x 15 bundles

        assert mifu_mode in ['all', 'aufbau', 'equipartition']
        assert (mifu_num_calibs >= 0) and (mifu_num_calibs < 20)
        
        output_file_list = []

        # How do you group bundles belonging to the same field?

        group_id = ('TARGNAME', 'PROGTEMP', 'OBSTEMP', 'IFU_DITHER')

        # Group the mIFU entries per field

        fields_dict = OrderedDict()

        for mifu_entry in mifu_entry_list:

            key = tuple(mifu_entry[col] for col in group_id)

            if key not in fields_dict.keys():
                fields_dict[key] = []

            fields_dict[key].append(mifu_entry)

        # Group the targets in each field per OB

        ob_nested_list = []

        for key in fields_dict.keys():

            # Get the entries in the field

            field_entry_list = fields_dict[key]

            # Guess the number of entries in the field

            num_entries_in_field = len(field_entry_list)

            # Create the group of mIFU entries per OB depending on the mode

            if (mifu_mode == 'aufbau') or (mifu_mode == 'equipartition'):

                # Maximum number of target entries per OB

                max_entries_per_ob = 20 - mifu_num_calibs

                if mifu_mode == 'aufbau':

                    # All the OBs (except the last one) will contain the maximum
                    # number of target entries plus the amount of extra targets
                    
                    num_entries_per_ob = max_entries_per_ob + mifu_num_extra

                elif mifu_mode == 'equipartition':

                    # Compute the number of OBs

                    num_obs = int(np.ceil(num_entries_in_field /
                                          max_entries_per_ob))

                    # Divide the targets to be the same (except the last one,
                    # due to rounding effects)

                    num_entries_per_ob = int(np.ceil(num_entries_in_field /
                                                     num_obs))

                else:

                    raise ValueError

                # Save the targets of each OB using Python slices
                # (do not worry, if stop is higher than the lenght of the list,
                #  this is allowed and the trick)
                
                start = 0

                while start < num_entries_in_field:

                    stop = start + num_entries_per_ob

                    ob_nested_list.append(field_entry_list[start:stop])

                    start = stop
                    
            elif mifu_mode == 'all':

                ob_nested_list.append(field_entry_list)

            else:

                raise ValueError

        # Proccess OB grouping

        logging.info(
            'Processing {} mIFU fields into {} OBs according to {} mode'.format(
                len(fields_dict), len(ob_nested_list), mifu_mode))

        for entry_group in ob_nested_list:
        
            output_file = self._process_ob(entry_group, xml_template,
                                           progtemp_dict, obstemp_dict,
                                           spatial_binning,
                                           output_dir=output_dir, prefix=prefix,
                                           suffix=suffix)
            
            output_file_list.append(output_file)

        return output_file_list

                    
    def generate_xmls(self, xml_template, progtemp_file=None,
                      obstemp_file=None, mifu_mode='all', mifu_num_calibs=2,
                      mifu_num_extra=0, spatial_binning=1, output_dir='',
                      prefix='', suffix='-t', pass_datamver=False):

        # Get the DATAMVER of the XML template

        xml_template_dom = xml.dom.minidom.parse(xml_template)
        xml_datamver = xml_template_dom.childNodes[0].getAttribute('datamver')

        # Get the dictionaries to interpret PROGTEMP and OBSTEMP, and their
        # DATAMVER values

        progtemp_datamver, progtemp_dict, forbidden_dict = get_progtemp_dict(
            filename=progtemp_file, assert_orb=True)

        obstemp_datamver, obstemp_dict = get_obstemp_dict(
            filename=obstemp_file)

        # Check DATAMVER of the IFU driver cat, the XML template, PROGTEMP file
        # and OBSTEMP file are consistent

        if self.datamver != xml_datamver:
            logging.critical(
               'DATAMVER mismatch ({} != {}) for XML template: '.format(
                    self.datamver, xml_datamver) +
                'Stop unless you are sure!')

            if pass_datamver == False:
                raise SystemExit(2)

        if self.datamver != progtemp_datamver:
            logging.critical(
               'DATAMVER mismatch ({} != {}) for PROGTEMP file: '.format(
                    self.datamver, progtemp_datamver) +
                'Stop unless you are sure!')

            if pass_datamver == False:
                raise SystemExit(2)

        if self.datamver != obstemp_datamver:
            logging.critical(
               'DATAMVER mismatch ({} != {}) for OBSTEMP file: '.format(
                    self.datamver, obstemp_datamver) +
                'Stop unless you are sure!')

            if pass_datamver == False:
                raise SystemExit(2)

        # Classify the entries into LIFU and mIFU
        
        lifu_entry_list = []
        mifu_entry_list = []

        for i, entry in enumerate(self.data):

            progtemp = entry['PROGTEMP']

            try:
                obsmode = get_obsmode_from_progtemp(progtemp,
                                                    progtemp_dict=progtemp_dict)
            except:
                obsmode = None

            if obsmode == 'LIFU':
                lifu_entry_list.append(entry)
            elif obsmode == 'mIFU':
                mifu_entry_list.append(entry)
            else:
                logging.warning('unexpected PROGTEMP in row {}: {}'.format(
                    i + 1, progtemp))

        # Generate the LIFU XMLs
                
        lifu_output_file_list = self._generate_lifu_xmls(
            lifu_entry_list, xml_template, progtemp_dict, obstemp_dict,
            spatial_binning, output_dir=output_dir, prefix=prefix,
            suffix=suffix)

        # Generate the mIFU XMLs
                
        mifu_output_file_list = self._generate_mifu_xmls(
            mifu_entry_list, xml_template, progtemp_dict, obstemp_dict,
            spatial_binning, mifu_mode=mifu_mode,
            mifu_num_calibs=mifu_num_calibs, mifu_num_extra=mifu_num_extra,
            output_dir=output_dir, prefix=prefix, suffix=suffix)
        
        # Get the list with the output files
        
        output_file_list = lifu_output_file_list + mifu_output_file_list
        
        return output_file_list


def create_xml_files(ifu_driver_cat_filename, output_dir, xml_template,
                     progtemp_file=None, obstemp_file=None,
                     mifu_mode='all', mifu_num_calibs=2, mifu_num_extra=0,
                     spatial_binning=1, targcat=None, prefix=None, suffix='-t',
                     pass_datamver=False, overwrite=False):
    """
    Create XML files with targets from an IFU driver cat.
    
    Parameters
    ----------
    ifu_driver_cat_filename : str
        A FITS file containing an IFU driver cat.
    output_dir : str
        Name of the directory which will containe the output XML files.
    xml_template : str
        A blank XML template to be populated with the information of the OBs.
    progtemp_file : str, optional
        A progtemp.dat file with the definition of PROGTEMP.
    obstemp_file : str, optional
        A obstemp.dat file with the definition of OBSTEMP.
    mifu_mode : {'all', 'aufbau', 'equipartition'}, optional
        Grouping mode for mIFU targets.
    mifu_num_calibs : int, optional
        Number of mIFU calibration stars and sky bundles to be considered in
        the mIFU grouping.
    mifu_num_extra : int, optional
        Number of extra mIFU targets to be included in the OBs when 'aufbau'
        grouping mode is used.
    spatial_binning : int, optional
        Number of pixels for the spatial binning.
    targcat : str, optional
        Filename of the FITS catalogue which will be submitted to WASP (if None,
        it will be derived from the input filename).
    prefix : str, optional
        Prefix to be used in the output files (it will be derived from
        ifu_driver_cat_filename if None is provided).
    suffix : str, optional
        Suffix to be used in the output files.
    pass_datamver : bool, optional
        Continue even if DATAMVER mismatch is detected.
    overwrite : bool, optional
        Overwrite the output FITS file.

    Returns
    -------
    output_file_list : list of str
        A list with the output XML files.
    """

    # Check that the input IFU driver cat exists and is a file

    assert os.path.isfile(ifu_driver_cat_filename)

    # Create an object with the IFU driver cat

    ifu_driver_cat = _IFUDriverCat(ifu_driver_cat_filename)

    # Remove the previous files if overwriting has been requested
    
    if overwrite == True:
        ifu_driver_cat.remove_xmls(output_dir=output_dir, prefix=prefix,
                                   suffix=suffix)

    # Create the XML files

    output_file_list = ifu_driver_cat.generate_xmls(
        xml_template, progtemp_file=progtemp_file, obstemp_file=obstemp_file,
        mifu_mode=mifu_mode, mifu_num_calibs=mifu_num_calibs,
        mifu_num_extra=mifu_num_extra, spatial_binning=spatial_binning,
        output_dir=output_dir, prefix=prefix, suffix=suffix,
        pass_datamver=pass_datamver)
    
    return output_file_list


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
             description='Create XML files with targets from an IFU driver cat')

    parser.add_argument('ifu_driver_cat',
                        help="""a FITS file containing an IFU driver cat""")

    parser.add_argument('--xml_template',
                        default=os.path.join('aux', 'BlankXMLTemplate.xml'),
                        help="""a blank XML template to be populated with the
                        information of the OBs""")

    parser.add_argument('--progtemp_file',
                        default=os.path.join('aux', 'progtemp.dat'),
                        help="""a progtemp.dat file with the definition of
                        PROGTEMP""")


    parser.add_argument('--obstemp_file',
                        default=os.path.join('aux', 'obstemp.dat'),
                        help="""a obstemp.dat file with the definition of
                        OBSTEMP""")

    parser.add_argument('--mifu_mode', default='all',
                        choices=['all', 'aufbau', 'equipartition', 'all'],
                        help="""grouping mode for mIFU targets: 'all' includes
                        all, 'aufbau' follows the aufbau principle,
                        'equipartition' divides up the bundles equally""")

    parser.add_argument('--mifu_num_calibs', default=2,
                        choices=range(20), type=int,
                        help="""number of mIFU calibration stars and sky bundles
                        to be considered in the mIFU grouping""")

    parser.add_argument('--mifu_num_extra', default=0, type=int,
                        help="""number of extra mIFU targets to be included in
                        the OBs when 'aufbau' grouping mode is used""")

    parser.add_argument('--spatial_binning', default=1, choices=range(1, 5),
                        type=int,
                        help='number of pixels for the spatial binning')

    parser.add_argument('--targcat', default=None, type=str,
                        help="""filename of the FITS catalogue which will be
                        submitted to WASP (if None, it will be derived from the
                        input filename)""")

    parser.add_argument('--outdir', dest='output_dir', default='output',
                        help="""name of the directory which will contain the
                        output XML files""")

    parser.add_argument('--prefix', dest='prefix', default=None,
                        help="""prefix to be used in the output files (it will
                        be derived from ifu_driver_cat if non provided)""")

    parser.add_argument('--pass_datamver', dest='pass_datamver',
                        action='store_true',
                        help='continue even if DATAMVER mismatch is detected')

    parser.add_argument('--overwrite', dest='overwrite', action='store_true',
                        help='overwrite the output files')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    if not os.path.exists(args.output_dir):
        logging.info('Creating the output directory')
        os.mkdir(args.output_dir)

    xml_template_dir = os.path.dirname(args.xml_template)
    if not os.path.exists(xml_template_dir):
        logging.info('Creating the directory of the blank XML template')
        os.mkdir(xml_template_dir)

    if not os.path.exists(args.xml_template):
        logging.info('Downloading the blank XML template')
        get_blank_xml_template(file_path=args.xml_template)

    if not os.path.exists(args.progtemp_file):
        logging.info('Downloading the progtemp file')
        get_progtemp_file(file_path=args.progtemp_file)

    if not os.path.exists(args.obstemp_file):
        logging.info('Downloading the obstemp file')
        get_obstemp_file(file_path=args.obstemp_file)

    create_xml_files(args.ifu_driver_cat, args.output_dir, args.xml_template,
                     progtemp_file=args.progtemp_file,
                     obstemp_file=args.obstemp_file,
                     mifu_mode=args.mifu_mode,
                     mifu_num_calibs=args.mifu_num_calibs,
                     mifu_num_extra=args.mifu_num_extra,
                     spatial_binning=args.spatial_binning, targcat=args.targcat,
                     prefix=args.prefix, pass_datamver=args.pass_datamver,
                     overwrite=args.overwrite)

