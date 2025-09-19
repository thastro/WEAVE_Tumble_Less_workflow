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
import os as _os
import re as _re
import tempfile as _tempfile

from workflow.utils.get_resources import get_progtemp_file as _get_progtemp_file


def get_progtemp_dict(filename=None, assert_orb=True, cache=False):
    """
    Get a dictionary to interpret progtemp values.

    Parameters
    ----------
    filename : str
        The name of the file formatted like progtemp.dat.
    assert_orb : bool, optional
        Assert whether the reduntand information in the file is not
        contradictory.
    cache : bool, optional
        If true, will use locally-cached version of progtemp.dat file, rather than
        retrieve from the internet
    
    Returns
    -------
    datamver : str
        The value of DATAMVER in the file.
    progtemp_dict : dict
        A dictionary with the information to interpret a progtemp value.
    forbidden_dict : dict
        A dictionary with the information of the forbidden progtem values


    """
    
    if (filename is None) and (cache):
        # set the filename, ensure it exists
        filename = _os.path.join(_os.path.dirname(_os.path.dirname(__file__)),'aux/progtemp.dat')
        if not _os.path.isfile(filename):
            logging.warning("User set cache=True but file {} does not exist - will download".format(filename))
            
            _get_progtemp_file(file_path=filename)            
        else:
            logging.warning("Using locally-cached copy of PROGTEMP file at {} - we advise disabling cache mode to ensure this file is up to date".format(
            filename))
            
    elif filename is None:
        filename = _os.path.join(_tempfile.mkdtemp(), 'progtemp.dat')
        _get_progtemp_file(file_path=filename)

    # Set of the regex expected to be found in the file

    datamver_regex = '^#DATAMVER (.+)$'
    header_regex = '^[0-9_]{5}:$'
    n_regex_option_a = (
        '^([0-9])____: '
        'observation:obs_mode=([A-Za-z]+) '
        'programme:spectrograph:red_Arm:resolution=([a-z]+) '
        'programme:spectrograph:red_Arm:VPH=([A-Z0-9]+) '
        'programme:spectrograph:blue_Arm:resolution=([a-z]+) '
        'programme:spectrograph:blue_Arm:VPH=([A-Z0-9]+)$'
    )
    n_regex_option_b = (
        '^([0-9])____: '
        'observation:obs_type=([A-Za-z]+) '
        'programme:spectrograph:red_Arm:resolution=([a-z]+) '
        'programme:spectrograph:red_Arm:VPH=([A-Z0-9]+) '
        'programme:spectrograph:blue_Arm:resolution=([a-z]+) '
        'programme:spectrograph:blue_Arm:VPH=([A-Z0-9]+)$'
    )
    orb_regex = (
        '^_([0-9])([0-9])([0-9])_: '
        '{([0-9]+) '
        'exposures:exposure:exp_time=([0-9]+) '
        'exposures:exposure:arm=both}$'
    )
    or_regex = (
        '^_([0-9])([0-9])__: '
        '{([0-9]+) '
        'exposures:exposure:exp_time=([0-9]+) '
        'exposures:exposure:arm=red}$'
        )
    ob_regex = (
        '^_([0-9])_([0-9])_: '
        '{([0-9]+) '
        'exposures:exposure:exp_time=([0-9]+) '
        'exposures:exposure:arm=blue}$'
    )
    i_regex = (
        '^____([0-9]): '
        'programme:spectrograph:red_Arm:binning_X=([0-9]+) '
        'programme:spectrograph:blue_Arm:binning_X=([0-9]+)$'
    )
    weave_forbidden_regex = '^WEAVE_FORBIDDEN:( [0-9_]{5})+$'
    pi_forbidden_regex = '^PI_FORBIDDEN:( [0-9_]{5})+$'

    # Create the variables to store the information

    datamver = None

    key_list = ['n', 'or', 'ob', 'i']

    progtemp_dict = {key: {} for key in key_list}

    forbidden_categories = ['weave', 'pi']

    forbidden_dict = {category: [] for category in forbidden_categories}

    # Read the file

    with open(filename, 'r') as f:
        lines = f.readlines()

    # Read each line and extract its info

    for line in lines:

        clean_line = line.strip()

        if _re.match(datamver_regex, clean_line):
            match = _re.match(datamver_regex, clean_line)
            datamver, = match.groups()
            datamver = datamver.strip()
        elif _re.match(header_regex, clean_line):
            pass
        elif (_re.match(n_regex_option_a, clean_line) or
              _re.match(n_regex_option_b, clean_line)):
            if _re.match(n_regex_option_a, clean_line):
                match = _re.match(n_regex_option_a, clean_line)
            elif _re.match(n_regex_option_b, clean_line):
                match = _re.match(n_regex_option_b, clean_line)
            n_char, obsmode, red_resolution, red_vph, blue_resolution, blue_vph = (
                match.groups())
            progtemp_dict['n'][n_char] = {
                'obsmode': obsmode,
                'red_resolution': red_resolution,
                'red_vph': red_vph,
                'blue_resolution': blue_resolution,
                'blue_vph': blue_vph 
            }
        elif _re.match(orb_regex, clean_line):
            pass
        elif _re.match(or_regex, clean_line):
            match = _re.match(or_regex, clean_line)
            o_char, r_char, red_num_exposures, red_exp_time = match.groups()
            progtemp_dict['or'][o_char + r_char] = {
                'red_num_exposures': int(red_num_exposures),
                'red_exp_time': int(red_exp_time)
            }
        elif _re.match(ob_regex, clean_line):
            match = _re.match(ob_regex, clean_line)
            o_char, b_char, blue_num_exposures, blue_exp_time = match.groups()
            progtemp_dict['ob'][o_char + b_char] = {
                'blue_num_exposures': int(blue_num_exposures),
                'blue_exp_time': int(blue_exp_time)
            }
        elif _re.match(i_regex, clean_line):
            match = _re.match(i_regex, clean_line)
            i_char, red_binning_x, blue_binning_x = match.groups()
            progtemp_dict['i'][i_char] = {
                'red_binning_x': int(red_binning_x),
                'blue_binning_x': int(blue_binning_x)
            }
        elif _re.match(weave_forbidden_regex, clean_line):
            forbidden_dict['weave'] = [
                '^{}'.format(substring).replace('_', '[0-9]')
                for substring in clean_line.split()[1:]
            ]
        elif _re.match(pi_forbidden_regex, clean_line):
            forbidden_dict['pi'] = [
                '^{}'.format(substring).replace('_', '[0-9]')
                for substring in clean_line.split()[1:]
            ]
        elif clean_line == '':
            pass
        else:
            raise ValueError

    # Read each line and extract its info

    if assert_orb == True:

        for line in lines:

            clean_line = line.strip()

            if _re.match(orb_regex, clean_line):

                match = _re.match(orb_regex, clean_line)
                o_char, r_char, b_char, num_exposures, exp_time = match.groups()

                red_num_exposures = (
                    progtemp_dict['or'][o_char + r_char]['red_num_exposures'])
                red_exp_time = (
                    progtemp_dict['or'][o_char + r_char]['red_exp_time'])

                assert red_num_exposures == int(num_exposures)
                assert red_exp_time == int(exp_time)

                blue_num_exposures = (
                    progtemp_dict['ob'][o_char + b_char]['blue_num_exposures'])
                blue_exp_time = (
                    progtemp_dict['ob'][o_char + b_char]['blue_exp_time'])

                assert blue_num_exposures == int(num_exposures)
                assert blue_exp_time == int(exp_time)

    return datamver, progtemp_dict, forbidden_dict


def get_progtemp_info(progtemp, progtemp_dict=None, add_datamver=True, cache=False):
    """
    Get the information from a given progtemp value.
 
    Parameters
    ----------
    progtemp : str
        A progtemp value.
    progtemp_dict : dict, optional
        A dictionary with the information to interpret a progtemp value. If
        None, it will retrive the information from Internet.
    add_datamver : bool, optional
        Add a keyword with the datamver value if the progtemp dictionary has
        been retrieved from Internet.
    cache : bool, optional
        If true, will use locally-cached version of progtemp.dat file, rather than
        retrieve from the internet

    
    Returns
    -------
    spectrograph_dict : dict
        A dictionary containing the information of the progtemp value.
    """
    # print(f"[DEBUG] progtemp={repr(progtemp)} ({type(progtemp)})")

    assert _re.match('^[0-9]{5}(\.[0-9]+(\+)?)?$', progtemp)
 
    spectrograph_dict = {}
 
    if progtemp_dict is None:
        
        datamver, progtemp_dict, forbidden_dict = get_progtemp_dict(cache=cache)

        if add_datamver == True:
            spectrograph_dict['datamver'] = datamver

    char_key_dict = {
        'n': progtemp[0],
        'or': progtemp[1] + progtemp[2],
        'ob': progtemp[1] + progtemp[3],
        'i': progtemp[4]
    }

    for char_key in char_key_dict.keys():

        aux_dict = progtemp_dict[char_key][char_key_dict[char_key]]

        for key in aux_dict.keys():
            spectrograph_dict[key] = aux_dict[key]

    return spectrograph_dict


def get_obsmode_from_progtemp(progtemp, progtemp_dict=None):
    """
    Get the obsmode value from a given progtemp.

    Parameters
    ----------
    progtemp : str
        A progtemp value.
    progtemp_dict : dict, optional
        A dictionary with the information to interpret a progtemp value. If
        None, it will retrive the information from Internet.
    
    Returns
    -------
    obsmode : {'MOS', 'LIFU', 'mIFU'}
        A progtemp value.
    """
    spectrograph_dict = get_progtemp_info(progtemp, progtemp_dict=progtemp_dict)

    obsmode = spectrograph_dict['obsmode']

    return obsmode

def get_resolution_from_progtemp(progtemp, progtemp_dict=None, cache=False):
    """
    Get the resolution value from a given progtemp.

    Parameters
    ----------
    progtemp : str
        A progtemp value.
    progtemp_dict : dict, optional
        A dictionary with the information to interpret a progtemp value. If
        None, it will retrive the information from Internet.
    cache : bool, optional
        If true, will use locally-cached version of progtemp.dat file, rather than
        retrieve from the internet

    
    Returns
    -------
    resolution : {'HR', 'LR'}
        A progtemp value.
    """

    spectrograph_dict = get_progtemp_info(progtemp, progtemp_dict=progtemp_dict,cache=cache)
    assert spectrograph_dict['red_resolution'] == spectrograph_dict['blue_resolution']

    _resolution = {'low':"LR",'high':"HR"}

    return _resolution[spectrograph_dict['blue_resolution']]


