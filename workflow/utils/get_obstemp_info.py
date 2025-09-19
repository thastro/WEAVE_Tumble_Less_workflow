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

from workflow.utils.get_resources import get_obstemp_file as _get_obstemp_file


def get_obstemp_dict(filename=None, cache=False):
    """
    Get a dictionary to interpret obstemp values.

    Parameters
    ----------
    filename : str
        The name of the file formatted like obstemp.dat.
    cache : bool, optional
        If true, will use locally-cached version of obsemp.dat file, rather than
        retrieve from the internet
    
    Returns
    -------
    datamver : str
        The value of DATAMVER in the file.
    obstemp_dict : dict
        A dictionary with the information to interpret an obstemp value.
    """

    if (filename is None) and (cache):
        # set the filename, ensure it exists
        filename = _os.path.join(_os.path.dirname(_os.path.dirname(__file__)),'aux/obstemp.dat')

        if not _os.path.isfile(filename):
            logging.warning("User set cache=True but file {} does not exist - will download".format(filename))
            _get_obstemp_file(file_path=filename)
        else:
            logging.warning("Using locally-cached copy of OBSTEMP file at {} - we advise disabling cache mode to ensure this file is up to date".format(
                filename))
    elif filename is None:
        filename = _os.path.join(_tempfile.mkdtemp(), 'obstemp.dat')
        _get_obstemp_file(file_path=filename)

    # Set of the regex expected to be found in the file

    datamver_regex = '^#DATAMVER (.+)$'
    header_regex = '^[A-Z_]{5}:$'
    seeing_max_regex = \
        '^([A-Z])____: observation:obsconstraints:seeing_max=([0-9.]+)$'
    transparency_min_regex = \
        '^_([A-Z])___: observation:obsconstraints:transparency_min=([0-9.]+)$'
    elevation_min_regex = \
        '^__([A-Z])__: observation:obsconstraints:elevation_min=([0-9.]+) ## airmass ([0-9.]+)$'
    moondist_min_regex = \
        '^___([A-Z])_: observation:obsconstraints:moondist_min=([0-9.]+)$'
    skybright_max_regex = \
        '^____([A-Z]): observation:obsconstraints:skybright_max=([0-9.]+)$'

    # Create the variables to store the information

    datamver = None

    key_list = [
        'seeing_max', 'transparency_min', 'elevation_min', 'moondist_min',
        'skybright_max'
    ]

    obstemp_dict = {key: {} for key in key_list}

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
        elif _re.match(seeing_max_regex, clean_line):
            match = _re.match(seeing_max_regex, clean_line)
            key, value = match.groups()
            obstemp_dict['seeing_max'][key] = value
        elif _re.match(transparency_min_regex, clean_line):
            match = _re.match(transparency_min_regex, clean_line)
            key, value = match.groups()
            obstemp_dict['transparency_min'][key] = value
        elif _re.match(elevation_min_regex, clean_line):
            match = _re.match(elevation_min_regex, clean_line)
            key, value, airmass = match.groups()
            obstemp_dict['elevation_min'][key] = value
        elif _re.match(moondist_min_regex, clean_line):
            match = _re.match(moondist_min_regex, clean_line)
            key, value = match.groups()
            obstemp_dict['moondist_min'][key] = value
        elif _re.match(skybright_max_regex, clean_line):
            match = _re.match(skybright_max_regex, clean_line)
            key, value = match.groups()
            obstemp_dict['skybright_max'][key] = value
        elif clean_line == '':
            pass
        else:
            raise ValueError

    return datamver, obstemp_dict


def get_obstemp_info(obstemp, obstemp_dict=None, add_datamver=True):
    """
    Get the information from a given obstemp value.

    Parameters
    ----------
    obstemp : str
        An obstemp value.
    obstemp_dict : dict, optional
        A dictionary with the information to interpret an obstemp value.
        If None, it will retrive the information from Internet.
    add_datamver : bool, optional
        Add a keyword with the datamver value if the obstemp dictionary has been
        retrieved from Internet.
    
    Returns
    -------
    obsconstraints_dict : dict
        A dictionary containing the information of the observation contraints.
    """
    
    assert _re.match('^[A-Z]{5}$', obstemp)

    obsconstraints_dict = {}

    if obstemp_dict is None:

        datamver, obstemp_dict = get_obstemp_dict()

        if add_datamver == True:
            obsconstraints_dict['datamver'] = datamver

    key_list = [
        'seeing_max', 'transparency_min', 'elevation_min', 'moondist_min',
        'skybright_max'
    ]
        
    for i, key in enumerate(key_list):
        obsconstraints_dict[key] = obstemp_dict[key][obstemp[i]]

    return obsconstraints_dict

