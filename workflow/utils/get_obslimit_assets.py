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

from workflow.utils.get_resources import get_obslimits_tool as _get_obslimits_tool
from workflow.utils.get_resources import get_configure_cfg as _get_configure_cfg



def get_obslimits_tool(filename=None, cache=False):
    """
    Function to retrieve the obs_limits.py file and provide a pointer to it

    Parameters
    ----------
    filename : str
        The name of the file formatted like obs_limits.py.
    
    Returns
    -------
    filepath : str
        The location on disk of the .py file

    """
    

    
    if (filename is None) and (cache):
        # set the filename, ensure it exists
        # filename = 'aux/obs_limits.py'
        filename = _os.path.join(_os.path.dirname(_os.path.dirname(__file__)),'aux/obs_limits.py')
        if not _os.path.isfile(filename):
            logging.warning("User set cache=True but file {} does not exist - will download".format(filename))
            _get_obslimits_tool(file_path=filename)
        else:
            logging.warning("Using locally-cached copy of obs_limits.py file at {} - we advise disabling cache mode to ensure this file is up to date".format(
                filename))
        
    elif filename is None:
        filename = _os.path.join(_tempfile.mkdtemp(), 'obs_limits.py')
        _get_obslimits_tool(file_path=filename)

    return filename

def get_configure_cfg(filename=None, assert_orb=True, cache=False):
    """
    Function to retrieve the configure.xml file and provide a pointer to it

    Parameters
    ----------
    filename : str
        The name of the file formatted like configure.xml.
    
    Returns
    -------
    filepath : str
        The location on disk of the .xml file

    """
    
    if (filename is None) and (cache):
        # set the filename, ensure it exists
        # filename = 'aux/configure.xml'
        filename = _os.path.join(_os.path.dirname(_os.path.dirname(__file__)),'aux/configure.xml')
        if not _os.path.isfile(filename):
            logging.warning("User set cache=True but file {} does not exist - will download".format(filename))
            _get_configure_cfg(file_path=filename)
        else:
            logging.warning("Using locally-cached copy of configure.xml file at {} - we advise disabling cache mode to ensure this file is up to date".format(
            filename))
        
    elif filename is None:
        filename = _os.path.join(_tempfile.mkdtemp(), 'configure.xml')
        _get_configure_cfg(file_path=filename)

    return filename




