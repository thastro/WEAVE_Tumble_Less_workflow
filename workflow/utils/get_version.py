#!/usr/bin/env python3

#
# Copyright (C) 2022 Cambridge Astronomical Survey Unit
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
from subprocess import check_output as _check_output


def get_version(path=None):
    """
    Get version of a git repository.
    
    Get version of a git repository if a directory is provided, otherwise, get
    the version corresponding to directory of this file.
    
    Returns
    -------
    result : str
        A string containing the version of the repository.
    """
    
    if path is None:
        git_repo_path = _os.path.dirname(_os.path.realpath(__file__))
    else:
        git_repo_path = path

    try:
        result = _check_output(
            ['git', 'describe', '--abbrev=8', '--dirty', '--always', '--tags'],
            cwd=git_repo_path).strip().decode('utf-8')
    except:
        logging.error(
            'error while trying to get the version of the repository ' +
            '{}: Have you missed to install git? '.format(git_repo_path) +
            'Does this directory contain a git repository?')
        
        raise Exception
    
    return result

