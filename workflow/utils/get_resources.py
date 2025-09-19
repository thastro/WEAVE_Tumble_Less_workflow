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


import os as _os
import urllib.request
import urllib.error
import urllib.parse
import urllib as _urllib
import logging

def get_master_cat(file_path='Master_CatalogueTemplate.fits',
                   url=('http://wasp.ast.cam.ac.uk/data_model/cats/weave/' +
                        'Master_CatalogueTemplate.fits')):
    """
    Download the latest version of the master catalogue.

    Parameters
    ----------
    file_path : str, optional
        The path used to save the downloaded file.
    url : str, optional
        The URL with the location to the file in Internet.

    Returns
    -------
    file_path : str
        The path used to save the downloaded file.
    """
    
    _urllib.request.urlretrieve(url, file_path)
    
    return file_path


def get_aladin_jar(file_path='Aladin.jar',
                   url='https://aladin.u-strasbg.fr/java/Aladin.jar'):
    """
    Download the latest version of Aladin.

    Parameters
    ----------
    file_path : str, optional
        The path used to save the downloaded file.
    url : str, optional
        The URL with the location to the file in Internet.

    Returns
    -------
    file_path : str
        The path used to save the downloaded file.
    """
    
    _urllib.request.urlretrieve(url, file_path)
    
    return file_path


def get_blank_xml_template(file_path='BlankXMLTemplate.xml',
                           url=('http://wasp.ast.cam.ac.uk/data_model/' +
                                'misc/BlankXMLTemplate.xml'), cache=False):
    """
    Download the latest version of the blank XML template.

    Parameters
    ----------
    file_path : str, optional
        The path used to save the downloaded file.
    url : str, optional
        The URL with the location to the file in Internet.

    Returns
    -------
    file_path : str
        The path used to save the downloaded file.
    """
    
    if (cache) and (_os.path.isfile(file_path)):
        return file_path
    _urllib.request.urlretrieve(url, file_path)
    
    return file_path


def get_progtemp_file(file_path='progtemp.dat',
                      url=('http://wasp.ast.cam.ac.uk/data_model/misc/' +
                           'progtemp.dat'), cache=True):
    """
    Download the latest version of progtemp.dat.

    Parameters
    ----------
    file_path : str, optional
        The path used to save the downloaded file.
    url : str, optional
        The URL with the location to the file in Internet.

    Returns
    -------
    file_path : str
        The path used to save the downloaded file.
    """
    
    if (cache) and (_os.path.isfile(file_path)):
        return file_path
    
    _urllib.request.urlretrieve(url, file_path)
    
    return file_path


def get_obstemp_file(file_path='obstemp.dat',
                     url=('http://wasp.ast.cam.ac.uk/data_model/misc/' +
                          'obstemp.dat'), cache=True):

    """
    Download the latest version of obstemp.dat.

    Parameters
    ----------
    file_path : str, optional
        The path used to save the downloaded file.
    url : str, optional
        The URL with the location to the file in Internet.

    Returns
    -------
    file_path : str
        The path used to save the downloaded file.
    """
    
    if (cache) and (_os.path.isfile(file_path)):
        return file_path
    _urllib.request.urlretrieve(url, file_path)
    
    return file_path

def get_configure_cfg(file_path='configure.xml',
                     url=('http://wasp.ast.cam.ac.uk/data_model/misc/' +
                          'configure.xml')):
    """
    Download the latest version of configure.xml.

    Parameters
    ----------
    file_path : str, optional
        The path used to save the downloaded file.
    url : str, optional
        The URL with the location to the file in Internet.

    Returns
    -------
    file_path : str
        The path used to save the downloaded file.
    """
    try:
        _urllib.request.urlretrieve(url, file_path)
    except urllib.error.URLError as ex:
        logging.error("Cannot access URL - are you online?")
        raise SystemExit(0)

    return file_path


def get_obslimits_tool(file_path='obs_limits.py',
                     url=('http://wasp.ast.cam.ac.uk/data_model/tools/' +
                          'obs_limits.py')):
    """
    Download the latest version of observability tool.

    Parameters
    ----------
    file_path : str, optional
        The path used to save the downloaded file.
    url : str, optional
        The URL with the location to the file in Internet.

    Returns
    -------
    file_path : str
        The path used to save the downloaded file.
    """
    
    try:
        _urllib.request.urlretrieve(url, file_path)
    except urllib.error.URLError as ex:
        logging.error("Cannot access URL - are you online?")
        raise SystemExit(0)

    return file_path




def _get_srvy_cat(survey, healpix_index, directory='',
                  server='http://wasp.ast.cam.ac.uk/',
                  clobber=False):
    """
    Download a catalogue of stars for the requested survey and HEALPix index.

    Parameters
    ----------
    survey : {'WD', 'WG'}
        The survey to be considered.
    healpix_index : int
        The desired HEALPix index.
    directory : str, optional
        The directory used to save the downloaded file.
    server : str, optional
        The URL of the server hosting WEAVE Operational Catalogues
    clobber : If the file already exists, should it be overwritten?

    Returns
    -------
    file_path : str
        The path used to save the downloaded file.
    """

    assert survey in ['WD', 'WG']
    filename = '{}_{}.fits'.format(survey, healpix_index)
    file_path = _os.path.join(directory, filename)

    if (_os.path.isfile(file_path)) and (not clobber):
        logging.info("Using cached copy of HEALPix ID {} for survey '{}'".format(
            healpix_index,survey))
        return file_path

    url = "{}{}/getHPID".format(_os.path.join(server,''),survey.lower())
    data = _urllib.parse.urlencode({'id': healpix_index}).encode('ascii')
    logging.info("Downloading copy of HEALPix ID {} for survey '{}'".format(
        healpix_index,survey))

    try:
        response = _urllib.request.urlopen(url, data=data)
    except _urllib.error.HTTPError as ex:
        if ex.code == 404:
            logging.error("HEALPix ID {} is not available for survey '{}'".format(
                healpix_index,survey))
        elif ex.code == 500:
            logging.error("There was an application error in server {} (HTTP/500 returned)".format(
                server))
        elif ex.code == 502:
            logging.error("Server {} is either over-capacity or the query function failed".format(
                server))
        elif ex.code == 503:
            logging.error("Server {} could not be reached".format(
                server))
        elif ex.code == 403:
            logging.error('Server {} did not grant permission to access Operational catalogues from your IP'.format(
                server))
            logging.error('Please ensure your IP is registered in the WASP IP whilelist: consult your STL or PI')
        else:
            print('The server {} returned with HTTP error {}'.format(
                server,ex.code))
        raise SystemExit(0)
    except urllib.error.URLError as ex:
        logging.error("Cannot access URL - are you online?")
        raise SystemExit(0)

    with open(file_path, 'wb') as f:
        contents = response.read()
        f.write(contents)
    
    return file_path


def get_guide_cat(healpix_index, directory='',
                  server='http://wasp.ast.cam.ac.uk/'):
    """
    Download a catalogue of guide stars for the requested HEALPix index.

    Parameters
    ----------
    healpix_index : int
        The desired HEALPix index.
    directory : str, optional
        The directory used to save the downloaded file.
    server : str, optional
        The URL of the server hosting WEAVE Operational Catalogues

    Returns
    -------
    file_path : str
        The path used to save the downloaded file.
    """

    file_path = _get_srvy_cat('WG', healpix_index, directory=directory, server=server)
    return file_path


def get_calib_cat(healpix_index, directory='',
                  server='http://wasp.ast.cam.ac.uk/'):
    """
    Download a catalogue of calib stars for the requested HEALPix index.

    Parameters
    ----------
    healpix_index : int
        The desired HEALPix index.
    directory : str, optional
        The directory used to save the downloaded file.
    server : str, optional
        The URL of the server hosting WEAVE Operational Catalogues

    Returns
    -------
    file_path : str
        The path used to save the downloaded file.
    """
    
    file_path = _get_srvy_cat('WD', healpix_index, directory=directory, server=server)
    
    return file_path

