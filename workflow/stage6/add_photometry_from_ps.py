#!/usr/bin/env python3

#
# Copyright (C) 2021 Yago Ascasibar and Scott C. Trager
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

"""
Download PS1 images and add photometry to an IFU catalogue.

For every source in the IFU catalogue, download Pan-STARRS (g, r, i) images and
fill the corresponding columns with circular aperture photometric measurements
within each fibre.

Authors: yago.ascasibar@uam.es, sctrager@astro.rug.nl
"""

import argparse
import logging
import os

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord, FK5, ICRS
from astropy.wcs import WCS, utils
from photutils.aperture import CircularAperture, aperture_photometry


def _get_stage6_version():
    """
    Get the current stage6 version from git.
    
    This function has been taken from FESWI.

    Returns
    -------
    stage6_version : str
        The version of the current stage6 pipeline from git. Returns empty
        string on failure.
    """

    logging.debug('getting stage6 version from git')

    try:

        # Get the version from a git command

        stage6_version = subprocess.check_output(
            ['git', '-C', os.path.dirname(os.path.realpath(__file__)),
             'describe', '--abbrev=8', '--dirty', '--always', '--tags'])

        # Decode the output if it is not already a string

        if type(stage6_version) is not str:
            stage6_version = stage6_version.decode('utf-8')

    except:
        logging.error('failure guessing the value for CS_VER')
        stage6_version = ''

    # Strip the string (it should contain a new line character at the end)

    stage6_version = stage6_version.strip()

    return stage6_version


def _getimages(ra, dec, size=240, filters="grizy"):
    """Query ps1filenames.py service to get a list of images.

    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """
    # <codecell> Code from PanSTARRS https://ps1images.stsci.edu/ps1image.html
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
           "&filters={filters}").format(**locals())
    table = Table.read(url, format='ascii')
    return table


def _geturl(ra, dec, size=240, output_size=None, filters="grizy",
           format="jpg", color=False):
    """
    Get URL for images in the table.

    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URL for single-filter grayscale images.
    Returns a string with the URL
    """
    if color and format == "fits":
        raise ValueError(
            "color images are available only for jpg or png formats")
    if format not in ("jpg", "png", "fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = _getimages(ra, dec, size=size, filters=filters)
    url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           "ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[np.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0, len(table)//2, len(table)-1]]
        for i, param in enumerate(["red", "green", "blue"]):
            url = url + "&{}={}".format(param, table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
    return url


def add_photometry_from_ps(ifu_cat, targclass='GALAXY', sky=True):
    """
    Download PS1 images and add photometry to an IFU catalogue.
    
    Creates a directory where one file per band is saved for every unique TARGID
    in the catalog.

    Parameters
    ----------
    ifu_cat : str
        A FITS file with an IFU catalogue.
    targclass : str, optional
        TARGCLASS for desired targets.
    """
    # <codecell> Specific WEAVE code

    PS1_PIXEL = 0.25  # arcsec
    LIFU_FIBRE_RADIUS = 1.305  # arcsec (2.6 arcsec diameter)
    LIFU_FIBRE_AREA = np.pi*LIFU_FIBRE_RADIUS**2

    mIFU_FIBRE_RADIUS = 1.305/2.  # arcsec (1.305 arcsec diameter)
    mIFU_FIBRE_AREA = np.pi*mIFU_FIBRE_RADIUS**2

    logging.info(
        'Downloading PS1 broadband imaging for {}'.format(ifu_cat))
    logging.info('Sky fibre conversion? {}'.format(sky))

    output_dir = os.path.splitext(ifu_cat)[0]+'_PS1'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        logging.info('Creating {}'.format(output_dir))

    with fits.open(ifu_cat, mode='update') as hdu_list:

        # Get the data from the catalogue

        data = hdu_list[1].data

        for galaxy in np.unique(data['TARGID']):
#            if galaxy != '':
#                print('galaxy: {}'.format(galaxy))
            if galaxy != '':
                rows = np.where(data['TARGID'] == galaxy)
                ra_Gaia = data['GAIA_RA'][rows]
                dec_Gaia = data['GAIA_DEC'][rows]
                coords_Gaia = SkyCoord(ra=ra_Gaia, dec=dec_Gaia,
                                       unit='deg', frame=ICRS())
                coords_PS1 = coords_Gaia.transform_to(FK5(equinox='J2000'))
                ra = coords_PS1.ra.value
                dec = coords_PS1.dec.value

                try:
                    hdu_list[1].data['PS_DR'][rows]='2'
                except:
                    try:
                        hdu_list[1].data['SDSS_DR'][rows]='2'
                    except:
                        raise ValueError('Template catalogue does not contain photometry information')
                # Center and size of the images:
                ra_min = np.min(ra)
                ra_max = np.max(ra)
                dec_min = np.min(dec)
                dec_max = np.max(dec)
                center_ra = (ra_min + ra_max)/2
                center_dec = (dec_min + dec_max)/2
                delta_ra = ra_max - ra_min
                delta_dec = dec_max - dec_min
                size_deg = np.max([delta_ra*np.sin(np.radians(center_dec)),
                                   delta_dec]) * 1.1  # add 10% buffer
                size_pix = int(size_deg * 3600. / PS1_PIXEL)
                logging.info('{}: center=({}, {}); size=({} deg, {} pix)'.format(
                    galaxy, center_ra, center_dec, size_deg, size_pix))

                oldphot = {}
                # Retrieve images:
                for band in "gir":
                    filename = '{}_{}.fits'.format(galaxy, band)
                    full_path = os.path.join(output_dir, filename)
                    if os.path.isfile(full_path):
                        logging.info('{} already exists'.format(filename))
                        hdu = fits.open(full_path)
                    else:
                        logging.info('Downloading {}'.format(filename))
                        fitsurl = _geturl(center_ra, center_dec, size=size_pix,
                                         filters=band, format="fits")
                        for url in fitsurl:
                            #band = url[-13]
                            # logging.debug(filename, band)
                            hdu = fits.open(url)
                            hdu.writeto(full_path, overwrite=True)

                    zero_point = 25 + 2.5*np.log10(hdu[0].header['exptime'])
                    wcs_info = WCS(hdu[0].header)
                    img = hdu[0].data

                    coords_pix = utils.skycoord_to_pixel(coords_PS1, wcs_info,
                                                         origin=0, mode='all')

                    # Configure>=3.3.7 writes out the six "dead" fibres into the first six rows; skip these
                    if hdu_list[1].data['IFU_SPAXEL'][6][0]!='m': #LIFU mode
    #                if hdu_list[1].data['IFU_SPAXEL'][0][0]!='m': #LIFU mode
                        logging.debug('LIFU mode')
                        radius=LIFU_FIBRE_RADIUS/PS1_PIXEL
                    else: # mIFU mode
                        logging.debug('mIFU mode')
                        radius=mIFU_FIBRE_RADIUS/PS1_PIXEL

                    logging.debug('radius={} pixels'.format(radius))

                    apertures = CircularAperture(np.array(coords_pix).transpose(),
                                                 r=radius)
                    phot_table = aperture_photometry(img, apertures)
                    # good_photo = np.where(phot_table['aperture_sum'] > 0)
                    bad_photo = np.where(phot_table['aperture_sum'] <= 0)
                    photo = zero_point - 2.5*np.log10(
                        phot_table['aperture_sum'].data)

                    # Update table:

                    colname = 'MAG_'+band.upper()
                    oldphot[band] = hdu_list[1].data[colname][rows]
                    hdu_list[1].data[colname][rows] = photo
                    # Scott's addition, 28.11.2020:
                    # set core spaxels with r>=25. to TARGUSE="S" -- 
                    # anything this faint in PanSTARRS imaging
                    # is indistinguishable from sky
                    if band=='r':
                        targuse=hdu_list[1].data['TARGUSE'][rows]
                        if sky:
                            # set TARGUSE properly -- I'm sure there's a cleverer way to do this!
                            logging.info('Changing TARGUSE, TARGCLASS')
                            newtarguse0=np.where(np.less_equal(hdu_list[1].data[colname][rows],25.) &
                                                 np.char.equal(targuse,'T'), 'T', 'S')
                        else:
                            logging.info('no targuse="S"')
                            newtarguse0=targuse
                        newtarguse1=np.where(np.char.equal(targuse,'R'), 'R', newtarguse0)
                        newtarguse2=np.where(np.char.equal(targuse,'C'), 'C', newtarguse1)
                        hdu_list[1].data['TARGUSE'][rows]=newtarguse2
                        # now set TARCLASS properly
                        targclass_aux = hdu_list[1].data['TARGCLASS'][rows]
                        newtargclass=np.where(np.char.equal(newtarguse2,'T') &
                                              np.char.not_equal(hdu_list[1].data['IFU_SPAXEL'][0],'S'),
                                              targclass,targclass_aux)
                        newtargclass2=np.where(np.char.equal(newtarguse2,'S'),'SKY',newtargclass)
                        hdu_list[1].data['TARGCLASS'][rows]=newtargclass2
                        # clean up photometry and header
                        try:
                            catname = np.where(np.char.equal(newtarguse2,'T'),
                                               'Integrated PanSTARRS1 photometry','')
                            hdu_list[1].data['OPTCAT'][rows]=catname
                            cat='OPTCAT'
                        except:
                            cat=None
                        for magband in "gri":
                            magcolname = 'MAG_'+magband.upper()
                            newphot0 = np.where(np.char.equal(newtarguse2,'S'),
                                                None, hdu_list[1].data[magcolname][rows])
                            newphot2 = np.where(np.logical_or(np.char.equal(newtarguse2,'R'),
                                                              np.char.equal(newtarguse2,'C')),
                                                oldphot[magband], newphot0)
                            hdu_list[1].data[magcolname][rows] = newphot2
                            if cat=='OPTCAT':
                                stage6_version = "1-{}".format(_get_stage6_version())
                                optcatver = np.where(np.char.equal(newtarguse2,'T'),
                                                     stage6_version, '')
                                hdu_list[1].data[cat+'_'+magcolname][rows] = newphot2
                                hdu_list[1].data['OPTCAT_DR'][rows] = optcatver
                                hdu_list[0].header[magcolname+'_CM']=cat+'_'+magcolname

                            try:
                                hdu_list[1].data['PS_'+magcolname][rows] = newphot2
                                hdu_list[0].header[magcolname+'_CM']='PS_'+magcolname
                            except:
                                try:
                                    hdu_list[1].data['SDSS_'+magcolname][rows] = newphot2
                                    hdu_list[0].header[magcolname+'_CM']='SDSS_'+magcolname
                                except:
                                    try:
                                        hdu_list[1].data['SDSS_FIBRE_'+magcolname][rows] = newphot2
                                        hdu_list[0].header[magcolname+'_CM']='SDSS_FIBRE_'+magcolname
                                    except:
                                        pass

                    # xx = hdu_list[1].data[colname]
                    # logging.debug('\n\n')
                    # logging.debug(len(rows[0]), rows)
                    # # logging.debug(len(good_photo[0]), good_photo)
                    # logging.debug(len(photo), photo)
                    # # logging.debug(len(xx), len(xx[rows][good_photo]), xx[rows][good_photo])
                    # logging.debug(xx[rows][0:10])

                    # Diagnostic plots:

                    # mu=zero_point-2.5*np.log10(img*LIFU_FIBRE_AREA/PS1_PIXEL**2)
                    # # plt.figure()
                    # plt.figure(figsize=(8, 8))
                    # plt.imshow(mu, origin='lower',
                    #            vmin=16, vmax=24, cmap='inferno_r')
                    # # plt.xlim(size_pix/2-100, size_pix/2+100)
                    # # plt.ylim(size_pix/2-100, size_pix/2+100)
                    # plt.colorbar()
                    # plt.show()
                    # plt.figure(figsize=(8, 8))
                    # # plt.imshow(mu, origin='lower',
                    #             # vmin=0, vmax=2, cmap='inferno_r')
                    # # plt.scatter(coords_pix[0], coords_pix[1], s=10)
                    # plt.scatter(coords_pix[0][good_photo],
                    #             coords_pix[1][good_photo], s=100,
                    #             c=photo, vmin=16, vmax=24, cmap='inferno_r')
                    # plt.xlim(size_pix/2-150, size_pix/2+150)
                    # plt.ylim(size_pix/2-150, size_pix/2+150)
                    # plt.colorbar().set_label('fibre '+band+' mag')
                    # plt.title(filename)
                    # plt.show()
                    fig = plt.figure(figsize=(12, 6))
                    ax = fig.add_axes([.1, .1, .4, .8])  # , projection=wcs_info)
                    sc = ax.scatter(ra, dec, s=4, c=photo,
                                    vmin=16, vmax=24, cmap='inferno_r')
                    ax.scatter(ra[bad_photo], dec[bad_photo], s=4, c='grey')
                    ax.grid('both')
                    ax.set_xlabel('RA (J2000)')
                    ax.set_ylabel('DEC (J2000)')
                    ax = fig.add_axes([.55, .1, .4, .8])  # , projection=wcs_info)
                    ax.axis('off')
                    sc = ax.scatter(ra, dec, s=40, c=photo,
                                    vmin=16, vmax=24, cmap='inferno_r')
                    ax.scatter(ra[bad_photo], dec[bad_photo], s=40, c='grey')
                    #ax.set_xlim(center_ra-.02, center_ra+.02)
                    #ax.set_ylim(center_dec-.02, center_dec+.02)
                    fig.suptitle(filename)
                    ax = fig.add_axes([.55, .1, .4, .025])
                    cb = fig.colorbar(sc, orientation="horizontal", cax=ax)
                    cb.set_label('fibre '+band+' mag')
                    plt.savefig(full_path[:-4]+'png')
                    # plt.show()
                    plt.close()

        # g = hdu_list[1].data['MAG_G']
        # logging.debug(len(rows[0]), rows)
        # logging.debug(len(good_photo[0]), good_photo)
        # logging.debug(len(photo), photo)
        # logging.debug(len(g), len(g[rows][good_photo]), g[rows][good_photo])

        logging.info('Saving {}'.format(ifu_cat))
        hdu_list.flush()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=
        'Download PS1 images and add photometry to an IFU catalogue')

    parser.add_argument('ifu_cat',
                        help='a FITS file containing an IFU catalogue')

    parser.add_argument('--targclass', default='GALAXY',
                        choices=['GALAXY', 'UNKNOWN', 'NEBULA', 'QSO', 'SKY',
                                 'STAR', 'STAR_BHB', 'STAR_CEP', 'STAR_EM',
                                 'STAR_EMP', 'STAR_FGK', 'STAR_IB', 'STAR_MLT',
                                 'STAR_MLUM', 'STAR_OB', 'STAR_BA', 'STAR_RRL',
                                 'STAR_VAR', 'STAR_WD', 'STAR_YSO'],
                        help='TARGCLASS for desired targets')

    parser.add_argument('--sky', action='store_false',
                        help='do not set spaxels with r>=25. to TARGUSE="S"')
    
    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    add_photometry_from_ps(args.ifu_cat, targclass=args.targclass, sky=args.sky)

