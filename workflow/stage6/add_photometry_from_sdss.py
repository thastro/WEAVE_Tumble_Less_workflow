#!/usr/bin/env python3

#
# Copyright (C) 2021 Scott C. Trager and Yago Ascasibar
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
Add SDSS photometry to an IFU FITS catalogue.

Add broadband photometry from SDSS (g,r,i) imaging to an existing IFU catalogue.
For every source in the IFU catalogue, fill the corresponding columns with
circular aperture photometric measurements within each fibre. Note that unlike
the PanSTARRS1 equivalent code by Yago, this code does NOT download SDSS
mosaics, as the web interface for these is quite clunky. See
https://dr12.sdss.org/mosaics for instructions and access to these mosaics.

Authors: sctrager@astro.rug.nl, yago.ascasibar@uam.es
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
from photutils import CircularAperture, aperture_photometry


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


def _nmgy2lups_1band(nmgy, band, _bvalues):
    # <codecell> following is taken from
    # https://github.com/esheldon/sdsspy/blob/master/sdsspy/util.py
    
    _log10 = np.log(10.0)
    # luptitudes = asinh mags
    b=_bvalues[band]
    lups = 2.5*(10.0-np.log10(b)) - 2.5*np.arcsinh(5.0*nmgy/b)/_log10
    # pogson magnitudes
    #lups = 22.5 - 2.5 * np.log10(nmgy)
    return lups


def add_photometry_from_sdss(ifu_cat, output_dir='output'):
    """
    Add SDSS photometry to an IFU FITS catalogue.
    
    Creates a directory where one file per band is saved for every unique TARGID
    in the catalog.

    Parameters
    ----------
    ifu_cat : str
        A FITS file with an IFU catalogue.
    output_dir : str, optional
        Directory where SDSS images will be saved.
    """
    # <codecell> Specific WEAVE (not necessarily LIFU?) code
    
    PS1_PIXEL = 0.25  # arcsec
    LIFU_FIBRE_RADIUS = 1.305  # arcsec (2.6 arcsec diameter)
    LIFU_FIBRE_AREA = np.pi*LIFU_FIBRE_RADIUS**2
    
    _bands = ['u','g','r','i','z']
    _bvalues=[1.4, 0.9, 1.2, 1.8, 7.4]
    
    logging.info(
        'Reading SDSS broadmand imaging for {}'.format(ifu_cat))

    with fits.open(ifu_cat, mode='update') as hdu_list:

        # Get the data from the catalogue

        data = hdu_list[1].data

        for galaxy in np.unique(data['TARGID']):
            rows = np.where(data['TARGID'] == galaxy)
            ra_Gaia = data['GAIA_RA'][rows]
            dec_Gaia = data['GAIA_DEC'][rows]
            coords_Gaia = SkyCoord(ra=ra_Gaia, dec=dec_Gaia,
                                   unit='deg', frame=ICRS())
            coords_SDSS = coords_Gaia.transform_to(FK5(equinox='J2000'))
            ra = coords_SDSS.ra.value
            dec = coords_SDSS.dec.value

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
            size_pix = int(size_deg * 3600. / SDSS_PIXEL)
            logging.info('{}: center=({}, {}); size=({} deg, {} pix)'.format(
                galaxy, center_ra, center_dec, size_deg, size_pix))

            oldphot = {}
            # Photometer images:
            for band in "gir":
                filename = '{}_{}.fits'.format(galaxy, band)
                full_path = os.path.join(output_dir, filename)
                if os.path.isfile(full_path):
                    logging.info('{} exists - good'.format(filename))
                    hdu = fits.open(full_path)
                else:
                    logging.info('Please ensure that SDSS fits file with name' +
                                 ' {} exits'.format(full_path))
                    sys.exit()
 
                wcs_info = WCS(hdu[0].header)
                img = hdu[0].data

                coords_pix = utils.skycoord_to_pixel(coords_SDSS, wcs_info,
                                                     origin=0, mode='all')
                apertures = CircularAperture(np.array(coords_pix).transpose(),
                                             r=LIFU_FIBRE_RADIUS/SDSS_PIXEL)
                phot_table = aperture_photometry(img, apertures)
                bad_photo = np.where(phot_table['aperture_sum'] <= 0)
                photo = _nmgy2lups_1band(phot_table['aperture_sum'].data,
                                         _bands.index(band), _bvalues)

                # Update table:

                colname = 'MAG_'+band.upper()
                oldphot[band] = hdu_list[1].data[colname][rows]
                hdu_list[1].data[colname][rows] = photo
                
                # Scott's addition, 28.11.2020:
                # set core spaxels with r>=25. to TARGUSE="S" -- 
                # anything this faint in PanSTARRS imaging
                # is indistinguishable from sky
                if band=='r':
                    logging.info('Changing TARGUSE, TARGCLASS')
                    targuse=hdu_list[1].data['TARGUSE'][rows]
                    newtarguse=np.where(
                        np.less_equal(hdu_list[1].data[colname][rows],25.) &
                        np.char.equal(targuse,'T'), 'T', 'S')
                    newtarguse2=np.where(
                        np.less_equal(hdu_list[1].data[colname][rows],25.) &
                        np.char.equal(targuse,'R'), 'R', newtarguse)
                    hdu_list[1].data['TARGUSE'][rows]=newtarguse2
                    targclass=hdu_list[1].data['TARGCLASS'][rows]
                    newtargclass=np.where(
                        np.char.equal(newtarguse2,'T') &
                        np.char.not_equal(hdu_list[1].data['IFU_SPAXEL'][0],
                        'S'),'GALAXY',targclass)
                    newtargclass2=np.where(
                        np.char.equal(newtarguse2,'S'),'SKY',newtargclass)
                    hdu_list[1].data['TARGCLASS'][rows]=newtargclass2
                    # clean up photometry and header
                    try:
                        catname = np.where(
                            np.char.equal(newtarguse2,'T'),
                            'Integrated PanSTARRS1 photometry',None)
                        hdu_list[1].data['OPTCAT'][rows]=catname
                        cat='OPTCAT'
                    except:
                        cat=None
                    for magband in "gri":
                        magcolname = 'MAG_'+magband.upper()
                        newphot = np.where(
                            np.char.equal(newtarguse2,'S'), None,
                            hdu_list[1].data[magcolname][rows])
                        newphot2 = np.where(
                            np.char.equal(newtarguse2,'R'), oldphot[magband],
                            newphot)
                        hdu_list[1].data[magcolname][rows] = newphot2
                        if cat=='OPTCAT':
                            stage6_version = "1:{}".format(_get_stage6_version)
                            optcatver = np.where(np.char.equal(newtarguse2,'T'),
                                                 stage6_version, None)
                            hdu_list[1].data[cat+'_'+magcolname][rows] = newphot2
                            hdu_list[1].data['OPTCAT_DR'][rows] = optcatver
                            hdu_list[0].header[magcolname+'_CM']=cat+'_'+magcolname
                        else:
                            try:
                                hdu_list[1].data['SDSS_FIBRE_'+magcolname][rows] = newphot2
                                hdu_list[0].header[magcolname+'_CM']='SDSS_FIBRE_'+magcolname
                            except:
                                try:
                                    hdu_list[1].data['PS_'+magcolname][rows] = newphot2
                                    hdu_list[0].header[magcolname+'_CM']='PS_'+magcolname
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

                # mu=zero_point-2.5*np.log10(img*LIFU_FIBRE_AREA/SDSS_PIXEL**2)
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
                ax = fig.add_axes([.5, 0, .5, 1])  # , projection=wcs_info)
                ax.axis('off')
                sc = ax.scatter(ra, dec, s=40, c=photo,
                                vmin=16, vmax=24, cmap='inferno_r')
                ax.scatter(ra[bad_photo], dec[bad_photo], s=40, c='grey')
                ax.set_xlim(center_ra-.02, center_ra+.02)
                ax.set_ylim(center_dec-.02, center_dec+.02)
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

    parser = argparse.ArgumentParser(
        description='Add SDSS photometry to an IFU FITS catalogue')

    parser.add_argument('ifu_cat',
                        help='a FITS file containing an IFU catalogue')

    parser.add_argument('--outdir', dest='output_dir', default='output',
                        help='directory where SDSS images will be saved')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    add_photometry_from_sdss(args.ifu_cat, output_dir=args.output_dir)

