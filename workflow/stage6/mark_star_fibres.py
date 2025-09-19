#!/usr/bin/env python3

#
# Copyright (C) 2021 Scott C. Trager
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
Mark Gaia eDR3 stars in a IFU FITS catalogue.

For a given IFU FITS catalogue, query the Gaia eDR3 database for stars near
individual fibres and then mark these fibres as TARGCLASS="STAR" and enter the
GAIA_MAG_X into the catalogues.

Author: sctrager@astro.rug.nl
"""

import argparse
import logging
import os

import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord, ICRS
import astropy.units as u
from astroquery.gaia import Gaia


def mark_star_fibres(ifu_cat, match_dist=2.5, parallax_over_error=3.0):
    """
    Mark Gaia eDR3 stars in a IFU FITS catalogue.
    
    For a given IFU FITS catalogue, query the Gaia eDR3 database for stars near
    individual fibres and then mark these fibres as TARGCLASS="STAR" and enter
    the GAIA_MAG_X into the catalogues.

    Parameters
    ----------
    ifu_cat : str
        Name of a FITS file containing a IFU catalogue.
    match_dist : float, optional
        The distance in (m/L)IFU fibre radii to search for nearby Gaia stars
        (in units of IFU fibre radii).
    parallax_over_error : float, optional
        Gaia sources are assumed to be real when parallax/error >
        parallax_over_error.
    """

    # check whether Gaia archive is accessable
    try:
        Gaia.ROW_LIMIT = 1  # Ensure the default row limit.
        coord = SkyCoord(ra=280, dec=-60, unit=(u.degree, u.degree), frame='icrs')
        radius = u.Quantity(0.01, u.deg)
        j = Gaia.cone_search_async(coord, radius)
    except:
        logging.error('Gaia server inaccessable; shutting down')
        import sys
        sys.exit()
    
    LIFU_FIBRE_RADIUS = 1.305 * u.arcsec  # arcsec (2.61 arcsec diameter)
    LIFU_LONG_AXIS = 90.39 * u.arcsec
    mIFU_FIBRE_RADIUS = 1.305/2. * u.arcsec # arcsec (1.305 arcsec diameter)
    mIFU_LONG_AXIS = 12.35 * u.arcsec
    
    logging.info(
        'Searching for Gaia eDR stars near fibres in {}'.format(ifu_cat))

    with fits.open(ifu_cat, mode='update') as hdu_list:

        Gaia_stars=Table()
        # loop over each target in the catalogue
        for galaxy in np.unique(hdu_list[1].data['TARGID']):
            print(galaxy)
            if galaxy != '':
                rows=np.where(hdu_list[1].data['TARGID'] == galaxy)
                logging.info(
                    'Searching for stars in the field of {}'.format(galaxy))

                if np.any(
                    np.char.startswith(hdu_list[1].data['IFU_SPAXEL'][rows],'m')):
                    ifu='mIFU'
                else:
                    ifu='LIFU'
                logging.debug('IFU mode: {}'.format(ifu))
                ra=hdu_list[1].data['GAIA_RA'][rows]
                dec=hdu_list[1].data['GAIA_DEC'][rows]
                central_ra=np.median(ra)
                central_dec=np.median(dec)
                logging.info('Central fibre of {} located at {},{}'.format(
                    galaxy, central_ra, central_dec))

                central_coord = SkyCoord(ra=central_ra, dec=central_dec,
                                         unit=(u.degree,u.degree), frame='icrs')
                # search radius for ALL stars in field: use IFU size + 25% for
                # dithers match radius is match_dist*fibre radius
                if ifu=='LIFU':
                    radius = 1.25 * LIFU_LONG_AXIS / 2.
                    match_radius=LIFU_FIBRE_RADIUS * match_dist
                elif ifu=='mIFU':
                    radius = 1.25 * mIFU_LONG_AXIS / 2.
                    match_radius=mIFU_FIBRE_RADIUS * match_dist
                else:
                    logging.error('wrong IFU mode, exiting')
                    sys.exit()
                # serach Gaia eDR3
                cone_search = Gaia.cone_search_async(
                    central_coord, radius, table_name='gaiaedr3.gaia_source')
                Gaia_sources = cone_search.get_results()
                logging.debug(Gaia_sources.pprint())
                Gaia_stars = vstack(
                    [Gaia_stars,
                     Gaia_sources[Gaia_sources['parallax_over_error'] >
                     parallax_over_error]])

        if len(Gaia_stars) == 0:
            logging.info('No appropriate Gaia stars present; exiting')
            hdu_list.flush()
            import sys ; sys.exit()
        logging.info(Gaia_stars.pprint())
        Gaia_stars.add_index('source_id')
        logging.info('Building coordinate masks')
        # construct star coordinates
        Gaia_star_coords = SkyCoord(ra=Gaia_stars['ra'], dec=Gaia_stars['dec'],
                                    unit=(u.degree,u.degree), frame='icrs')
        # construct fibre coordinates
        fibre_coords = SkyCoord(ra=hdu_list[1].data['GAIA_RA'],
                                dec=hdu_list[1].data['GAIA_DEC'],
                                unit=(u.degree,u.degree), frame='icrs')
        # compute separation of each fibre from each star found
        fibre_separations = np.array(
            [fibre_coords.separation(Gaia_star_coord)
             for Gaia_star_coord in Gaia_star_coords]) * u.degree
        # is there a star within match_radis?
        nearby_star = np.logical_or.reduce(np.less_equal(fibre_separations,
                                                         match_radius))
        # which is the nearest star (even if it isn't within match_radius)?
        which_nearby_star = np.argmin(fibre_separations,axis=0)

        # change TARGCLASS and TARGNAME
        mask=np.logical_and(nearby_star,
            np.logical_and(np.char.not_equal(hdu_list[1].data['TARGUSE'],'R'),
            np.char.not_equal(hdu_list[1].data['TARGUSE'],'C')))
        mask=np.logical_and(mask,
            np.char.not_equal(hdu_list[1].data['TARGSRVY'],'WD'))
        newtarguse = np.where(mask,'T',hdu_list[1].data['TARGUSE'])
        hdu_list[1].data['TARGUSE']=newtarguse
        newtargclass = np.where(mask,'STAR',hdu_list[1].data['TARGCLASS'])
        hdu_list[1].data['TARGCLASS']=newtargclass
        logging.info('Updated TARGCLASS')
        # get Gaia photometry for each star
        gaia_g = np.where(nearby_star,
                          Gaia_stars['phot_g_mean_mag'][which_nearby_star],
                          None)
        gaia_bp = np.where(nearby_star,
                          Gaia_stars['phot_bp_mean_mag'][which_nearby_star],
                          None)
        gaia_rp = np.where(nearby_star,
                          Gaia_stars['phot_rp_mean_mag'][which_nearby_star],
                          None)
        hdu_list[1].data['GAIA_MAG_G'] = gaia_g
        hdu_list[1].data['GAIA_MAG_BP'] = gaia_bp
        hdu_list[1].data['GAIA_MAG_RP'] = gaia_rp
        logging.info('Updated Gaia photometry')

        logging.info('Saving {}'.format(ifu_cat))
        hdu_list.flush()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Mark Gaia eDR3 stars in a IFU FITS catalogue')

    parser.add_argument('ifu_cat',
                        help='a FITS file containing an IFU catalogue')

    parser.add_argument('--match_dist', default=2.5,
                        help="""Star matching distance in units of IFU fibre
                        radii""")

    parser.add_argument('--parallax_over_error', default=3.0,
                        help="""Gaia sources are assumed to be real when
                        parallax/error > parallax_over_error""")

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    mark_star_fibres(args.ifu_cat, match_dist=args.match_dist,
                     parallax_over_error=args.parallax_over_error)

