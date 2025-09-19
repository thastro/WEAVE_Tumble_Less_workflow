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


"""
Retrieve of WEAVE guide and calib stars.

The authors of this module are:

- David Murphy (dmurphy@ast.cam.ac.uk),
  Cambridge Astronomical Survey Unit (CASU, IoA).
- Luis Peralta de Arriba (lperalta@ast.cam.ac.uk),
  Cambridge Astronomical Survey Unit (CASU, IoA).

Notes
----- 
The dependencies of this module are:

- Python core
- numpy
- matplotlib
- astropy
- astropy_healpix
"""


import logging
import tempfile as _tempfile

import numpy as _np
import sys
import os
from copy import deepcopy
from astropy import coordinates as _coordinates
from astropy.io import fits as _fits
from astropy.table import Column as _Column
from astropy.table import Table as _Table
from astropy.table import vstack as _vstack
from astropy_healpix import HEALPix as _HEALPix
import matplotlib as _mpl
_mpl.use('Agg')
import matplotlib.pyplot as _plt
from matplotlib.patches import Ellipse as _Ellipse

from workflow.utils.get_resources import get_guide_cat as _get_guide_cat
from workflow.utils.get_resources import get_calib_cat as _get_calib_cat
from workflow.utils.get_obslimit_assets import (get_obslimits_tool,
                                               get_configure_cfg)



class DefaultDir:
    def __enter__(self):
        _file = __file__.replace('/./','/')
        loc = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(_file))),'stage3/input/')
        return loc

    def __exit__(self,exc_type, exc_val, exc_tb):
        return

class _AuxStars:


    def __init__(self, ra, dec, obsmode, resolution, dither_size_arcsec=10.0, nside=32):

        self.star_type = ''
        self.ra = ra
        self.dec = dec
        self.obsmode = obsmode
        self.resolution = resolution
        self.dither_size_arcsec = dither_size_arcsec
        self.nside = nside
        
        self._assert_args()

        self._set_geometry()

        self.healpix_indices = None
        self.full_table = None
        self.useful_table = None
        self.useful_maglim_table = None
        self.selected_pa = None
        self.selected_table = None

        
    def _assert_args(self):
        
        if not ((self.ra >= 0.0) and (self.ra < 360.0)):
            logging.error('invalid ra value: {}'.format(self.ra))
            raise ValueError
        
        if not ((self.dec >= -90.0) and (self.dec < 90.0)):
            logging.error('invalid dec value: {}'.format(self.dec))
            raise ValueError

        if self.obsmode not in ['MOS', 'LIFU', 'mIFU']:
            logging.error('invalid obsmode: {}'.format(obsmode))
            raise ValueError
        
        
    def _set_geometry(self):

        # NB: The dither size is substracted twice to avoid edge issues

        if (self.obsmode == 'LIFU'):

            # Set the position of the guide camera with respect to the central
            # fibre of the LIFU bundle: the offset is described as movement
            # towards the North and East directions when the position angle is 0
            # (so negative values correspond to offsets towards to South and
            # West directions)

            self.cam_north_offset = -27.7 / 60.0
            self.cam_east_offset = 0.0 / 60.0

            cam_offset_dist = _np.hypot(self.cam_north_offset,
                                        self.cam_east_offset)

            # Set the effective radious of the field of view of the guide camera
            
            cam_width  =  4.0 / 60.0
            cam_height = 3.75 / 60.0
            
            self.fov_radious = (min(cam_width, cam_height) / 2 -
                                2 * self.dither_size_arcsec / 3600)

            # Set the limits for looking for stars
            
            self.ra_min = (self.ra - ((cam_offset_dist + self.fov_radious) /
                                      _np.cos(_np.deg2rad(self.dec))))
            self.ra_max = (self.ra + ((cam_offset_dist + self.fov_radious) /
                                      _np.cos(_np.deg2rad(self.dec))))
            self.dec_min = self.dec - (cam_offset_dist + self.fov_radious)
            self.dec_max = self.dec + (cam_offset_dist + self.fov_radious)

            import astropy.units as uu
            from astropy.coordinates import SkyCoord
            
        
            
        else:

            # Set the geometry of the guiding camera to None

            self.cam_north_offset = None
            self.cam_east_offset = None

            # Set the radious of the field of view

            self.fov_radious = 1.0 - 2 * self.dither_size_arcsec / 3600
            
            # Set the limits for looking for stars

            self.ra_min = self.ra - (self.fov_radious / 
                                     _np.cos(_np.deg2rad(self.dec)))
            self.ra_max = self.ra + (self.fov_radious / 
                                     _np.cos(_np.deg2rad(self.dec)))
            self.dec_min = self.dec - self.fov_radious
            self.dec_max = self.dec + self.fov_radious

        # Fix edge issues for declination (for right ascension is not needed)
        
        if self.dec_min < -90.0:
            self.dec_min = -90.0

        if self.dec_max > 90.0:
            self.dec_max = 90.0


    def _set_healpix_indices(self, num_sample=20):

        if self.healpix_indices is None:

            hp = _HEALPix(nside=self.nside, order='nested',
                          frame=_coordinates.ICRS())

            healpix_indices = []

            dec_min = self.dec_min
            dec_max = self.dec_max

            if dec_min < -90.0:
                dec_min = -90.0
            if dec_max > 90.0:
                dec_max = 90.0

            for ra in _np.linspace(self.ra_min, self.ra_max, num_sample):
                for dec in _np.linspace(dec_min, dec_max, num_sample):

                    healpix_index = hp.skycoord_to_healpix(
                        _coordinates.SkyCoord(ra, dec, unit='deg'))

                    if healpix_index not in healpix_indices:
                        healpix_indices.append(healpix_index)

            healpix_indices.sort()

        self.healpix_indices = healpix_indices


    def _get_cat(self, healpix_index, directory):

        raise NotImplementedError

        

    def _retrieve_cats(self,cache=False,prefilter=True):

        assert self.healpix_indices is not None

        table_list = []

        run_fn = DefaultDir
        if not cache:
            run_fn = _tempfile.TemporaryDirectory
        
        for healpix_index in self.healpix_indices:
            with run_fn() as dest_dir:
                file_path = self._get_cat(healpix_index, dest_dir)

                with _fits.open(file_path) as hdu_list:
                    aux_table = _Table(hdu_list[1].data)

            table_list.append(aux_table)

        full_table = _vstack(table_list)
        if (self.obsmode == 'LIFU') and (prefilter):
            weave_fov = 1.0
            # make early cuts here to slim down the potential number of sources
            centre_coord = _coordinates.SkyCoord(self.ra, self.dec, unit='deg')
            coords = _coordinates.SkyCoord(full_table['GAIA_RA'],full_table['GAIA_DEC'],unit='deg')
            dists = centre_coord.separation(coords).deg

            cam_offset_dist = _np.hypot(self.cam_north_offset, self.cam_east_offset)
            min_rad = (cam_offset_dist - self.fov_radious) * 0.95
            max_rad = (cam_offset_dist + self.fov_radious) * 1.05
            
            # mask = ((table['DISTANCE'] >= min_radious) *
            #         (table['DISTANCE'] <= max_radious))

            
            mask = (dists >= min_rad) * (dists <= max_rad)
            full_table = full_table[mask]
            logging.info(
                'LIFU pre-filtering of targets from {} candidates to {}'.format(
                    len(dists),len(full_table)))
        
        return full_table


    def _add_distance_and_angle_to_table(self, full_table):

        num_rows = len(full_table)

        centre_coord = _coordinates.SkyCoord(self.ra, self.dec, unit='deg')

        distance_list = []
        angle_list = []

        full_table_coords = _coordinates.SkyCoord(full_table['GAIA_RA'],
                                                  full_table['GAIA_DEC'],
                                                  unit='deg')

        distance_list = centre_coord.separation(full_table_coords).deg
        angle_list = centre_coord.position_angle(full_table_coords).wrap_at(
            _coordinates.Angle(180, unit='deg')).deg

        # for row in full_table:

        #     row_coord = _coordinates.SkyCoord(row['GAIA_RA'], row['GAIA_DEC'],
        #                                       unit='deg')

        #     distance_list.append(centre_coord.separation(row_coord).deg)
        #     angle_list.append(centre_coord.position_angle(row_coord).wrap_at(
        #         _coordinates.Angle(180, unit='deg')).deg)

        distance_column = _Column(distance_list, name='DISTANCE')
        angle_column = _Column(angle_list, name='ANGLE')

        full_table.add_column(distance_column)
        full_table.add_column(angle_column)

        return full_table


    def _set_full_table(self,cache=False):

        if self.full_table is None:
            self._set_healpix_indices()
            full_table = self._retrieve_cats(cache=cache)
            full_table = self._add_distance_and_angle_to_table(full_table)

            self.full_table = full_table


    def _get_stars_in_plate_a_b(self, table):

        mask = (table['DISTANCE'] <= self.fov_radious)

        fov_table = table[mask]

        return fov_table

    
    def _get_stars_in_guide_cam_ring(self, table):

        cam_offset_dist = _np.hypot(self.cam_north_offset, self.cam_east_offset)

        min_radious = cam_offset_dist - self.fov_radious
        max_radious = cam_offset_dist + self.fov_radious

        mask = ((table['DISTANCE'] >= min_radious) *
                (table['DISTANCE'] <= max_radious))

        ring_table = table[mask]

        return ring_table


    def _set_useful_table(self, mag_lim=99.0):
        if self.obsmode != 'LIFU':

            useful_table = self._get_stars_in_plate_a_b(self.full_table)

        else:
            
            useful_table = self._get_stars_in_guide_cam_ring(self.full_table)

        self.useful_table = useful_table
        

    def _select_targprog(self, table, targprogs):
        mask = _np.array([_targprog in targprogs for _targprog in table['TARGPROG']])
        selected_table = table[mask]
        non_selected_table = table[_np.logical_not(mask)]

        return selected_table, non_selected_table
        
        
    def _select_central_stars_in_plate_a_b(self, table, num_central_stars=1,
                                           mag_lim=99.):

        # Get the central stars if there was no mag limit
        argsort = _np.argsort(table['DISTANCE'])
        selected_indices = argsort[0:num_central_stars]
        mask = _np.array([(i in selected_indices)
                          for i in range(len(table))])

        # Identify all sources fainter than mag_lim (no matter where)
        all_maglim_rejected = table[_np.logical_not(table['GAIA_MAG_G'] <= mag_lim)]
        central_sel_no_maglim = table[mask]
        central_maglim_rejected = central_sel_no_maglim[central_sel_no_maglim['GAIA_MAG_G'] > mag_lim]
        

        # Now do for a mag-filtered subset table
        mag_mask = (table['GAIA_MAG_G'] <= mag_lim)
        _table = table[mag_mask]
        _argsort = _np.argsort(_table['DISTANCE'])
        _selected_indices = _argsort[0:num_central_stars]
        _mask = _np.array([(i in _selected_indices)
                          for i in range(len(_table))])
        all_non_central_maglim = _table[_np.logical_not(_mask)]

        central_star_table = _table[_mask]
        non_selected_table = _vstack([all_maglim_rejected,all_non_central_maglim])

        if mag_lim < 99:
            # add in the centrals that are in 0:num_central_stars, but filtered by mag
            # self.useful_maglim_table = central_maglim_rejected 
            
            
            if self.useful_maglim_table:
                self.useful_maglim_table = _vstack([self.useful_maglim_table, central_sel_no_maglim])
            else:
                # self.useful_maglim_table = central_sel
                self.useful_maglim_table = central_sel_no_maglim
                                                   
        # central_star_table = _table[_mask]
        # non_selected_table = table[_np.logical_not(mask)]

        return central_star_table, non_selected_table


    def _get_index_of_nearest_angle(self, ref_angle, angle_array,
                                    selected_indices=None,
                                    ranklist=False):

        index = None

        wrapped_angle_list = [
            _coordinates.Angle(angle, unit='deg').wrap_at(
                _coordinates.Angle(ref_angle + 180, unit='deg')).deg
            for angle in angle_array]

        ranking = {}
        
        for i, angle_i in enumerate(wrapped_angle_list):

            if (selected_indices == None) or (i not in selected_indices):

                angle_diff = _np.abs(angle_i - ref_angle)
                ranking[angle_diff] = i

                if (index == None) or (angle_diff < min_angle_diff):
                    index = i
                    min_angle_diff = angle_diff

        if ranklist:
            rank = list(ranking)
            rank.sort()
            rank = _np.array([ranking[r] for r in rank])
            return rank
        return index

    
    def _select_stars_in_ring_of_plate_a_b(self, table, num_stars,
                                           min_cut=0.9, max_cut=1.0,
                                           mag_lim=99.,
                                           rescap=True):

        min_radious = min_cut * self.fov_radious
        max_radious = max_cut * self.fov_radious

        mask = ((table['DISTANCE'] >= min_radious) *
                (table['DISTANCE'] <= max_radious))

        mag_maglim_reject_mask = mask * table['GAIA_MAG_G'] > mag_lim    
        ring_table_maglim_rejected = table[mag_maglim_reject_mask]
        mag_mask = table['GAIA_MAG_G'] <= mag_lim
        mask *= mag_mask
        ring_table = table[mask]

        if mag_lim < 99:
            if self.useful_maglim_table:
                # self.useful_maglim_table = _vstack([self.useful_maglim_table, table[mag_mask]])
                self.useful_maglim_table = _vstack([self.useful_maglim_table, ring_table_maglim_rejected])
            else:
                self.useful_maglim_table = table[mag_mask]

        if num_stars is None:

            selected_table = ring_table

        elif num_stars > 0:
            if len(ring_table) < num_stars:
                # There are fewer candidate stars than the number requested by the user - just include all the stars, no need to subsample
                selected_table = ring_table
            
            else:
                selected_indices = []

                # what is the resolution to subsample at?
                max_samples = num_stars
                resstr = ''
                if (num_stars > 360) and (rescap):
                    max_samples = 360
                    resstr = ' (capped due to rescap=True: requested num_stars will NOT be returned)'
                elif (num_stars > 1.8*360):
                    resstr = ' (high resolution: will take time to complete)'
                sampling_angles = _np.linspace(0, 360, max_samples, endpoint=False)
                dtheta = sampling_angles[1]            

                    
                logging.warning("Candidate stars in annulus ({}) exceeds requested num_stars ({}). Will azimuthally subsample at {:1.2f} deg resolution{}".format(
                    len(ring_table),num_stars,dtheta,resstr))
                    

                for angle in sampling_angles:
                    index = self._get_index_of_nearest_angle(
                        angle, ring_table['ANGLE'],
                        selected_indices=selected_indices)

                    if index is not None:
                        selected_indices.append(index)

                mask = _np.array([(i in selected_indices)
                                  for i in range(len(ring_table))])

                selected_table = ring_table[mask]

        else:

            selected_table = ring_table[0:0]

        return selected_table


    def _get_guide_cam_coord(self, pa):

        cam_offset_dist = _coordinates.Angle(
            _np.hypot(self.cam_north_offset, self.cam_east_offset), unit='deg')
        cam_pa_angle = _coordinates.Angle(
            _np.rad2deg(_np.arctan2(self.cam_north_offset,
                                    -self.cam_east_offset)) +
            pa - 90, unit='deg')

        centre_coord = _coordinates.SkyCoord(self.ra, self.dec, unit='deg')
        cam_coord = centre_coord.directional_offset_by(cam_pa_angle,
                                                       cam_offset_dist)

        return cam_coord


    def _select_stars_in_guide_cam(self, table, pa_request, num_stars_request,mag_lim=99.):

        # Get the coordinate of the guide cam for the requested PA

        cam_coord = self._get_guide_cam_coord(pa_request)

        # Compute the distances to the centre of the guide cam

        cam_distance_list = []

        for row in table:

            row_coord = _coordinates.SkyCoord(row['GAIA_RA'], row['GAIA_DEC'],
                                              unit='deg')

            cam_distance_list.append(cam_coord.separation(row_coord).deg)

        cam_distance_array = _np.array(cam_distance_list)

        # Get a table with the guide stars in the field of view of the guide
        # camera

        
        mask = (cam_distance_array <= self.fov_radious)

        pre_selected_table = table[mask]

            # filter out any from the selection that are not in the mag limit
            selected_maglim_rejected = pre_selected_table[_np.logical_not(pre_selected_table['GAIA_MAG_G'] <= mag_lim)]

            if self.useful_maglim_table:
                self.useful_maglim_table = _vstack([self.useful_maglim_table, selected_maglim_rejected])
            else:
                self.useful_maglim_table = selected_maglim_rejected
        
            mag_mask = table['GAIA_MAG_G'] <= mag_lim
            mask *= mag_mask
            pre_selected_table = table[mask]


        # Select the nearest star to the centre of the guide camera if it is in
        # the field of view. Otherwise, return an empty table

        if _np.count_nonzero(mask) > 0:

            argsort = _np.argsort(cam_distance_array[mask])

            if num_stars_request is not None:
                num_stars = num_stars_request
            else:
                num_stars = len(pre_selected_table)

            selected_table = pre_selected_table[argsort][0:num_stars]

        else:

            selected_table = pre_selected_table[0:0]

        return selected_table


    def _select_stars_in_near_guide_cam_position(self, table, pa_request,
                                                 num_stars_request,
                                                 mag_lim=99., ha_range=2.0,
                                                 ranklist=True, cache=False):

        # ranklist: provide a ranked list of closest, and go through this list until
        # a selection is within limits

               
        if len(table) > 0:
            obslim_py = get_obslimits_tool(cache=cache)
            sys.path = [os.path.dirname(obslim_py)] + sys.path
            from obs_limits import ObsLimits

            xml_configure = get_configure_cfg(cache=cache) 
            pa_check = ObsLimits(cfg=xml_configure)                    

            cam_pa_at_request = _coordinates.Angle(
                _np.rad2deg(_np.arctan2(self.cam_north_offset,
                                        -self.cam_east_offset)) +
                pa_request - 90, unit='deg').wrap_at(
                    _coordinates.Angle(180, unit='deg')).deg

            argmins = self._get_index_of_nearest_angle(cam_pa_at_request,
                                                       table['ANGLE'],ranklist=True)
            selected_pa = None
            ranking = 0
            for _argmin in argmins:
                ranking += 1
                cam_pa_at_zero = _coordinates.Angle(
                    _np.rad2deg(_np.arctan2(self.cam_north_offset,
                                            -self.cam_east_offset)) -
                    90, unit='deg').wrap_at(
                        _coordinates.Angle(180, unit='deg')).deg

                this_pa = _coordinates.Angle(
                    table['ANGLE'][_argmin] - cam_pa_at_zero, unit='deg').wrap_at(
                        _coordinates.Angle(180, unit='deg')).deg

                limit_check = pa_check.limit_check(ha_range*15., self.dec, this_pa, detail=True, verbose=False)
                pa_valid = limit_check[0]
                checks = limit_check[1]
                if pa_valid:
                    selected_pa = this_pa
                    logging.info("Closest viable alternative PA ({:.2f}) was ranked {}/{}".format(selected_pa,ranking,len(argmins)))
                    break
            assert selected_pa, "There were no viable guide stars that permit observation within rotator limits for this Declination and Hour Angle range"

            selected_table = self._select_stars_in_guide_cam(table, selected_pa,
                                                             num_stars_request,
                                                             mag_lim=mag_lim)

        else:

            selected_pa = pa_request
            selected_table = table
            
        return selected_pa, selected_table


    def _set_selected_table(self, num_stars_request, pa_request=None,
                            num_central_stars=1, min_cut=0.9, max_cut=1.0,
                            targprog=None,mag_lim=99.,cache=False,ha_range=2.0):

        useful_table = self.useful_table
        
        
        if targprog:
            # filter according to value required in targprog
            useful_table, rejected = self._select_targprog(useful_table,targprog)

        # If there is not request, let's start with a position angle of 0.0

        if pa_request not in [None, _np.nan]:
            pa_request = pa_request
        else:
            pa_request = 0.0

        if self.obsmode != 'LIFU':

            # For MOS/mIFU, get central stars and other ones in the outer ring
            # of the plate
            
            selected_pa = pa_request

            central_star_table, non_selected_table = (
                self._select_central_stars_in_plate_a_b(
                    useful_table, num_central_stars=num_central_stars,mag_lim=mag_lim))

            if num_stars_request is not None:
                num_ring_stars = num_stars_request - num_central_stars
            else:
                num_ring_stars = None

            ring_stars_table = self._select_stars_in_ring_of_plate_a_b(
                non_selected_table, num_ring_stars,
                min_cut=min_cut, max_cut=max_cut, mag_lim=mag_lim)

            selected_table = _vstack([central_star_table, ring_stars_table])

        else:

            # For LIFU, try to get a star in the central in the guide cam, and
            # if not possible, choose the nearest available star to the desired
            # position

            selected_table = self._select_stars_in_guide_cam(useful_table,
                                                             pa_request,
                                                             num_stars_request,
                                                             mag_lim=mag_lim)
            obslim_py = None
            if len(selected_table) > 0:
                # There are guide stars at the PA requested - good!
                selected_pa = pa_request
            else:
                # No guides in guidecam FOV - go hunting for another PA nearby
                selected_pa, selected_table = (
                    self._select_stars_in_near_guide_cam_position(
                        useful_table, pa_request, num_stars_request,
                        mag_lim=mag_lim, ha_range=ha_range, cache=cache))

                # Now need to do a check on the new PA, to ensure it is within rotator
                # limits based on the Declination and a canonical HA range
                obslim_py = get_obslimits_tool(cache=cache)
                sys.path = [os.path.dirname(obslim_py)] + sys.path
                # sys.path.append(os.path.dirname(obslim_py))
                from obs_limits import ObsLimits
                
                xml_configure = get_configure_cfg(cache=cache) 
                pa_check = ObsLimits(cfg=xml_configure)                    
                limit_check = pa_check.limit_check(ha_range*15., self.dec, selected_pa, detail=True, verbose=False)
                pa_valid = limit_check[0]
                checks = limit_check[1]
                if not pa_valid:
                    pa_ranges = pa_check.pa_range(ha_range*15.,self.dec,verbose=False,explicit=True)
                    range_plural = ''
                    viable_ranges = ""
                    for _range in pa_ranges:
                        viable_ranges += " [{:.1f},{:.1f}]".format(_range[0],_range[1])
                        if _range != pa_ranges[0]:
                            range_plural = 's'

                    
                    # pamin, pamax = pa_check.pa_range(ha_range,_dec,verbose=False)
                    logging.error("There were no guide stars available at the requested PA ({:.1f}). The closest viable PA selected by the workflow ({:.1f}) exceeds the rotator limits, and is therefore unobservable".format(pa_request,selected_pa))
                    logging.error("This arises because your *original* pointing at Declination {:.3f} with {}hr HA angle range has changed PA".format(self.dec,ha_range))
                    if obslim_py:
                        logging.error("Please use the obs_limits.py tool (local copy: {}) to check the legal PA for your pointing".format(obslim_py))
                    else:
                        logging.error("Please use the obs_limits.py tool to check the legal PA for your pointing")
                    logging.error("We advise users to identify the viable PA range and select a PA at the mid point of this.")
                    logging.error("Viable PA range{}:{}".format(range_plural,viable_ranges))
                    logging.error("This XML should be re-generated from stage1, re-specifying the PA to ensure the rotator limits are not exceeded")

                    raise SystemExit(1)

        self.selected_table = selected_table
        self.selected_pa = float(selected_pa)

        
    def _print_summary(self,mag_lim):

        if self.full_table is not None:
            logging.info(
                '{} {} stars have been retrieved for ({}, {})'.format(
                    len(self.full_table), self.star_type, self.ra, self.dec))

        if self.selected_table is not None:
            logging.info(
                '{} {} stars are useful for ({}, {})'.format(
                    len(self.useful_table), self.star_type, self.ra, self.dec))

        if self.selected_table is not None:
            if self.obsmode != 'LIFU':
                logging.info(
                    '{} {} stars has been selected for ({}, {})'.format(
                        len(self.selected_table), self.star_type,
                        self.ra, self.dec))
            else:
                logging.info(
                    '{} {} stars has been selected for ({}, {})'.format(
                        len(self.selected_table), self.star_type,
                        self.ra, self.dec) +
                    ' PA = {} deg'.format(self.selected_pa))

            if (not len(set(self.selected_table['GAIA_DR']))) and (mag_lim < 99):
                # a magnitude limit was used for this selection:
                # issue a warning it might have been too strict?
                logging.warning(
                    "A faint {} magnitude cut of {:.2f} was used for a request that returned no targets: consider relaxing the constraint".format(
                        self.star_type,mag_lim))
            if not len(self.selected_table):
                logging.warning("No {} targets in selection".format(self.star_type))
                return

            assert len(set(self.selected_table['GAIA_DR'])) == 1
            
            gaia_dr_value = str(self.selected_table['GAIA_DR'][0])
            
            logging.info(
                'Notice that GAIA_DR value for the selected stars is {}'.format(
                    gaia_dr_value))

    def _fix_ra_for_plotting(self, ra_array):

        mask_below_ra_min = (ra_array < self.ra_min)

        if _np.count_nonzero(mask_below_ra_min) > 0:
            ra_array[mask_below_ra_min] = ra_array[mask_below_ra_min] + 360

        mask_above_ra_max = (ra_array > self.ra_max)

        if _np.count_nonzero(mask_above_ra_max) > 0:
            ra_array[mask_above_ra_max] = ra_array[mask_above_ra_max] - 360

        return ra_array

                
    def _plot_table(self, table, ax, color, ringcolor=None, label=None):

        if (table is not None) and (len(table) > 0):

            # Get the coordinates from the table

            ra_array = table['GAIA_RA']
            dec_array = table['GAIA_DEC']

            # Fix coordinates where there are edge problems

            ra_array = self._fix_ra_for_plotting(ra_array)

            # Plot the coordinates

            ax.scatter(ra_array, dec_array, c=color, edgecolors=ringcolor, label=label)


    def _plot_useful_region_in_plate(self, ax):
    
        center = (self.ra, self.dec)

        width = self.ra_max - self.ra_min
        height = self.dec_max - self.dec_min

        ellipse = _Ellipse(center, width, height,
                           facecolor='lightgrey', edgecolor='none',
                           linewidth=3, zorder=-10)

        ax.add_patch(ellipse)

                
    def _plot_cut_region_in_plate(self, ax, cut):

        if cut is not None:

            center = (self.ra, self.dec)

            width = (self.ra_max - self.ra_min) * cut
            height = (self.dec_max - self.dec_min) * cut

            ellipse = _Ellipse(center, width, height,
                            facecolor='none', edgecolor='k',
                            linewidth=1, linestyle='--', zorder=-8)

            ax.add_patch(ellipse)

                
    def _plot_useful_region_for_guide_cam(self, ax):
    
        center = (self.ra, self.dec)

        out_width = self.ra_max - self.ra_min
        out_height = self.dec_max - self.dec_min

        out_ellipse = _Ellipse(center, out_width, out_height,
                               facecolor='lightgrey', edgecolor='none',
                               linewidth=0, zorder=-10)

        ax.add_patch(out_ellipse)

        in_width = ((self.ra_max - self.ra_min) -
                    2 * 2 * self.fov_radious / _np.cos(_np.deg2rad(self.dec)))
        in_height = ((self.dec_max - self.dec_min) -
                     2 * 2 * self.fov_radious)

        _ellipse = _Ellipse(center, in_width, in_height,
                               facecolor='w', edgecolor='none',
                               linewidth=0, zorder=-9)

        ax.add_patch(_ellipse)


    def _plot_guide_cam_fov(self, ax, pa, linestyle='-', color='k', zorder=-5):

        if pa not in [None, _np.nan]:

            guide_coord = self._get_guide_cam_coord(pa)

            ra = guide_coord.ra.deg
            dec = guide_coord.dec.deg

            if ra < self.ra_min:
                ra += 360

            if ra > self.ra_max:
                ra -= 360

            center = (ra, dec)

            width = 2 * self.fov_radious / _np.cos(_np.deg2rad(self.dec))
            height = 2 * self.fov_radious

            ellipse = _Ellipse(center, width, height,
                               facecolor='none', edgecolor=color,
                               linewidth=1, linestyle=linestyle, zorder=zorder)

            ax.add_patch(ellipse)

            ax.scatter([ra], [dec], c=color, marker='x', zorder=zorder)
    
                
    def _plot(self, plot_filename, pa_request='',
              num_stars='', num_stars_request='',
              num_central_stars='', min_cut='', max_cut='', mag_lim=99.):

        fig, ax = _plt.subplots(
            subplot_kw={'position': [0.03, 0.13, 0.75, 0.75]})

        ax.set_title(
            '{}\n{} stars at ({:.5f}, {:.5f})'.format(
                self.obsmode, self.star_type.capitalize(), self.ra, self.dec))

        ax.set_xlabel('RA (deg)')
        ax.set_ylabel('Dec (deg)')

        text_str = ('num_stars:\n' +
                    '   {}\n'.format(num_stars) +
                    'num_stars_request:\n' +
                    '   {}\n\n'.format(num_stars_request))

        if self.obsmode != 'LIFU':
            self._plot_useful_region_in_plate(ax)
            self._plot_cut_region_in_plate(ax, min_cut)
            self._plot_cut_region_in_plate(ax, max_cut)

            text_str += ('num_central_stars:\n' +
                         '   {}\n'.format(num_central_stars) +
                         'min_cut:\n' +
                         '   {}\n'.format(min_cut) +
                         'max_cut:\n' +
                         '   {}\n'.format(max_cut))                         

        else:
            self._plot_useful_region_for_guide_cam(ax)
            self._plot_guide_cam_fov(ax, pa_request,
                                     linestyle=':', color='blue')
            self._plot_guide_cam_fov(ax, self.selected_pa,
                                     linestyle='--', color='k')

            if type(self.selected_pa) == float:
                selected_pa_str = '{:.3f}'.format(self.selected_pa)
            else:
                selected_pa_str = '{}'.format(self.selected_pa)

            if type(pa_request) == float:
                pa_request_str = '{:.3f}'.format(pa_request)
            else:
                pa_request_str = '{}'.format(pa_request)

            text_str += ('selected_pa:\n' +
                         '   {} deg\n'.format(selected_pa_str) +
                         'pa_request:\n' +
                         '   {} deg\n'.format(pa_request_str))

        if mag_lim < 99:
            text_str += ('mag_cut:\n' +
                         '   {:.2f}\n'.format(mag_lim))

                         
        ax.text(1.05, 1.0, text_str, transform=ax.transAxes, fontsize=12,
                va='top', ha='left')

        self._plot_table(self.full_table, ax, 'red', label="All")
        self._plot_table(self.useful_table, ax, 'darkorange', label="Useful")
        if (self.useful_maglim_table) and (len(self.useful_maglim_table) != len(self.useful_table)):
            self._plot_table(self.useful_maglim_table, ax, 'orange', ringcolor='green', label="Mag cut") 

        self._plot_table(self.selected_table, ax, 'green', label="Selected")

        ax.scatter([self.ra], [self.dec], c='k', marker='+')

        ax.legend(bbox_to_anchor=(1.04, 0), loc="lower left", borderaxespad=0)
        ax.set_xlim([self.ra_max, self.ra_min])
        ax.set_ylim([self.dec_min, self.dec_max])
        ax.set_aspect(1. / _np.cos(_np.deg2rad(self.dec)))

        fig.savefig(plot_filename)

                
    def get_table(self, selection='filter', verbose=True, plot_filename=None,
                  useful_table_filename=None, pa_request=None,
                  num_stars_request=None, num_central_stars=1, min_cut=0.9,
                  max_cut=1.0,cache=False, mag_lim=99., ha_range=2.0,
                  calib_mixres=False):

        # Assert the input values

        assert selection in ['filter', 'useful', 'full']
        
        if (num_stars_request) is not None:
            assert num_stars_request > 0
        
        if pa_request not in [None, _np.nan]:
            assert (pa_request > -180.0) or (pa_request < 180.0)

        # Set the needed tables
        self._set_full_table(cache=cache)

        if selection in ['useful', 'filter']:

            self._set_useful_table(mag_lim=mag_lim)

            if selection == 'filter':
                targprog_lookup = {}
                metric = {}
                metric['WG'] = self.obsmode
                metric['WD'] = self.resolution
                targsrvy = self.useful_table[0]['TARGSRVY']
                targprog_lookup['WD'] = {"HR":["HR"],"LR":["LR"]}
                targprog_lookup['WG'] = {"mIFU":["ANY"],"MOS":["ANY"],"LIFU":["LIFU","ANY"]}
                    
                targprog = targprog_lookup[targsrvy][metric[targsrvy]]

                self._set_selected_table(num_stars_request,
                                         pa_request=pa_request,
                                         num_central_stars=num_central_stars,
                                         min_cut=min_cut, max_cut=max_cut,
                                         targprog=targprog, mag_lim=mag_lim,
                                         cache=cache, ha_range=ha_range) 

                if (calib_mixres) and (targsrvy == 'WD'):
                    # switch the retrieval to LR only:
                    targprog_lookup['WD']['HR'] = ['LR']

                    # copy the existing HR table, then reset
                    selected_tableHR = deepcopy(self.selected_table)
                    self.selected_table = None
                    
                    targprog = targprog_lookup[targsrvy][metric[targsrvy]]
                    _num_stars_request = num_stars_request
                    if _num_stars_request != None:
                        _num_stars_request = _num_stars_request - len(selected_tableHR)
                    if (_num_stars_request == None) or (_num_stars_request > 0):
                        self._set_selected_table(_num_stars_request,
                                                 pa_request=pa_request,
                                                 num_central_stars=num_central_stars,
                                                 min_cut=min_cut, max_cut=max_cut,
                                                 targprog=targprog, mag_lim=mag_lim,
                                                 cache=cache, ha_range=ha_range) 
                        _selected_table = _vstack([selected_tableHR,self.selected_table])
                        logging.info("Using LR ({}) and HR ({}) datasets, a total of {} targets are available in selection".format(
                            len(selected_tableHR),len(self.selected_table),len(_selected_table)))
                        self.selected_table = _selected_table

                    elif (_num_stars_request != None) and (_num_stars_request <= 0):
                        logging.warning("calib_mixres has been set, but the number of HR calibs is at your requested num_stars_request limit ({}).".format(num_stars_request))
                        logging.warning("If you want to also include LR targets, re-run this function without a num_stars_request value.")

        # Choose the output depending on the selection method

        if selection == 'full':
            pa = pa_request
            table = self.full_table
        elif selection == 'useful':
            pa = pa_request
            table = self.useful_table
        elif selection == 'filter':
            pa = self.selected_pa
            table = self.selected_table

        # Print a summary and plot a figure if requested

        if verbose == True:
            self._print_summary(mag_lim)

        if plot_filename is not None:

            num_stars = len(table)

            self._plot(plot_filename, pa_request=pa_request,
                       num_stars=num_stars, num_stars_request=num_stars_request,
                       num_central_stars=num_central_stars,
                       min_cut=min_cut, max_cut=max_cut, mag_lim=mag_lim)
        
        if ((useful_table_filename is not None) and
            (selection in ['useful', 'filter'])):
            
            self.useful_table.write(useful_table_filename)

        return pa, table

    
class GuideStars(_AuxStars):
    """
    Handle the retrieval of WEAVE guide stars.
    
    This class provides a mechanism for querying and retrieving guide star
    targets for the construction of WEAVE 'protofields'.
    
    Parameters
    ----------
    ra : float
        The right ascension (decimal degrees) of either the central spaxel or
        central FoV.
    dec : float
        The declination (decimal degrees) of either the central spaxel or
        central FoV.
    obsmode : str
        Either LIFU, mIFU or MOS.
    resolution : str
        Either LR or HR
    dither_size_arcsec : float, optional
        Maximum allowed size for dithering patterns.
    nside : int, optional
        Override the default HEALPix nside value. Will likely end in tears.
    """

    def __init__(self, ra, dec, obsmode, resolution, dither_size_arcsec=10.0, nside=32):

        super().__init__(ra, dec, obsmode, resolution,
                         dither_size_arcsec=dither_size_arcsec, nside=nside)

        self.star_type = 'guide'


    def _get_cat(self, healpix_index, directory):
        file_path = _get_guide_cat(healpix_index, directory)

        return file_path


    def get_table(self, selection='filter', verbose=True, plot_filename=None,
                  useful_table_filename=None, pa_request=None,
                  num_stars_request=None, num_central_stars=1, min_cut=0.0,
                  max_cut=0.9,cache=False, mag_lim=99., ha_range=2.0):
        """
        Get a table with a set of selected stars with the requested conditions.

        Parameters
        ----------
        selection : {'filter', 'useful', 'full'}
            Selection to be made in the output table: 'full' returns all the
            stars near to the coordinates, 'useful' returns all the stars which
            could be used for guiding at that coordinates (ignoring any request
            for the position angle), 'filter' returns a table following the
            prescriptions given in the other keywords of this method.
        verbose : bool
            Print a summary or not.
        plot_filename : str, optional
            Name of the file to save a figure with the full table, the useful
            table and the selected table.
        useful_table_filename : str, optional
            Name of the file to save the useful table.
        pa_request : float, optional
            The position angle (degrees) of rotation (it can be different to
            zero only for LIFU). np.nan and None are also valid values.
        num_stars_request : int, optional
            Maximum number of guide stars in the output. None means no limit.
        num_central_stars : int, optional
            Number of stars near to centre to be selected (only used when
            obsmode is not LIFU).
        min_cut : float, optional
            Minimum cut factor to be used for the non-central stars (only used
            when obsmode is not LIFU).
        max_cut : float, optional
            Maximum cut factor to be used for the non-central stars (only used
            when obsmode is not LIFU).
        cache : bool, optional
            Do not download catalogue files if they aready exist
        mag_lim : float, optional
            Apply a faint limit cut in Gaia G to the available stars

        Returns
        -------
        pa : float
            The selected position angle.
        table : astropy.table.Table
            A table containing the guide stars.
        """

        pa, table = super().get_table(selection=selection, verbose=verbose,
                                      plot_filename=plot_filename,
                                      useful_table_filename=useful_table_filename,
                                      pa_request=pa_request,
                                      num_stars_request=num_stars_request,
                                      num_central_stars=num_central_stars,
                                      min_cut=min_cut, max_cut=max_cut,
                                      cache=cache, mag_lim=mag_lim, ha_range=ha_range)

        return pa, table


class CalibStars(_AuxStars):
    """
    Handle the retrieval of WEAVE calib stars.
    
    This class provides a mechanism for querying and retrieving calib star
    targets for the construction of WEAVE 'protofields'.
    
    Parameters
    ----------
    ra : float
        The right ascension (decimal degrees) of either the central spaxel or
        central FoV.
    dec : float
        The declination (decimal degrees) of either the central spaxel or
        central FoV.
    obsmode : str
        Either LIFU, mIFU or MOS.
    resolution : str
        Either LR or HR
    dither_size_arcsec : float, optional
        Maximum allowed size for dithering patterns.
    nside : int, optional
        Override the default HEALPix nside value. Will likely end in tears.
    """


    def __init__(self, ra, dec, obsmode, resolution, dither_size_arcsec=10.0, nside=32, res=None):

        assert obsmode != 'LIFU'

        super().__init__(ra, dec, obsmode, resolution,
                         dither_size_arcsec=dither_size_arcsec, nside=nside)

        self.star_type = 'calib'


    def _get_cat(self, healpix_index, directory):
        file_path = _get_calib_cat(healpix_index, directory)

        return file_path


    def get_table(self, selection='filter', verbose=True, plot_filename=None,
                  useful_table_filename=None, num_stars_request=None,
                  num_central_stars=0, min_cut=0.0, max_cut=0.9, cache=False,
                  mag_lim=99.,calib_mixres=False):
        """
        Get a table with a set of selected stars with the requested conditions.

        Parameters
        ----------
        selection : {'filter', 'useful', 'full'}
            Selection to be made in the output table: 'full' returns all the
            stars near to the coordinates, 'useful' returns all the stars which
            which are in the field of view at that coordinates, 'filter' returns
            a table following the prescriptions given in the other keywords of
            this method.
        verbose : bool
            Print a summary or not.
        plot_filename : str, optional
            Name of the file to save a figure with the full table, the useful
            table and the selected table.
        useful_table_filename : str, optional
            Name of the file to save the useful table.
        num_stars_request : int, optional
            Maximum number of guide stars in the output. None means no limit.
        num_central_stars : int, optional
            Number of stars near to centre to be selected.
        min_cut : float, optional
            Minimum cut factor to be used for the non-central stars.
        max_cut : float, optional
            Maximum cut factor to be used for the non-central stars.
        cache : bool, optional
            Do not download catalogue files if they aready exist
        mag_lim : float, optional
            Apply a faint limit cut in Gaia G to the available stars
        calib_mixres : bool, optional
            Permit selection of LR (fainter) calibration targets under a HR field
            where no (or fewer) HR calibration targets are available.

        Returns
        -------
        table : astropy.table.Table
            A table containing the calib stars.
        """
        pa, table = super().get_table(selection=selection, verbose=verbose,
                                      plot_filename=plot_filename,
                                      useful_table_filename=useful_table_filename,
                                      pa_request=None,
                                      num_stars_request=num_stars_request,
                                      num_central_stars=num_central_stars,
                                      min_cut=min_cut, max_cut=max_cut,
                                      cache=cache, mag_lim=mag_lim,
                                      calib_mixres=calib_mixres)

        return table


if __name__ == '__main__':

    # Get a parser

    import argparse

    parser = argparse.ArgumentParser(
        description='Get WEAVE guide and calib stars')

    parser.add_argument('type_stars', choices=['guide', 'calib'],
                        help='type of stars to be searched')
    parser.add_argument('obsmode', choices=['MOS', 'LIFU', 'mIFU'],
                        help='obsmode to be considered in the search')
    parser.add_argument('central_ra', type=float,
                        help='RA in degrees of the centre of the FoV')
    parser.add_argument('central_dec', type=float,
                        help='Dec in degrees of the centre of the FoV')
    parser.add_argument('--res', choices=['LR','HR'],
                        help='The resolution of the observation (calib mode only)')
    parser.add_argument('--pa_request', default=-1,
                        help="""Requested PA in degress for LIFU observations;
                        -1 means no limit""")
    parser.add_argument('--lifu_num_guide_stars_request', default=1,
                        help="""maximum number of LIFU guide stars in the
                        output; -1 means no limit""")
    parser.add_argument('--mifu_num_guide_stars_request', default=-1,
                        help="""maximum number of mIFU guide stars in the
                        output; -1 means no limit""")
    parser.add_argument('--mifu_num_central_guide_stars', default=1, type=int,
                        help="""number of mIFU guide stars near to centre to be
                        selected""")
    parser.add_argument('--mifu_min_guide_cut', default=0.3, type=float,
                        help="""minimum cut factor to be used for the
                        non-central mIFU guide stars""")
    parser.add_argument('--mifu_max_guide_cut', default=0.9, type=float,
                        help="""maximum cut factor to be used for the
                        non-central mIFU guide stars""")
    parser.add_argument('--num_calib_stars_request', default=-1, type=int,
                        help="""maximum number of calib stars in the output;
                        -1 means no limit""")
    parser.add_argument('--num_central_calib_stars', default=0, type=int,
                        help="""number of calib stars near to centre to be
                        selected""")
    parser.add_argument('--min_calib_cut', default=0.0, type=float,
                        help="""minimum cut factor to be used for the
                        non-central calib stars""")
    parser.add_argument('--max_calib_cut', default=1.0, type=float,
                        help="""maximum cut factor to be used for the
                        non-central calib stars""")
    parser.add_argument('--cache', action='store_true',
                        help='Do not download catalogue files that already exist')
    parser.add_argument('--guide_maglim', default=99.0, type=float,
                        help="""Maximum Gaia G magnitude to be used for
                        the guide stars""")
    parser.add_argument('--calib_maglim', default=99.0, type=float,
                        help="""Maximum Gaia G magnitude to be used for
                        the calibration stars""")
    parser.add_argument('--calib_mixres', action='store_true',
                        help='Permit selection of LR WD in HR fields')
    parser.add_argument('--no_verbose', dest='verbose', action='store_false',
                        help='deactivate verbose output')

    # Get the arguments from the parser

    args = parser.parse_args()
    # In arguments which accept None value, do casting or transform None values
    
    if args.pa_request > -1:
        pa_request = float(args.pa_request)
    else:
        pa_request = None

    if args.lifu_num_guide_stars_request > -1:
        lifu_num_guide_stars_request = int(args.lifu_num_guide_stars_request)
    else:
        lifu_num_guide_stars_request = None
    
    if args.mifu_num_guide_stars_request > -1:
        mifu_num_guide_stars_request = int(args.mifu_num_guide_stars_request)
    else:
        mifu_num_guide_stars_request = None
    
    if args.num_calib_stars_request > -1:
        num_calib_stars_request = int(args.num_calib_stars_request)
    else:
        num_calib_stars_request = None

    # Set the logging level
    
    logging.basicConfig(level=logging.INFO)
    
    # Choose a filename for the output plot
    
    if (args.obsmode == 'LIFU') and (pa_request is not None):
    
        if pa_request is not None:
            pa_request_str = '{:.3f}'.format(pa_request)
        else:
            pa_request_str = 'none'
    
        plot_filename = '{}s_lifu_{:.5f}_{:.5f}_{}.png'.format(
            args.type_stars, args.central_ra, args.central_dec, pa_request_str)
    else:
        plot_filename = '{}s_{:.5f}_{:.5f}.png'.format(
            args.type_stars, args.central_ra, args.central_dec)

    # Get the requested type of stars

    if args.type_stars == 'guide':

        if args.obsmode == 'LIFU':
            num_stars_request = lifu_num_guide_stars_request
        else:
            num_stars_request = mifu_num_guide_stars_request
    
        guide_stars = GuideStars(args.central_ra, args.central_dec,
                                 args.obsmode)

        actual_pa, stars_table = guide_stars.get_table(
            verbose=args.verbose, plot_filename=plot_filename,
            pa_request=pa_request, num_stars_request=num_stars_request,
            num_central_stars=args.mifu_num_central_guide_stars,
            min_cut=args.mifu_min_guide_cut, max_cut=args.mifu_max_guide_cut,
            cache=args.cache,mag_lim=args.guide_maglim)

    elif args.type_stars == 'calib':
        
        assert args.obsmode != 'LIFU'
        res_index = next((i for i, action in enumerate(parser._actions) if action.dest == 'res'), None)
        opts = parser._actions[res_index].choices
        assert args.res != None, "You must explicitly define a resolution ({})".format(opts)
        
        calib_stars = CalibStars(args.central_ra, args.central_dec,
                                 args.obsmode, args.res)
        stars_table = calib_stars.get_table(
            verbose=args.verbose, plot_filename=plot_filename,
            num_stars_request=num_calib_stars_request,
            num_central_stars=args.num_central_calib_stars,
            min_cut=args.min_calib_cut, max_cut=args.max_calib_cut,
            cache=args.cache, mag_lim=args.calib_maglim,
            calib_mixres=args.calib_mixres)

    else:

        raise ValueError
    
    # Print the most interesting columns of the returned table
    
    stars_table.keep_columns(
        ['CNAME', 'TARGPROG', 'GAIA_RA', 'GAIA_DEC', 'DISTANCE', 'ANGLE'])
    
    for line in str(stars_table).split('\n'):
        logging.info(line)
    
    # Report the plot file created
    
    logging.info('Plot available at file {}'.format(plot_filename))

