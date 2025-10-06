import os
import pandas as pd
import numpy as np

# import matplotlib.pyplot as plt
# from matplotlib.patches import Circle

import logging
import time
import xml.dom.minidom

from astropy.table import Table, vstack
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u

from workflow.utils.classes import OBXML
from workflow.utils.classes._guide_and_calib_stars import GuideStars as _GuideStars
from workflow.utils.classes._guide_and_calib_stars import CalibStars as _CalibStars


from workflow.utils.get_resources import (get_blank_xml_template,
                                          get_progtemp_file,
                                          get_obstemp_file)
# from workflow.utils import get_version # Currently not used. will likely be re-used later
from workflow.utils.get_obstemp_info import (get_obstemp_info,
                                             get_obstemp_dict)

from workflow.utils.get_progtemp_info import (get_progtemp_dict,
                                              get_progtemp_info,
                                              get_resolution_from_progtemp)


class TLOBXML(OBXML):
    """ Modified OBXML Class for WEAVE MOS XML creation

    This class builds upon the IFU/MOS workflow OBXML class, adds new functionality
    necessary to allow for WEAVE-TL (Tumble-Less) OBs. Unless specifically
    overwritten by the new attributes below, all functionality is inherited from
    the base OBXML class.
    """

    def create_and_prepare_fields(self, field_ra, field_dec, central_target_table):

        field_template = self.fields.getElementsByTagName('field')[0].cloneNode(True)

        field_coord_table = central_target_table[['GAIA_ID', 'GAIA_RA', 'GAIA_DEC',
                                                  'TARGPROG', 'TARGID']]

        targprog = field_coord_table[(field_coord_table['GAIA_RA'] == field_ra) &
                                     (field_coord_table['GAIA_DEC'] == field_dec)
                                     ]['TARGPROG'].value[0]

        targid = field_coord_table[(field_coord_table['GAIA_RA'] == field_ra) &
                                   (field_coord_table['GAIA_DEC'] == field_dec)
                                   ]['TARGID'].value[0]

        for node in field_template.childNodes:
            if node.nodeType is xml.dom.minidom.Node.COMMENT_NODE:
                field_template.removeChild(node)

        # Remove existing targets from the template to start fresh
        for node in field_template.getElementsByTagName('target'):
            field_template.removeChild(node)

        # Set field attributes including the 'order'
        field_attrib_dict = {
            'Dec_d': field_dec,
            'RA_d': field_ra,
            'order': ""
        }
        self._set_attribs(field_template, field_attrib_dict)

        # Set targid atribute for configure tl
        field_template.setAttribute('targid', targid)
        field_template.setAttribute('targprog', targprog)
        # print('initial', field_template.toxml())
        return field_template, targprog

    def add_table_as_targets(self, field_element, table):

        """Add science targets to the XML file filling in target attributes"""

        col_to_attrib_target_dict = {
                                     'cname': 'CNAME',
                                     'targcat': 'TARGCAT',
                                     'targclass': 'TARGCLASS',
                                     'targdec': 'GAIA_DEC',
                                     'targepoch': 'GAIA_EPOCH',
                                     'targid': 'GAIA_ID',
                                     'targname': 'TARGNAME',
                                     'targparal': 'GAIA_PARAL',
                                     'targpmdec': 'GAIA_PMDEC',
                                     'targpmra': 'GAIA_PMRA',
                                     'targprio': 'TARGPRIO',
                                     'targprog': 'TARGPROG',
                                     'targra': 'GAIA_RA',
                                     'targsrvy': 'TARGSRVY',
                                     'targuse': 'TARGUSE'
                                    }

        col_to_attrib_photometry_dict = {
                                         'emag_bp': 'GAIA_MAG_BP_ERR',
                                         'emag_gg': 'GAIA_MAG_G_ERR',
                                         'emag_rp': 'GAIA_MAG_RP_ERR',
                                         'mag_bp': 'GAIA_MAG_BP',
                                         'mag_gg': 'GAIA_MAG_G',
                                         'mag_rp': 'GAIA_MAG_RP'
                                        }

        target_template = self.fields.getElementsByTagName('target')[0].cloneNode(True)

        # Check if fibreid is present in the target element remove if thats the case
        # This happens because we are cloning the first target element, which will be a
        # guide star target element
        if target_template.hasAttribute('fibreid'):
            target_template.removeAttribute('fibreid')

        for row in table:
            target = target_template.cloneNode(True)
            target_attrib_dict = {key: row[col] for key, col in
                                  col_to_attrib_target_dict.items()}
            # print(target_attrib_dict)
            self._set_attribs(target, target_attrib_dict)
            photometry = target.getElementsByTagName('photometry')[0]
            photometry_attrib_dict = {key: row[col] for key, col in
                                      col_to_attrib_photometry_dict.items()}
            self._set_attribs(photometry, photometry_attrib_dict)

            field_element.appendChild(target)

        # Append the fully prepared field element to the fields
        self.fields.appendChild(field_element)

    def set_tl_fields(self, field_ra, field_dec):

        field_list = self.fields.getElementsByTagName('field')
        field_template = field_list[0].cloneNode(True)

        for node in field_template.getElementsByTagName('target'):
            field_template.removeChild(node)

        for node in field_template.childNodes:
            if node.nodeType is xml.dom.minidom.Node.COMMENT_NODE:
                field_template.removeChild(node)

        # Get a target element and save it as template
        target_list = self.fields.getElementsByTagName('target')
        target_template = None

        for target in target_list:
            if target.getAttribute('targuse') == 'T':
                target_template = target.cloneNode(True)
                break

        assert target_template is not None

        # Clean the fields element
        for node in self.fields.getElementsByTagName('field'):
            self.fields.removeChild(node)
        order = ""  # self.first_science_order

        field = field_template.cloneNode(True)
        field_attrib_dict = {
            'Dec_d': field_dec,
            'RA_d': field_ra,
            'order': order
            }

        self._set_attribs(field, field_attrib_dict)
        self.central_ra = field_ra
        self.central_dec = field_dec

        # Add a targets to the field
        target = target_template.cloneNode(True)

        field.appendChild(target)

        # Add the field to fields
        self.fields.appendChild(field)

    def remove_configure_outputs(self):

        # Remove some attributes from the observation element
        configure_attrib_list = ['configure_version', 'seed']

        for attrib in configure_attrib_list:
            if attrib in self.configure.attributes.keys():
                self.configure.removeAttribute(attrib)

        # Remove some attributes from the target elements
        target_attrib_list = [
            'automatic', 'configid', 'fibreid', 'ifu_pa', 'ifu_spaxel',
            'targx', 'targy'
        ]

        # Remove target attributes from field template
        for node in self.dom.getElementsByTagName('target'):
            for attrib in target_attrib_list:
                if attrib in node.attributes.keys():
                    node.removeAttribute(attrib)
                    # Get rid of "%%%" values in photometry
                    for photometry in node.getElementsByTagName('photometry'):
                        for attr_name in photometry.attributes.keys():
                            if photometry.getAttribute(attr_name) == '%%%':
                                photometry.setAttribute(attr_name, '')

        # Remove some elements
        elem_list = [
            'optical_axis', 'distortion_coefficients',
            'telescope', 'focal_plane_map', 'hour_angle_limits'
        ]

        self._remove_elem_list_by_name(elem_list)

    def clean_target_list(self):

        num_removed_targets = 0

        for field in self.fields.getElementsByTagName('field'):
            for target in field.getElementsByTagName('target'):
                if str(target.getAttribute('targcat')) == '%%%':
                    num_removed_targets += 1
                    field.removeChild(target)
                else:
                    target.setAttribute('ifu_spaxel', '')
                    target.setAttribute('ifu_pa', '')
                    target.setAttribute('automatic', '0')

                if str(target.getAttribute('targuse')) == 'S':
                    for photometry in target.getElementsByTagName('photometry'):
                        target.removeChild(photometry)

        if num_removed_targets > 0:
            logging.info('{} non-allocated targets removed'.format(
                         num_removed_targets))

    def _set_epoch(self, epoch):
        self.conditions = self.dom.getElementsByTagName('conditions')[0]

        conditions_dict = {'epoch': '{:.2f}'.format(epoch.byear),
                           'ha': '', 'pressure': '', 'relative_humidity': '',
                           'temperature': '', 'tlr': ''}
        self._set_attribs(self.conditions, conditions_dict)

    def _add_tl_guide_table_to_target_element(self, table, field_template, assert_targuse=None,
                                              add_cnames=True, cache=False):

        """Adds guide stars to the XML file filling in target attributes"""

        col_to_attrib_target_dict = {
                                     'cname': 'CNAME',
                                     'targcat': 'TARGCAT',
                                     'targclass': 'TARGCLASS',
                                     'targdec': 'GAIA_DEC',
                                     'targepoch': 'GAIA_EPOCH',
                                     'targid': 'TARGID',
                                     'targparal': 'GAIA_PARAL',
                                     'targpmdec': 'GAIA_PMDEC',
                                     'targpmra': 'GAIA_PMRA',
                                     'targprio': 'TARGPRIO',
                                     'targprog': 'TARGPROG',
                                     'targra': 'GAIA_RA',
                                     'targsrvy': 'TARGSRVY',
                                     'targuse': 'TARGUSE',
                                     'fibreid': 'FIBREID'  # Added fibre id for pre-configured
                                                           # guide fibres
                                    }

        col_to_attrib_photometry_dict = {
                                         'emag_bp': 'GAIA_MAG_BP_ERR',
                                         'emag_g': 'MAG_G_ERR',
                                         'emag_gg': 'GAIA_MAG_G_ERR',
                                         'emag_i': 'MAG_I_ERR',
                                         'emag_r': 'MAG_R_ERR',
                                         'emag_rp': 'GAIA_MAG_RP_ERR',
                                         'mag_bp': 'GAIA_MAG_BP',
                                         'mag_g': 'MAG_G',
                                         'mag_gg': 'GAIA_MAG_G',
                                         'mag_i': 'MAG_I',
                                         'mag_r': 'MAG_R',
                                         'mag_rp': 'GAIA_MAG_RP',
                                        }

        self.assert_errors = 0

        target_template = self.fields.getElementsByTagName('target')[0].cloneNode(True)
        # print('LAST', table[['GAIA_ID', 'GAIA_RA', 'GAIA_DEC', 'TARGID', 'FIBREID']])

        for row in table:
            target = target_template.cloneNode(True)
            target.setAttribute('fibreid', row['FIBREID'])

            target_attrib_dict = {key: row[col] for key, col in
                                  col_to_attrib_target_dict.items()}

            self._set_attribs(target, target_attrib_dict)

            photometry = target.getElementsByTagName('photometry')[0]
            photometry_attrib_dict = {key: row[col] for key, col in
                                      col_to_attrib_photometry_dict.items()}
            self._set_attribs(photometry, photometry_attrib_dict)

            field_template.appendChild(target)
            # print(target.toxml())
        self.fields.appendChild(field_template)

        if self.assert_errors:
            logging.warning(
                f"There were {self.assert_errors} instances of missing keys "
                f"in the target attribute lookup"
            )

    def _set_tl_guide_stars(self, actual_pa, guides_table, field_template,
                            add_cnames=True, cache=False):

        # Update the PA value if needed

        pa = self._get_pa()
        obsmode = self._get_obsmode()

        # If the observing mode is MOS or mIFU, and the requested pa is NaN (%%%)
        # then we should respect this, and allow Configure to set the optimum PA

        if actual_pa != pa:
            if (obsmode == 'LIFU') or (not np.isnan(pa)):

                if not np.isnan(pa):
                    logging.info('Requested value for PA has NOT been adopted')
                self._set_pa(actual_pa)

        # Add the guide stars to the XML file

        if len(guides_table) == 0:
            logging.error('There are no guide stars available')
            raise SystemExit(2)
        # print('matched guide table per field:')
        # print(guides_table[['GAIA_ID', 'GAIA_RA', 'GAIA_DEC', 'TARGID', 'FIBREID']])
        # print('---------------------------------------')
        self._add_tl_guide_table_to_target_element(guides_table,
                                                   field_template,
                                                   assert_targuse='G',
                                                   add_cnames=add_cnames,
                                                   cache=cache)

    def _get_guide_stars(self, plot_filename=None, useful_table_filename=None,
                         mos_num_stars_request=8, mos_num_central_stars=1,
                         mos_min_cut=0.8, mos_max_cut=0.95, dyn_range_lim=2.0,
                         cache=False, mag_lim=99.):

        obsmode = self._get_obsmode()

        pa = self._get_pa()
        # max_guide = self._get_max_guide()
        resolution = get_resolution_from_progtemp(self._get_progtemp(), cache=False)

        guide_stars = _GuideStars(self.central_ra, self.central_dec, obsmode, resolution)

        actual_pa, guides_table = guide_stars.get_table(
            verbose=True, plot_filename=plot_filename,
            useful_table_filename=useful_table_filename, pa_request=pa,
            num_stars_request=mos_num_stars_request,
            num_central_stars=mos_num_central_stars,
            min_cut=mos_min_cut, max_cut=mos_max_cut,
            cache=cache, mag_lim=mag_lim)

        # Trim full table to MOS suitable candidates
        useful_table = guides_table[guides_table['TARGPROG'] == 'ANY']

        # Do central / annulus filtering with revised list
        central_star_table, non_selected_table = (
                guide_stars._select_central_stars_in_plate_a_b(
                    useful_table,
                    num_central_stars=mos_num_central_stars))

        # Deleted the MOS approach as fibre reachability is already handeled by TL field
        # creation process
        selected_table = vstack([central_star_table, guides_table])

        return actual_pa, selected_table

    def _add_tl_guide_stars(self, guide_star_frame, field_template, plot_filename=None,
                            useful_table_filename=None,
                            mos_num_stars_request=8, mos_num_central_stars=1,
                            mos_min_cut=0.3, mos_max_cut=1.0, dyn_range_lim=2.0,
                            cache=False, mag_lim=99., verbose=False):

        actual_pa, guides_table = self._get_guide_stars(
            plot_filename=plot_filename,
            useful_table_filename=useful_table_filename,
            mos_num_stars_request=mos_num_stars_request,
            mos_num_central_stars=mos_num_central_stars,
            mos_min_cut=mos_min_cut, mos_max_cut=mos_max_cut,
            dyn_range_lim=dyn_range_lim, cache=cache,
            mag_lim=mag_lim)

        # compare guide stars to guide_star csv
        # if guide stars are in csv, add them to the xml
        guides_star_frame_WG = guides_table.to_pandas()
        # print(self.guide_star_batch['GAIA_ID'])
        # print(guides_star_frame_WG['GAIA_ID'].head(5))
        # Check if assignd guide stars are in the query selection
        match_frame = guides_star_frame_WG[guides_star_frame_WG['GAIA_ID'].isin(
                                           self.guide_star_batch['GAIA_ID'])]
        # print(match_frame[['GAIA_ID', 'GAIA_RA', 'GAIA_DEC']])

        # Drop duplicates (central guide star) TODO: Fix this, so drop duplicates is not needed
        match_frame = match_frame.drop_duplicates(subset='GAIA_ID',
                                                  keep='first').reset_index(drop=True)
        # print(match_frame[['GAIA_ID', 'GAIA_RA', 'GAIA_DEC']]),
        if verbose:
            print('Input guide stars:')
            print(self.guide_star_batch[['GAIA_ID', 'fieldId',
                                        'guide', 'fibreId']].sort_values(by='fibreId'))
        match_frame = pd.merge(match_frame, self.guide_star_batch[['GAIA_ID', 'fibreId']],
                               on='GAIA_ID', how='left')

        # match_frame = match_frame.drop(columns=['GAIA_ID'])
        match_frame = match_frame.rename(columns={'fibreId': 'FIBREID'})

        if verbose:
            print('Allocated guide stars:')
            print(match_frame[['GAIA_ID', 'TARGID',
                               'GAIA_RA', 'GAIA_DEC',
                               'FIBREID']].sort_values(by='FIBREID'))

        matched_table = Table.from_pandas(match_frame)
        self._set_tl_guide_stars(actual_pa, matched_table, field_template)

    def _get_calib_stars(self, resolution, plot_filename=None,
                         useful_table_filename=None,
                         num_stars_request=20, num_central_stars=0,
                         min_cut=0.0, max_cut=1.0,
                         cache=False, mag_lim=99., calib_mixres=False):

        obsmode = self._get_obsmode()
        resolution = get_resolution_from_progtemp(self._get_progtemp(), cache=False)

        calib_stars = _CalibStars(self.central_ra, self.central_dec, obsmode, resolution)

        calibs_table = calib_stars.get_table(
            verbose=True, plot_filename=plot_filename,
            useful_table_filename=useful_table_filename,
            num_stars_request=num_stars_request,
            num_central_stars=num_central_stars,
            min_cut=min_cut, max_cut=max_cut,
            cache=cache, mag_lim=mag_lim, calib_mixres=calib_mixres)

        # Do central / annulus filtering with revised list
        central_star_table, non_selected_table = (
                calib_stars._select_central_stars_in_plate_a_b(
                    calibs_table,
                    num_central_stars=num_central_stars))

        if num_stars_request is not None:
            num_ring_stars = num_stars_request - num_central_stars
        else:
            num_ring_stars = None

        ring_stars_table = calib_stars._select_stars_in_ring_of_plate_a_b(
            non_selected_table, num_ring_stars,
            min_cut=min_cut, max_cut=max_cut)

        selected_table = vstack([central_star_table, ring_stars_table])

        return selected_table

    def _add_calib_stars(self, resolution,
                         plot_filename=None, useful_table_filename=None,
                         num_stars_request=2, num_central_stars=0, min_cut=0.2,
                         max_cut=0.4, add_cnames=True, cache=False,
                         mag_lim=99., calib_mixres=False):

        obsmode = self._get_obsmode()

        if obsmode == 'LIFU':
            # No calibration stars needed!
            return

        calibs_table = self._get_calib_stars(resolution,
                                             plot_filename=plot_filename,
                                             useful_table_filename=useful_table_filename,
                                             num_stars_request=num_stars_request,
                                             num_central_stars=num_central_stars,
                                             min_cut=min_cut, max_cut=max_cut,
                                             cache=cache, mag_lim=mag_lim,
                                             calib_mixres=calib_mixres)

        self._set_calib_stars(calibs_table, add_cnames=add_cnames)

    def add_tl_guide_and_calib_stars(self,
                                     resolution,
                                     guide_star_frame,
                                     field_template,
                                     guide_plot_filename=None,
                                     guide_useful_table_filename=None,
                                     mos_num_guide_stars_request=8,
                                     mos_num_central_guide_stars=1,
                                     mos_min_guide_cut=0.8,
                                     mos_max_guide_cut=1.0,
                                     dyn_range_lim=2.0,
                                     calib_plot_filename=None,
                                     calib_useful_table_filename=None,
                                     num_calib_stars_request=20,
                                     num_central_calib_stars=0,
                                     min_calib_cut=0.0,
                                     max_calib_cut=1.0,
                                     calib_mixres=False,
                                     cache=False,
                                     guide_mag_lim=99.,
                                     calib_mag_lim=99.,
                                     calib_stars=False,
                                     verbose=False,
                                     add_cnames=False):

        self._add_tl_guide_stars(guide_star_frame, field_template,
                                 plot_filename=guide_plot_filename,
                                 useful_table_filename=guide_useful_table_filename,
                                 mos_num_stars_request=mos_num_guide_stars_request,
                                 mos_num_central_stars=mos_num_central_guide_stars,
                                 mos_min_cut=mos_min_guide_cut,
                                 mos_max_cut=mos_max_guide_cut,
                                 dyn_range_lim=dyn_range_lim,
                                 cache=cache,
                                 mag_lim=guide_mag_lim, verbose=verbose)

        if calib_stars:
            self._add_calib_stars(resolution=resolution,
                                  plot_filename=calib_plot_filename,
                                  useful_table_filename=calib_useful_table_filename,
                                  num_stars_request=num_calib_stars_request,
                                  num_central_stars=num_central_calib_stars,
                                  min_cut=min_calib_cut,
                                  max_cut=max_calib_cut,
                                  calib_mixres=calib_mixres,
                                  cache=cache,
                                  mag_lim=calib_mag_lim,
                                  add_cnames=add_cnames)

    def _set_guide_star_batch(self, guide_star_frame):
        self.guide_star_batch = guide_star_frame.reset_index(drop=True).astype(
                                                                              {'GAIA_ID': str})

    def set_tl_exposures(self, expcal_file, exposure_time, progtemp, num_science_exposures=1):

        def read_exp_calib_file(filename):
            import json
            with open(filename, 'r') as f:
                exp_calib_data = json.load(f)
            return exp_calib_data

        def generate_order_list(num_science_exposures):
            # Initial calib exposures
            order_list = ["1", "1", "2", "2"]
            num_fields = self.fields.getElementsByTagName('field').length

            if num_fields > 1:
                # Science exposures starting from order 3
                science_orders = [str(3 + i) for i in range(int(num_science_exposures))]
                # Final calib exposures with order number equal to last science order + 1
                final_order_number = str(3 + num_science_exposures)
            # Special case for one field and one science exposure
            elif (num_fields == 1) and (num_science_exposures == 1):
                science_orders = [str(3)]
                final_order_number = str(3 + (num_science_exposures))
            # Special case for one field and more than one science exposure
            elif num_fields == 1:
                science_orders = [str(3 + i) for i in range(int(num_science_exposures))]
                final_order_number = str(3 + int(num_science_exposures))

            else:
                raise Exception("No fields found in the XML")

            order_list.extend(science_orders)
            order_list.extend([final_order_number, final_order_number])

            # print(order_list, science_orders)
            return order_list, science_orders

        expcal_data = read_exp_calib_file(expcal_file)

        # Generate the dynamic order_list
        order_list, science_orders = generate_order_list(num_science_exposures)

        # Get all exposure elements as a list
        original_exposures = list(self.exposures.getElementsByTagName('exposure'))

        # Collect 'science' exposures and calib exposures
        science_exposures = [exp for exp in original_exposures
                             if exp.getAttribute('type') == 'science']
        calibration_exposures = [exp for exp in original_exposures
                                 if exp.getAttribute('type') != 'science']

        # Ensure there are enough calib exposures
        total_calibration_needed = (
            order_list.count("1") +
            order_list.count("2") +
            order_list.count(str(3 + num_science_exposures))
        )

        if len(calibration_exposures) < total_calibration_needed:
            raise Exception(
                f"Not enough calibration exposures. Required: {total_calibration_needed}, "
                f"Available: {len(calibration_exposures)}")

        # Ensure there is at least one science exposure to copy from
        if len(science_exposures) == 0:
            raise Exception("No 'science' exposure found to copy from.")

        # Prepare lists to hold the exposures
        new_exposures = []

        # Index trackers for exposures
        calibration_index = 0
        science_index = 0

        # Assign exposures based on the order_list
        for order_num in order_list:
            if order_num in ["1", "2"]:
                # Assign calibration exposures for orders "1" and "2"
                if calibration_index >= len(calibration_exposures):
                    raise Exception("Not enough calibration exposures for initial orders.")
                exposure = calibration_exposures[calibration_index]
                calibration_index += 1
            elif int(order_num) >= 3 and int(order_num) <= (2 + num_science_exposures):
                # Assign science exposures
                if science_index < len(science_exposures):
                    exposure = science_exposures[science_index]
                    science_index += 1
                else:
                    # If we need more science exposures, copy from the first one
                    exposure = science_exposures[0].cloneNode(deep=True)
                # Update science exposure attributes
                exposure.setAttribute('speed', 'slow')
                exposure.setAttribute('exp_time', str(int(exposure_time)))
                exposure.setAttribute('arm', 'both')
            else:
                # Assign calibration exposures for final orders
                if calibration_index >= len(calibration_exposures):
                    raise Exception("Not enough calibration exposures for final orders.")
                exposure = calibration_exposures[calibration_index]
                calibration_index += 1

            # Set the order attribute
            exposure.setAttribute('order', order_num)

            # Update calibration exposure attributes if necessary
            if exposure.getAttribute('type') != 'science':
                _type = exposure.getAttribute('type')
                _arm = exposure.getAttribute('arm')
                _exp_data = expcal_data[progtemp[0]]
                _binning = progtemp[-1]

                if _type in _exp_data:
                    exp_match = str(_exp_data[_type][_binning][_arm])
                    exposure.setAttribute('exp_time', str(int(exp_match)))

                    # Copy additional attributes
                    for key in _exp_data[_type][_binning]:
                        if key not in ['red', 'blue']:
                            if exposure.getAttribute(key) == '%%%':
                                exposure.setAttribute(key, _exp_data[_type][_binning][key])

            new_exposures.append(exposure)

        # Remove all existing exposures
        while self.exposures.hasChildNodes():
            self.exposures.removeChild(self.exposures.firstChild)

        # Append the new exposures to the 'exposures' element
        for exposure in new_exposures:
            self.exposures.appendChild(exposure)

        self.science_orders = science_orders

    def clone_tl_fields(self, num_clones=1):
        # Create a static copy of the initial field list
        field_list = self.fields.getElementsByTagName('field')
        # print('NUM fields', len(field_list))
        # print('CLONES per field', num_clones)
        if len(field_list) > 1:
            for field in field_list:
                # Find the reference point to insert after the current field
                next_sibling = field.nextSibling
                for i in range(int(num_clones - 1)):
                    # Create a fresh clone for each new field
                    field_template = field.cloneNode(True)
                    # Insert the clone directly after the original field
                    if next_sibling:
                        self.fields.insertBefore(field_template, next_sibling)
                    else:
                        self.fields.appendChild(field_template)
        # Single fields with single exposures must not be cloned!
        # Order is set to empty string to be treated as single exposure
        elif (len(field_list) == 1) and (num_clones == 1):
            field_list[0].setAttribute('order', '')
        elif (len(field_list) == 1) and (num_clones > 1):
            for i in range(int(num_clones - 1)):
                field_template = field_list[0].cloneNode(True)
                self.fields.appendChild(field_template)

    def set_tl_field_orders_for_configure(self):
        field_list = self.fields.getElementsByTagName('field')
        for i, field in enumerate(field_list):
            field.setAttribute('order', str(i + 1))

    def set_tl_field_orders(self):
        field_list = self.fields.getElementsByTagName('field')
        if len(field_list) > 1:
            field_order = 3
            for field in field_list:
                field.setAttribute('order', str(field_order))
                field_order += 1
        elif len(field_list) == 1:
            field_list[0].setAttribute('order', '')

    def set_coords_in_observation(self):
        field_list = self.fields.getElementsByTagName('field')
        ra = field_list[0].getAttribute('RA_d')
        dec = field_list[0].getAttribute('Dec_d')
        self.observation.setAttribute('RA_d', ra)
        self.observation.setAttribute('Dec_d', dec)

        return float(ra), float(dec)

    def set_field_offsets_o(self, ra0, dec0):
        field_list = self.fields.getElementsByTagName('field')
        offsets_ra = []
        offsets_dec = []

        for field in field_list:
            ra1 = float(field.getAttribute('RA_d'))
            dec1 = float(field.getAttribute('Dec_d'))

            delta_ra = ((ra1 - ra0) * 3600) * np.cos(np.radians(dec0))
            delta_dec = (dec1 - dec0) * 3600

            offsets_ra.append(np.round(delta_ra, decimals=1))
            offsets_dec.append(np.round(delta_dec, decimals=1))

        offsets_ra_str = " ".join(map(str, offsets_ra))
        offsets_dec_str = " ".join(map(str, offsets_dec))

        # Decide with Sergio what value we should take. Default == 0
        # self.dithering.setAttribute('apply_dither', "1")

        offsets = self.observation.getElementsByTagName('offsets')[0]
        offsets.setAttribute('offset_step_ra', offsets_ra_str)
        offsets.setAttribute('offset_step_dec', offsets_dec_str)

        return np.sqrt(delta_ra**2 + delta_dec**2)/3600

    def set_field_offsets(self, ra0, dec0):
        """
        Calculate offsets in RA and Dec (in arcseconds) for each 'field' element
        using Astropy's spherical_offsets_to, then store them in the XML.
        Finally, return the separation in degrees for the last field.
        """
        # Prepare lists to store offsets
        offsets_ra = []
        offsets_dec = []

        # Create the reference SkyCoord (in degrees)
        c0 = SkyCoord(ra0*u.deg, dec0*u.deg)

        field_list = self.fields.getElementsByTagName('field')

        for field in field_list:
            # Parse RA/Dec for this field (assumed in degrees)
            ra1 = float(field.getAttribute('RA_d'))
            dec1 = float(field.getAttribute('Dec_d'))

            # Create a SkyCoord for the current field
            c1 = SkyCoord(ra1*u.deg, dec1*u.deg)

            # Use spherical_offsets_to to get the RA/Dec offsets in the local tangent plane
            d_ra, d_dec = c0.spherical_offsets_to(c1)  # returns two Angle objects

            # Convert offsets to arcseconds
            d_ra_arcsec = d_ra.to(u.arcsec).value
            d_dec_arcsec = d_dec.to(u.arcsec).value

            # Round to 1 decimal place
            offsets_ra.append(np.round(d_ra_arcsec, 1))
            offsets_dec.append(np.round(d_dec_arcsec, 1))

        # Convert offset arrays to space-separated strings
        offsets_ra_str = " ".join(map(str, offsets_ra))
        offsets_dec_str = " ".join(map(str, offsets_dec))

        # Store in XML attributes
        offsets = self.observation.getElementsByTagName('offsets')[0]
        offsets.setAttribute('offset_step_ra', offsets_ra_str)
        offsets.setAttribute('offset_step_dec', offsets_dec_str)

        # Return separation in degrees for the *last* field processed
        # (mimicking np.sqrt(delta_ra**2 + delta_dec**2)/3600 for the last field)
        if field_list:
            return c0.separation(c1).deg
        else:
            print("No fields found in the XML")

    def set_ob_class(self, ob_class):
        self.observation.setAttribute('ob_class', ob_class)

    def set_casu_id(self, casu_id):
        self.observation.setAttribute('casu_id', casu_id)


class TLDriverCat:

    def __init__(self, filename, guide_star_csv, expcal_json='./resources/expcal.json',
                 targcat=None, xml_template=None, progtemp_file=None, obstemp_file=None):

        file_list = list(np.array(filename, ndmin=1))

        # Save the filename
        self.filenames = file_list

        # Set the default prefix to be used for output XMLs

        self.default_prefix = os.path.splitext(
            os.path.basename(self.filenames[0]))[0].replace('-mos_driver_cat', '')

        if self.default_prefix[-1] != '-':
            self.default_prefix = self.default_prefix + '-'

        # Set the targcat value
        if targcat is None:
            self.targcat = os.path.basename(self.filenames[0]).replace(
                '-mos_driver_cat', '')
        else:
            self.targcat = targcat

        # Read and save the information into the MOS driver catalogue
        header_keywords = ['DATAMVER', 'TRIMESTE', 'CAT_MAIL',
                           'CAT_CC', 'TACALLOC', 'TACID']

        self.header_values = {}
        for key in header_keywords:
            self.header_values[key] = []

        data = []

        for filename in self.filenames:
            logging.info('Reading {}'.format(filename))
            with fits.open(filename) as hdu_list:
                for key in header_keywords:
                    try:
                        self.header_values[key].append(hdu_list[0].header[key])

                    except Exception as e:
                        logging.warning(f'"{key}" FITS header keywords not present -' +
                                        'submission may fail WASP'
                                        f"(Exception: {e})")
                        self.header_values[key] = ['']

                data.append(Table(hdu_list[1].data))
                # print(Table(hdu_list[1].data).columns)

        self.data = vstack(data)  # Merge catalogues

        # Set guide_star_csv pandas dataframe
        self.guide_star_frame = pd.read_csv(guide_star_csv, dtype={'NAME': str})

        # Multiple input catalogue consistency checks
        if len(self.filenames) > 1:
            for key in ['TRIMESTE', 'TACALLOC', 'TACID']:
                try:
                    assert (np.array(self.header_values[key]) ==
                            self.header_values[key][0]).all()

                except Exception as e:
                    logging.warning(f'"{key}" FITS Header keywords' +
                                    'do not all match:' + f'{self.header_values[key]}'
                                    f"(Exception: {e})")

        if xml_template is None:
            self.xml_template = 'BlankXMLTemplate.xml'
            get_blank_xml_template(self.xml_template, cache=False)
        else:
            self.xml_template = xml_template

        xml_template_dom = xml.dom.minidom.parse(self.xml_template)
        self.xml_datamver = xml_template_dom.childNodes[0].getAttribute('datamver')

        if progtemp_file is None:
            progtemp_file = 'progtemp.dat'
            get_progtemp_file(progtemp_file, cache=False)

        self.progtemp_datamver, self.progtemp_dict, self.forbidden_dict = get_progtemp_dict(
            filename=progtemp_file, assert_orb=True, cache=False)

        if obstemp_file is None:
            obstemp_file = 'obstemp.dat'
            get_obstemp_file(obstemp_file, cache=False)

        self.obstemp_datamver, self.obstemp_dict = get_obstemp_dict(
            filename=obstemp_file, cache=False)

        if expcal_json is not None:
            self.expcal_json = expcal_json

# #############################################################################################
# ########################## Newly added tl ob creation main function #########################

    def process_tl_ob(self,
                      batch_id,
                      target_csv,
                      progtemp,
                      obstemp,
                      plate,
                      h_suffix,
                      exposure_time=1200,
                      calib_stars=True,
                      field_name='WTL_test',
                      epoch=None,
                      report_verbosity=0,
                      linkedgroup_priority=None,
                      linkedgroup='',
                      obsgroup=None,
                      obsgroup_validity=None,
                      radius=1,
                      output_dir='./xmls/',
                      prefix=None,
                      suffix='',
                      avoidance_list=False,
                      targsrvy_dict=None,
                      max_sky=None,
                      max_guide=None,
                      max_calibration=None,
                      num_guide_stars_request=None,
                      dyn_range_lim=3.0,
                      num_central_guide_stars=1,
                      min_guide_cut=0.2,
                      max_guide_cut=0.98,
                      num_calib_stars_request=100,
                      num_central_calib_stars=0,
                      min_calib_cut=0.,
                      max_calib_cut=1.0,
                      calib_mixres=False,
                      gaia=True,
                      tycho=True,
                      twomass=True,
                      lga=True,
                      avoidance_maglim=16.0,
                      avoidance_scale=1,
                      cache=False,
                      guide_mag_lim=10.5,
                      calib_mag_lim=99.,
                      verbose=False,
                      ob_class=None,
                      casu_id=None,
                      # pa=None
                      ):

        self.report_verbosity = int(np.clip(report_verbosity, 0, 1))

        current_time = Time(time.time(), format='unix')
        if epoch is not None:
            self.epoch_time = Time(epoch, format='byear')
        else:
            self.epoch_time = Time(time.time(), format='unix')

        # Verify observing epoch not radiusin past
        assert self.epoch_time >= current_time

        # Guess some parameters which depends on PROGTEMP
        spectrograph_dict = get_progtemp_info(progtemp,
                                              progtemp_dict=None,  # self.progtemp_dict,
                                              cache=False)

        # Change science exposure time for twilight survey. This is a temporary fix!
        spectrograph_dict['red_exp_time'] = exposure_time
        spectrograph_dict['blue_exp_time'] = exposure_time

        obsmode = spectrograph_dict['obsmode']

        assert obsmode in ['MOS']

        assert (spectrograph_dict['red_resolution'] ==
                spectrograph_dict['blue_resolution'])
        resolution = spectrograph_dict['red_resolution']

        red_vph = spectrograph_dict['red_vph']
        blue_vph = spectrograph_dict['blue_vph']

        assert (spectrograph_dict['red_num_exposures'] ==
                spectrograph_dict['blue_num_exposures'])

        # num_science_exposures = spectrograph_dict['red_num_exposures']

        assert (spectrograph_dict['red_exp_time'] ==
                spectrograph_dict['blue_exp_time'])
        # science_exptime = spectrograph_dict['red_exp_time']

        assert (spectrograph_dict['red_binning_x'] ==
                spectrograph_dict['blue_binning_x'])
        binning_x = spectrograph_dict['red_binning_x']

        # Set the spatial binning from the input of this method
        binning_y = 1

        chained = ('+' in progtemp)

        # Set unique OB names for each trimester
        observation_name = self.default_prefix[:-1] + '_' + str(batch_id)
        observation_name = observation_name.replace('-', '_')

        if linkedgroup_priority is not None:
            logging.warning("'linkedgroup_priority' functionality not yet \
                            implemented in OCS. Submission with filled values \
                            may result in WASP Submission Failure")

        if linkedgroup != '':
            linkedgroup_priority = '0'
        else:
            linkedgroup_priority = ''

        if obsgroup == '':
            obsgroup_validity = ''
        else:
            logging.warning("'obsgroup' functionality not yet \
                implemented in OCS. Overriding submission values.")

            obsgroup = ''
            obsgroup_validity = ''

        assert h_suffix in ['N', 'S'], "Invalid hemisphere suffix. Use 'N' or 'S'."

        # Set the position angle based on hemisphere
        if h_suffix == 'N':
            pa = '180'
        elif h_suffix == 'S':
            pa = '0'

        # Prints
        print('-----------------------------------------------------------------------')
        print(f'Processing TL OB batch-ID: {batch_id}; field plate: {plate}; pa: {pa}')
        print('-----------------------------------------------------------------------')

        targsrvy_list = list(set([target['TARGSRVY'] for target in self.data]))
        targsrvy_list = list(filter(None, targsrvy_list))
        targsrvy_list.sort()

        if len(targsrvy_list) > 1:
            logging.warning('Multiple TARGSRVY values found in the input catalogue:'
                            f' {targsrvy_list}')

        # Guess some parameters from obstemp
        obsconstraints_dict = get_obstemp_info(obstemp,
                                               obstemp_dict=self.obstemp_dict)

        elevation_min = obsconstraints_dict['elevation_min']
        moondist_min = obsconstraints_dict['moondist_min']
        seeing_max = obsconstraints_dict['seeing_max']
        skybright_max = obsconstraints_dict['skybright_max']
        transparency_min = obsconstraints_dict['transparency_min']

        # Create XML files for TL-CONFIGURE
        central_target_table = Table.read(target_csv)
        field_batch = central_target_table[central_target_table['batchId'] == batch_id]

        # num_fields_in_batch = len(np.unique(field_batch['fieldId'], return_counts=True))
        num_fields_in_batch = len(field_batch['fieldId'])

        assert num_fields_in_batch >= 1
        # print('NUMSC before', num_science_exposures)
        # print('NUMF before', num_fields_in_batch)

        # exposure_times_fields = num_science_exposures * num_fields_in_batch
        # print('NUMSC', exposure_times_fields)
        # print('NUMF', num_fields_in_batch)

        field_centres = [SkyCoord(ra*u.degree, dec*u.degree, frame='icrs')
                         for ra, dec in zip(field_batch['GAIA_RA'], field_batch['GAIA_DEC'])]

        targprogs = set(field_batch['TARGPROG'])
        tarprog_mask = np.array([targprog in targprogs
                                 for targprog in self.data['TARGPROG']])
        # print(targprogs)
        # print(self.data[['TARGPROG', 'TARGID', 'TARGUSE']])
        # print(self.data[self.data['TARGUSE'] == 'T'])

        batch_data = self.data[tarprog_mask]
        # print(batch_data[['GAIA_RA', 'TARGUSE']])
        field_ids = field_batch['fieldId'].tolist()

        batch_xml = TLOBXML(self.xml_template)

        self._def_max_calibration = int(batch_xml.configure.getAttribute('max_calibration'))
        self._def_max_guide = int(batch_xml.configure.getAttribute('max_guide'))
        self._def_max_sky = int(batch_xml.configure.getAttribute('max_sky'))
        self._def_pa = batch_xml.observation.getAttribute('pa')

        # Remove the elements and attributes which are configure outputs
        batch_xml.remove_configure_outputs()

        # Remove the non-used elements
        batch_xml.remove_non_used_elements()

        # Set the attributes of the root element
        # ifu_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + '/ifu'
        # ifu_version = 'WIW{}'.format(get_version(ifu_path))

        # mos_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + '/mos'
        # mos_version = 'WMW{}'.format(get_version(mos_path))
        # comment = '{} {}'.format(ifu_version,  mos_version)

        comment = '{} {}'.format('WIWv0.6-dirty',  'WMWv0.1-17-g89d9d5cd-dirty')

        root_attrib_dict = {
            'author': self.header_values['CAT_MAIL'][0],
            'cc_report': ','.join(self.header_values['CAT_CC']),
            'comment': comment,
            'report_verbosity': self.report_verbosity,
            'datamver': self.header_values['DATAMVER'][0]
        }

        batch_xml.set_root_attrib(root_attrib_dict)

        # Set the contents of the spectrograph element
        batch_xml.set_spectrograph(binning_x=binning_x, binning_y=binning_y,
                                   resolution=resolution, red_vph=red_vph,
                                   blue_vph=blue_vph)

        # Set the attributes of the observation element
        observation_attrib_dict = {
            'chained': chained,
            'name': observation_name,
            'obstemp': obstemp,
            'obs_mode': obsmode,
            'progtemp': progtemp,
            'trimester': self.header_values['TRIMESTE'][0],
            'linkedgroup': linkedgroup,
            'linkedgroup_priority': linkedgroup_priority,
            'tac_alloc': self.header_values['TACALLOC'][0],
            'tac_id': '',
            'obsgroup': obsgroup,
            'obsgroup_validity': obsgroup_validity,
        }

        assert pa is not None

        if pa is not None:
            observation_attrib_dict['pa'] = self._def_pa

        batch_xml.set_observation(observation_attrib_dict)

        # Add sequence attribute into <observation> element and get rid of %%%
        batch_xml.observation.setAttribute('sequence', 'WTLScience')
        batch_xml.observation.setAttribute('RA_d', '')
        batch_xml.observation.setAttribute('Dec_d', '')
        batch_xml.observation.setAttribute('pa', '')

        if max_sky is not None:
            logging.warning(
                f"max_sky value of {max_sky} differs from template default of "
                f"{self._def_max_sky} - may result in WASP Submission Failure")
        else:
            max_sky = self._def_max_sky

        if max_guide is not None:
            logging.warning(
                f"max_guide value of {max_guide} differs from template default of "
                f"{self._def_max_guide} -  may result in WASP Submission Failure")
        else:
            max_guide = self._def_max_guide

        if max_calibration is not None:
            logging.warning(
                f"max_calibration value of {max_calibration} differs from "
                f"template default of {self._def_max_calibration} "
                "- may result in WASP Submission Failure")
        else:
            max_calibration = self._def_max_calibration

        assert plate in ['PLATE_A', 'PLATE_B']

        num_sky_fibres = max_sky
        plate_max = {'PLATE_A': 960, 'PLATE_B': 940}

        configure_attrib_dict = {
            'max_calibration': max_calibration,
            'max_guide': max_guide,
            'max_sky': max_sky,
            'num_sky_fibres': num_sky_fibres,
            'plate': plate
        }

        batch_xml._set_epoch(self.epoch_time)

        batch_xml.set_configure(configure_attrib_dict)

        # Set the attributes of the obsconstraints element
        obsconstraints_attrib_dict = {
            'elevation_min': elevation_min,
            'moondist_min': moondist_min,
            'seeing_max': seeing_max,
            'skybright_max': skybright_max,
            'transparency_min': transparency_min
        }

        batch_xml.set_obsconstraints(obsconstraints_attrib_dict)

        dithering_attrib_dict = {
            'apply_dither': 0,
            'dither_version': ""
        }
        batch_xml.set_dithering(dithering_attrib_dict)

        offsets = batch_xml.dom.getElementsByTagName('offsets')[0]
        batch_xml._set_attribs(offsets, {'offset_step_ra': "",
                                         'offset_step_dec': ""})

        # Set the contents of the survey element
        batch_xml.set_surveys(targsrvy_list, plate_max[plate])

        if targsrvy_dict is not None:
            survey_nodes = batch_xml.surveys.getElementsByTagName('survey')

            for survey in survey_nodes:
                sname = survey.getAttribute('name')
                try:
                    survey.setAttribute('max_fibres',
                                        str(np.minimum(plate_max[plate],
                                                       targsrvy_dict[f"{sname}_max_fibres"])))
                    survey.setAttribute('priority',
                                        str(np.clip(targsrvy_dict[f"{sname}_priority"],
                                            0.0, 1.0)))
                except Exception as e:
                    logging.warning(
                        f"TARGSRVY {sname} values not present in dict, leaving as default "
                        f"(Error: {e})")
                    raise

        # Probably not needed
        # self.catmask = field_centre.separation(catalog) < radius*u.deg
        # self.prog_obs_cut = ((self.data['PROGTEMP'] == progtemp) *
        #                     (self.data['OBSTEMP'] == obstemp))

        if pa is not None:
            batch_xml.observation.setAttribute('pa', pa)
        else:
            batch_xml.observation.setAttribute('pa', self._def_pa)

        # Set the contents of the fields element
        for field_centre, field_id in zip(field_centres, field_ids):
            batch_xml.central_ra = field_centre.ra.deg
            batch_xml.central_dec = field_centre.dec.deg

            field_element, targprog = batch_xml.create_and_prepare_fields(
                batch_xml.central_ra,
                batch_xml.central_dec,
                central_target_table
            )

            batch_xml.fields.removeChild(batch_xml.fields.firstChild)

            targets_targprogs = batch_data['TARGPROG']
            targets_mask = [target_targprog == targprog
                            for target_targprog in targets_targprogs]

            # print('batch data', batch_data[['TARGNAME', 'TARGPROG']])
            filtered_batch_data = batch_data[targets_mask]
            # print('filtered batch data', filtered_batch_data[['TARGNAME', 'TARGPROG']])

            # Add targets to the field
            batch_xml.add_table_as_targets(field_element, filtered_batch_data)

            # Set guide star batch variable (pandas Dataframe) using the current field id
            current_field = self.guide_star_frame[self.guide_star_frame['fieldId'] == field_id]

            batch_xml._set_guide_star_batch(current_field)
            # Add guide stars to the field (Calibration stars will be added later)

            batch_xml.add_tl_guide_and_calib_stars(
                resolution,
                guide_star_frame=batch_xml.guide_star_batch,
                field_template=field_element,
                guide_plot_filename=None,
                guide_useful_table_filename=None,
                mos_num_guide_stars_request=num_guide_stars_request,
                mos_num_central_guide_stars=num_central_guide_stars,
                mos_min_guide_cut=min_guide_cut,
                mos_max_guide_cut=max_guide_cut,
                dyn_range_lim=dyn_range_lim,
                calib_plot_filename=None,
                calib_useful_table_filename=None,
                num_calib_stars_request=num_calib_stars_request,
                num_central_calib_stars=num_central_calib_stars,
                min_calib_cut=min_calib_cut,
                max_calib_cut=max_calib_cut,
                calib_mixres=calib_mixres,
                cache=cache,
                guide_mag_lim=guide_mag_lim,
                calib_mag_lim=calib_mag_lim,
                calib_stars=calib_stars,
                verbose=verbose
            )

        batch_xml.clean_target_list()

        for child in list(batch_xml.fields.childNodes):
            if not child.childNodes:
                batch_xml.fields.removeChild(child)
            elif ((child.getAttribute('RA_d') == '%%%') and
                  (child.getAttribute('Dec_d') == '%%%')):
                batch_xml.fields.removeChild(child)

        # Set field orders for configure
        batch_xml.set_tl_field_orders_for_configure()

        # For testing purposes likely to be removed in final version together with var:ob_class
        if ob_class is not None:
            batch_xml.set_ob_class(ob_class)

        # For testing purposes likely to be removed in final version together with var:casu_id
        if casu_id is not None:
            batch_xml.set_casu_id(casu_id)

        # Write the OB XML to a file
        output_path = self._get_tl_output_path(field_name, output_dir=output_dir,
                                               prefix=prefix, suffix=suffix,
                                               h_suffix=h_suffix)

        batch_xml.write_xml(output_path)

        """
        if avoidance_list:
            from skyavoidance import gen_avoidance_list, write_avoidance_list

            avoidance_list = gen_avoidance_list(ob_xml, gaia, tycho,
                                                twomass, lga,
                                                avoidance_maglim,
                                                avoidance_scale)

            path = os.path.splitext(output_path)[0]+'.sky.xml'
            avoidance_xml = write_avoidance_list(avoidance_list, path,
                                                 self.xml_datamver)
        """
        return output_path

    def _get_tl_output_path(self, field_name, output_dir='', prefix=None,
                            suffix='', h_suffix=''):

        # Set the prefix of the output file if it has not been provided
        if prefix is None:
            prefix = self.default_prefix
        split_prefix = prefix.split('_')
        prefix = split_prefix[0] + '_' + f'p{split_prefix[-2]}' + '_'

        # print(split_prefix)
        # Choose the first filename which does not exist
        index = 0
        output_path = ''

        while (index < 1) or os.path.exists(output_path):
            index += 1
            output_basename = '{}{}_{:02d}_b{}_{}.xml'.format(
                prefix, field_name, index, str(int(suffix)).zfill(3), h_suffix)
            output_path = os.path.join(output_dir, output_basename)
        return output_path

    def finalize_tl_ob(self, config_ob_xml, progtemp, output_dir,
                       num_science_exposures=1,
                       science_exptime=100,
                       field_name='WTLS',
                       prefix=None,
                       xml_verbosity=False,
                       drop_targinfo_from_field=False):

        batch_xml = TLOBXML(config_ob_xml)
        # set the report_verbosity attribute to 1 in <weave> element
        if xml_verbosity:
            batch_xml.root.setAttribute('report_verbosity', '1')

        for target in batch_xml.fields.getElementsByTagName('target'):
            if target.getAttribute('targname') == "%%%":
                target.setAttribute('targname', '')

        first_field = batch_xml.fields.getElementsByTagName('field')[0]
        batch_id = first_field.getAttribute('targprog').split('_')[0]
        observation_name = self.default_prefix[:-1] + '_' + str(batch_id)
        observation_name = observation_name.replace('-', '_')

        num_fields = batch_xml.fields.getElementsByTagName('field').length
        exposure_times_fields = num_science_exposures * num_fields
        # print('BATCH #:', batch_id)

        # Clone the fields to enable the OCS to process multiple exposures per field
        batch_xml.clone_tl_fields(num_clones=num_science_exposures)
        batch_xml.set_tl_field_orders()

        # Set exposures here to use final field count for <exposures>
        batch_xml.set_tl_exposures(expcal_file=self.expcal_json,
                                   exposure_time=science_exptime,
                                   progtemp=progtemp,
                                   num_science_exposures=exposure_times_fields)

        # Set the initial coordinates of the observation and calculate the field offsets
        initial_ra, intitial_dec = batch_xml.set_coords_in_observation()
        offset = batch_xml.set_field_offsets(initial_ra, intitial_dec).round(2)
        offset = str(offset).replace('.', '_')

        # Add "M" for "merged" to the OB name. Note that single fields will also be named
        # with the "merged" suffix.
        batch_xml.observation.setAttribute('name', observation_name + '_M')

        configure_element = batch_xml.dom.getElementsByTagName('configure')[0]
        configure_element.setAttribute('seed', '')

        # This below is currently a fix to enable WASP validation of the OBs
        # Need to discuss if we want to add twilight_configure_version attribute,
        # or agree on a tl-configure version number recognized by WASP.
        # Rename the twilight_configure_version attribute to configure_version:
        configure_element.setAttribute(
            'configure_version',
            configure_element.getAttribute('twilight_configure_version')
        )
        configure_element.removeAttribute('twilight_configure_version')

        # get configure version value and remove the "t" character
        configure_version = configure_element.getAttribute('configure_version')
        configure_version = configure_version.split('t')[0]
        configure_element.setAttribute('configure_version', configure_version)

        # get epoch which is in the subelement conditions of the batch_xml.observation elemnt
        epoch = (
            batch_xml.observation.getElementsByTagName('conditions')[0].getAttribute('epoch')
        )

        epoch = str(np.round(float(epoch), 3)).replace('.', '_')

        pa_value = batch_xml.observation.getAttribute('pa')
        if pa_value == '0':
            hemisphere = '_S'
        elif pa_value == '180':
            hemisphere = '_N'

        if drop_targinfo_from_field:
            for field in batch_xml.fields.getElementsByTagName('field'):
                field.removeAttribute('targprog')
                field.removeAttribute('targid')

        field_name = field_name + hemisphere

        h_suffix = f"e{epoch}_exp{num_science_exposures}_t{science_exptime}_o{offset}"

        # Write the OB XML to a file
        output_path = self._get_tl_output_path(field_name, output_dir=output_dir,
                                               prefix=prefix, suffix=batch_id,
                                               h_suffix=h_suffix)

        batch_xml.write_xml(output_path)


# #############################################################################################

# ############## Newly added functions to create target/guidestar catalogues ##################


def create_central_target_csv(catalogue_csv, plate, hemisphere=None, ra_col='ra', de_col='dec',
                              mag_col='phot_g_mean_mag', src_id='dr3_source_id',
                              output_dir='./output/'):

    # Read and preprocess the catalogue
    if type(catalogue_csv) is pd.DataFrame:
        cat_frame = catalogue_csv
    elif isinstance(catalogue_csv, str):
        cat_frame = pd.read_csv(catalogue_csv).astype({src_id: str})

    cat_frame = cat_frame[cat_frame['guide'] == 1]
    cat_frame = cat_frame.dropna(subset=['fibreId'])
    cat_frame = cat_frame[[src_id, ra_col, de_col, mag_col, 'GAIA_ID',
                           'guide', 'fibreId', 'batchId',
                           'fieldId', 'TARGID', 'batchScore']]

    cat_frame[['guide', 'fibreId', 'batchScore']] = cat_frame[['guide',
                                                               'fibreId',
                                                               'batchScore']].astype(int)
    if type(catalogue_csv) is pd.DataFrame:
        full_frame = catalogue_csv
    elif isinstance(catalogue_csv, str):
        # Read the full catalogue to count guide == 2 entries
        full_frame = pd.read_csv(catalogue_csv)

    # Count the number of entries with guide == 2 for each fieldId
    count_guide_2 = full_frame[full_frame['guide'] == 2
                               ].groupby(
                                'fieldId').size().reindex(cat_frame['fieldId'].unique(),
                                                          fill_value=0)

    # Map the counts back to the original DataFrame
    cat_frame['count_guide_2'] = cat_frame['fieldId'].map(count_guide_2)

    # Sort the DataFrame by batchId and then by the count of guide == 2 in descending order
    cat_frame = cat_frame.sort_values(by=['batchId', 'count_guide_2', 'fieldId'],
                                      ascending=[True,
                                                 False,
                                                 True]).drop(columns='count_guide_2')

    # Apply assign_targprog function
    cat_frame = assign_targprog(catalogue_df=cat_frame)

    # Rename columns
    # cat_frame = cat_frame.rename(columns={src_id: 'NAME', ra_col: 'RA', de_col: 'DEC',
    #                                     mag_col: 'GAIA_MAG_G'})

    # Ensure the output directory exists
    # output_dir = output_dir
    # output_dir.mkdir(parents=True, exist_ok=True)

    # Save to CSV
    if hemisphere is not None:
        csv_name = f'central_target_list_{plate}_{hemisphere}.csv'
    else:
        csv_name = f'central_target_list_{plate}.csv'

    cat_frame.to_csv(output_dir + csv_name, index=False)


def create_guide_stars_csv(catalogue_csv, plate, hemisphere=None, ra_col='ra', de_col='dec',
                           src_id='dr3_source_id', output_dir='./output/'):

    if type(catalogue_csv) is pd.DataFrame:
        cat_frame = catalogue_csv
    elif isinstance(catalogue_csv, str):
        cat_frame = pd.read_csv(catalogue_csv).astype({src_id: str})

    cat_frame = cat_frame[cat_frame['guide'].isin([1, 2])]
    cat_frame = cat_frame.dropna(subset=['fibreId'])
    cat_frame = cat_frame[[src_id, ra_col, de_col, 'GAIA_ID',
                           'guide', 'fibreId', 'batchId', 'fieldId']]

    cat_frame[['guide', 'fibreId', 'batchId', 'fieldId']] = cat_frame[['guide',
                                                                       'fibreId',
                                                                       'batchId',
                                                                       'fieldId']].astype(int)
    # cat_frame = cat_frame.rename(columns={src_id: 'NAME', ra_col: 'RA',
    #                                     de_col: 'DEC'})

    if hemisphere is not None:
        csv_name = f'guide_stars_{plate}_{hemisphere}.csv'
    else:
        csv_name = f'guide_stars_{plate}.csv'

    cat_frame.to_csv(output_dir + csv_name, index=False)


def assign_obstemp(catalogue_df,
                   obstemp):
    """
    Assign correct OBSTEMP to science targets and guide stars (sourceFlag=9).

    Args:
        catalogue_df (pd.DataFrame): DataFrame containing the catalogue
        obstemp (str): OBSTEMP value to assign

    Returns:
        pd.DataFrame: Updated DataFrame with assigned OBSTEMP
    """
    catalogue_df['obstemp'] = obstemp
    catalogue_df.loc[catalogue_df['sourceFlag'] == 9, 'obstemp'] = 'XXXXX'
    return catalogue_df


def assign_progtemp(catalogue_df,
                    progtemp):
    """
    Assign correct PROGTEMP to science targets and guide stars (sourceFlag=9).

    Args:
        catalogue_df (pd.DataFrame): DataFrame containing the catalogue
        progtemp (str): PROGTEMP value to assign

    Returns:
        pd.DataFrame: Updated DataFrame with assigned PROGTEMP
    """

    catalogue_df['progtemp'] = str(progtemp)
    catalogue_df.loc[catalogue_df['sourceFlag'] == 9, 'progtemp'] = '00000.0'
    return catalogue_df


def assign_targprog(catalogue_df):
    """
    Assigns TARGPROG for every target based on batchId and fieldId. This specific TARGPROG
    will only be present in the XML files, and remains unchanged in the input catalogue!

    Example: batchId = 1 and fieldId = 42 -> TARGPROG = '1_42', dtype=str

    Args:
        catalogue_df (pd.DataFrame): DataFrame containing the catalogue

    Returns:
        pd.DataFrame: Updated DataFrame with assigned TARGPROG
    """

    batch_ids = catalogue_df['batchId'].astype(int)
    field_ids = catalogue_df['fieldId'].astype(int)
    catalogue_df.loc[:, 'TARGPROG'] = batch_ids.astype(str) + '_' + field_ids.astype(str)
    return catalogue_df


def assign_target_attributes(catalogue_df,
                             src_id='dr3_source_id',
                             conditions=None):
    """
    Annotate catalogue rows with WEAVE target attributes (`targclass`, `targuse`,
    `targprio`, `targname`, `targid`) based on flags and optional priority rules.
    Currently this function limits itself to only creating targclass = "STAR" targets.
    If necessary this will be changed in the future.

    By default:
      - `targclass` is 'STAR' for all, overridden to 'SKY' where `sourceFlag == -1`.
      - `targuse` is 'T' (science), overridden to 'G' for guides (`sourceFlag == 9`)
        and 'S' for sky (`sourceFlag == -1`).
      - `targprio` is 1.0 globally, but can be overridden per-value using `conditions`.
      - `targname` and `targid` are populated from the column named by `src_id`.

    Args:
        catalogue_df (pd.DataFrame): Input table that must include:
            - `sourceFlag` (int): classification flag (e.g., 9 = guide, -1 = sky).
            - the identifier column named by `src_id`.
            - any columns referenced in `conditions` (if provided).
        src_id (str, optional): Column used to fill `targname` and `targid`.
            Defaults to 'dr3_source_id'.
        conditions (dict[str, dict[Any, float]] | None, optional): Optional mapping to
            assign custom priorities. Structure:
            `{ column_name: { value: priority, ... }, ... }`.
            For each `(column_name, value)`, rows where `catalogue_df[column_name] == value`
            will receive `targprio = priority`.

    Returns:
        pd.DataFrame: The input DataFrame with added/updated columns:
            - `targclass` ∈ {'STAR', 'SKY'}
            - `targuse` ∈ {'T', 'G', 'S'}
            - `targprio` (float)
            - `targname`, `targid` (copied from `src_id`)

    Notes:
        - Priority assignments in `conditions` overwrite the default `targprio` (1.0)
          and may overwrite one another if multiple rules match the same row; the last
          matching rule applied wins due to iteration order.
        - This function modifies a copy-on-write-aware DataFrame in place w.r.t. logical
          assignment semantics; if you need to avoid mutating the original frame, pass
          `catalogue_df.copy()`.

    Raises:
        KeyError: If `src_id` is not a column in `catalogue_df`, or if any column named
          in `conditions` is missing.
    """

    catalogue_df['targclass'] = 'STAR'
    catalogue_df.loc[catalogue_df['sourceFlag'] == -1, 'targclass'] = 'SKY'

    catalogue_df['targuse'] = 'T'
    catalogue_df.loc[catalogue_df['sourceFlag'] == 9, 'targuse'] = 'G'
    catalogue_df.loc[catalogue_df['sourceFlag'] == -1, 'targuse'] = 'S'

    # Default priority
    catalogue_df['targprio'] = 1.0

    # Set priorities according to conditions
    if conditions:
        for column, value_priority_map in conditions.items():
            for value, priority in value_priority_map.items():
                mask = catalogue_df[column] == value
                catalogue_df.loc[mask, 'targprio'] = priority

    # Assign GAIA source id to targname and targid
    catalogue_df['targname'] = catalogue_df[src_id]
    catalogue_df['targid'] = catalogue_df[src_id]

    return catalogue_df


def create_plate_fits_file(infile,
                           coverage_df,
                           plate='A',
                           hemisphere='N',
                           outpath='./',
                           src_id_col='CNAME',  # src id column name in FITS table to match on
                           coverage_targprog_col='TARGPROG'):
    """
    Create a plate-specific FITS file by filtering and updating the target table,
    while preserving the primary HDU and all other extensions.

    The function:
      - keeps only rows in extension 1 where `src_id_col` is present in
        `coverage_df[src_id_col]`,
      - overwrites (or creates) `TARGPROG` in extension 1 using
        `coverage_df[coverage_targprog_col]`, matched by `src_id_col`,
      - appends “sky” rows from `coverage_df` where `sourceFlag == -1` and the
        `src_id_col` value is blank (these are not matched by ID),
      - writes the result as `<outpath>/<infile_base>_<plate>_<hemisphere>.fits`
        with FITS checksums,
      - and returns the first non-empty unique `PROGTEMP` and `OBSTEMP` strings
        found among non-sky rows.

    If `plate` is not None, helper CSVs for guide stars and central targets are also
    written alongside the FITS using `coverage_df`.

    Args:
        infile (str): Path to the input FITS file. The primary HDU and all
            non-target-table extensions are copied unchanged.
        coverage_df (pd.DataFrame): Coverage catalogue used to filter and annotate the
            target table. Must include `src_id_col`, `coverage_targprog_col`, and (for
            sky rows) `sourceFlag`.
        plate (str, optional): Plate identifier ('A' or 'B'). If not None, CSV helper
            files for guides/central targets are written. Defaults to 'A'.
        hemisphere (str, optional): Hemisphere code ('N' or 'S') passed to the helper
            CSV writers. Defaults to 'N'.
        outpath (str, optional): Output directory for the new FITS and helper CSVs.
            Defaults to './'.
        src_id_col (str, optional): Column name used to match rows (e.g. 'CNAME').
            Defaults to 'CNAME'.
        coverage_targprog_col (str, optional): Column in `coverage_df` containing the
            desired `TARGPROG` values. Defaults to 'TARGPROG'.

    Returns:
        tuple[str, str]:
            - `progtemp_str`: First non-empty unique `PROGTEMP` value among non-sky rows.
            - `obstemp_str`: First non-empty unique `OBSTEMP` value among non-sky rows.

    Notes:
        - Matching is performed on a whitespace-stripped string version of `src_id_col`.
        - When appending sky rows, the function conforms their columns and dtypes to the
          FITS table schema (including byte/unicode handling).
        - Only extension 1 (the main target table) is filtered/updated; all other HDUs
          are preserved as-is.
        - The output filename is derived from `infile` (base name) plus
          `_<plate>_<hemisphere>.fits`.

    Raises:
        TypeError: If `coverage_df` is not a pandas DataFrame.
        KeyError: If required columns (`src_id_col`, `coverage_targprog_col`) are missing
            from `coverage_df`, or if `src_id_col` is missing from the FITS table.
    """

    # Helper functions
    def _series_to_clean_str(s):
        return (s.astype('string')
                 .str.strip()
                 .fillna(''))

    def _to_bytes_array(values, target_dtype):
        return np.asarray([str(v).encode('utf-8') for v in values], dtype=target_dtype)

    def _coerce_to_dtype(values, target_dtype):
        """Coerce a 1D array-like 'values' into the numpy dtype of a FITS column."""
        kind = target_dtype.kind
        if kind == 'S':  # fixed-length bytes
            return _to_bytes_array(values, target_dtype)
        if kind == 'U':  # unicode
            return np.asarray([str(v) for v in values], dtype=target_dtype)
        if kind == 'f':  # float
            arr = pd.to_numeric(values, errors='coerce').to_numpy(dtype=float)
            return arr.astype(target_dtype, copy=False)
        if kind in ('i', 'u'):  # ints
            arr = pd.to_numeric(values, errors='coerce').fillna(0).to_numpy(dtype=np.int64)
            return arr.astype(target_dtype, copy=False)
        if kind == 'b':  # bool
            return np.asarray(values, dtype=bool).astype(target_dtype, copy=False)
        try:
            return np.asarray(values, dtype=target_dtype)
        except Exception:
            return np.asarray([str(v) for v in values], dtype=target_dtype)

    def _unique_nonempty_from_table(table, colname):
        """Return first non-empty unique string value from table[colname],
        preferring rows with sourceFlag != -1 if that column exists."""
        if len(table) == 0 or colname not in table.colnames:
            return ""
        col = table[colname]
        # decode/strip to strings
        if col.dtype.kind == 'S':
            vals = np.char.strip(np.char.decode(col, 'utf-8', errors='ignore'))
        else:
            vals = np.char.strip(col.astype(str))
        nonempty = vals != ''

        if 'sourceFlag' in table.colnames:
            try:
                sf = np.array(table['sourceFlag'])
                mask = nonempty & (sf != -1)
            except Exception:
                mask = nonempty
        else:
            mask = nonempty

        if not np.any(mask):
            return ""
        uniq = np.unique(vals[mask])
        return uniq[0] if len(uniq) > 0 else ""

    # ############################ Main Logic ############################

    infile_name = infile.split("/")[-1].split(".")[0]

    if plate is not None:
        create_guide_stars_csv(coverage_df, plate=plate, ra_col='GAIA_RA',
                               de_col='GAIA_DEC', src_id=src_id_col, hemisphere=hemisphere,
                               output_dir=outpath)
        create_central_target_csv(coverage_df, plate=plate, ra_col='GAIA_RA',
                                  de_col='GAIA_DEC', mag_col='GAIA_MAG_G',
                                  hemisphere=hemisphere,
                                  src_id=src_id_col, output_dir=outpath)

    if not isinstance(coverage_df, pd.DataFrame):
        raise TypeError("coverage_df must be a pandas DataFrame")

    # Normalise coverage data used for matching/mapping
    cov_cnames = _series_to_clean_str(coverage_df[src_id_col])
    cov_targprog = _series_to_clean_str(coverage_df[coverage_targprog_col])

    keep_set = set(cov_cnames.tolist())
    cname_to_targprog = dict(zip(cov_cnames.tolist(), cov_targprog.tolist()))

    # Identify sky rows to append: sourceFlag == -1 and CNAME missing/blank
    sky_mask = (coverage_df.get('sourceFlag', pd.Series([], dtype='Int64')) == -1)
    blank_cname_mask = _series_to_clean_str(coverage_df[src_id_col]).eq('')
    append_sky_df = coverage_df[sky_mask & blank_cname_mask].copy()

    ext = 1

    with fits.open(infile, mode='readonly', memmap=False) as hdul:
        new_hdul = fits.HDUList([hdul[0].copy()])

        for i in range(1, len(hdul)):
            if i != ext:
                new_hdul.append(hdul[i].copy())
                continue

            tbl_hdu = hdul[i]
            data = tbl_hdu.data
            if data is None:
                new_hdul.append(tbl_hdu.copy())
                continue

            tbl = Table(data)

            # CNAME column to bytes
            cname_col_data = tbl[src_id_col]
            if cname_col_data.dtype.kind == 'S':
                norm_cname = np.char.strip(np.char.decode(cname_col_data,
                                                          'utf-8', errors='ignore'))
            else:
                norm_cname = np.char.strip(cname_col_data.astype(str))

            mask_keep = np.isin(norm_cname, list(keep_set))
            tbl_filtered = tbl[mask_keep]

            # If there are kept rows, map/overwrite TARGPROG from coverage_df
            if len(tbl_filtered) > 0:
                if tbl_filtered[src_id_col].dtype.kind == 'S':
                    norm_cname_filtered = np.char.strip(
                        np.char.decode(tbl_filtered[src_id_col], 'utf-8', errors='ignore')
                    )
                else:
                    norm_cname_filtered = np.char.strip(tbl_filtered[src_id_col].astype(str))

                mapped_targprog = np.array(
                    [cname_to_targprog.get(cn, "") for cn in norm_cname_filtered],
                    dtype=object
                )

                # Overwrite TARGPROG in the FITS table.
                # (essential for correct guide fibre allocation,
                # original input FITS not changed!)
                if 'TARGPROG' not in tbl_filtered.colnames:
                    has_bytes = any(tbl_filtered[c].dtype.kind == 'S'
                                    for c in tbl_filtered.colnames)
                    if has_bytes:
                        width = max(16, int(mapped_targprog.astype(str).map(len).max()
                                            if len(mapped_targprog) else 16))
                        dtype = f'S{width}'
                        tbl_filtered['TARGPROG'] = _to_bytes_array(mapped_targprog,
                                                                   np.dtype(dtype))
                    else:
                        tbl_filtered['TARGPROG'] = mapped_targprog.astype(str)
                else:
                    if tbl_filtered['TARGPROG'].dtype.kind == 'S':
                        enc = _to_bytes_array(mapped_targprog, tbl_filtered['TARGPROG'].dtype)
                        tbl_filtered['TARGPROG'] = enc
                    else:
                        tbl_filtered['TARGPROG'] = mapped_targprog.astype(str)

            # Append sky rows (sourceFlag == -1, blank CNAME) at the end
            if len(append_sky_df) > 0:
                schema_ref = tbl_filtered if len(tbl_filtered) > 0 else tbl
                names = list(schema_ref.colnames)
                dtypes = [schema_ref[c].dtype for c in names]
                n_new = len(append_sky_df)

                col_data = {}
                for c, dt in zip(names, dtypes):
                    if c in append_sky_df.columns:
                        src_vals = append_sky_df[c]
                    elif c == 'TARGPROG' and coverage_targprog_col in append_sky_df.columns:
                        src_vals = append_sky_df[coverage_targprog_col]
                    elif c == src_id_col:
                        src_vals = pd.Series([''] * n_new)
                    else:
                        src_vals = pd.Series([np.nan] * n_new)

                    if len(src_vals) != n_new:
                        src_vals = pd.Series(list(src_vals)
                                             + [np.nan] * (n_new - len(src_vals)))[:n_new]

                    col_data[c] = _coerce_to_dtype(src_vals, dt)

                sky_tbl = Table(col_data)

                if len(tbl_filtered) > 0:
                    tbl_out = vstack([tbl_filtered, sky_tbl],
                                     join_type='exact',
                                     metadata_conflicts='silent')
                else:
                    tbl_out = sky_tbl

            else:
                tbl_out = tbl_filtered

            # Rebuild BinTableHDU with same header/name
            new_tbl = fits.BinTableHDU(
                data=tbl_out.as_array(),
                header=tbl_hdu.header.copy(),
                name=tbl_hdu.name
            )
            new_hdul.append(new_tbl)

        # Write with checksums
        outfile = f"{outpath}{infile_name}_{plate}_{hemisphere}.fits"
        new_hdul.writeto(outfile, overwrite=True, checksum=True)

    # Extract PROGTEMP / OBSTEMP from non-sky, non-empty rows
    out_tbl = Table(new_hdul[ext].data)
    progtemp_str = _unique_nonempty_from_table(out_tbl, 'PROGTEMP')
    obstemp_str = _unique_nonempty_from_table(out_tbl, 'OBSTEMP')

    return progtemp_str, obstemp_str


if __name__ == '__main__':

    # Example usecase:

    plate = 'A'
    hemisphere = 'N'

    # Create fits input catalogue, central_target_list.csv and guide_stars.csv
    progtemp, obstemp = create_plate_fits_file(infile='./WS2025XX-XXX.fits',
                                               coverage_df='./tl_coverage_plate_A_N.csv',
                                               plate=plate,
                                               hemisphere=hemisphere,
                                               src_id_col='CNAME',
                                               outpath='./')

    # Initialise TLDriverCat
    tl_driver = TLDriverCat(filename=f'./WS2025XX-XXX_{plate}_{hemisphere}.fits',
                            guide_star_csv=f'./guide_stars_{plate}_{hemisphere}.csv',
                            expcal_json='./resources/expcal.json')

    # Set paths for central_target_list.csv
    target_csv_path = f'./central_target_list_{plate}_{hemisphere}.csv'

    # Read in as pandas dataframe to get unique field batch ids:
    central_target_df = pd.read_csv(target_csv_path)
    print(central_target_df.shape[0], 'central guide stars found!')
    field_batch_ids = central_target_df['batchId'].unique()

    # Run creation of CONFIGURE xmls at './output_data/xmls'
    for batch_id in field_batch_ids:
        tl_driver.process_tl_ob(batch_id, target_csv_path,
                                h_suffix=hemisphere,
                                progtemp=progtemp,
                                obstemp=obstemp,
                                plate=f'PLATE_{plate}',
                                calib_stars=False,
                                output_dir=f'./xmls/plate_{plate}/{hemisphere}/',
                                guide_mag_lim=16.,
                                calib_mag_lim=20.0,
                                min_guide_cut=0.1,
                                max_guide_cut=0.98,
                                suffix=batch_id,
                                verbose=False,
                                ob_class='science',
                                field_name='DEMO')
