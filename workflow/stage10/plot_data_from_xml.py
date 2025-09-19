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


import argparse
import logging
import os.path
import re
from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
from matplotlib.patches import Ellipse, FancyArrow, FancyArrowPatch
from matplotlib.backends.backend_pdf import PdfPages
from astropy import coordinates

from workflow.utils.get_data_from_xmls import \
    get_obs_mode, get_coord_of_first_field, get_data_per_field_of_one_xml


class _IFUPlot():

    def __init__(self, xml_file, output_dir='', sim=False, verbose=False):

        # Save the input XML filename

        self.xml_file = xml_file

        # Save the simulation and verbose mode

        self.sim = sim
        self.verbose = verbose

        # Set the filename for the PDF file and create its object

        xml_basename_wo_ext = os.path.splitext(os.path.basename(xml_file))[0]

        if self.sim is True:
            sim_str = '-sim'
        else:
            sim_str = ''

        if self.verbose is True:
            verbose_str = '-verbose'
        else:
            verbose_str = ''

        self.pdf_filename = os.path.join(
            output_dir, xml_basename_wo_ext + sim_str + verbose_str + '.pdf')

        self.pdf = PdfPages(self.pdf_filename)

        # Get the obs_mode and assert that it is an IFU XML

        self.obs_mode = get_obs_mode(xml_file)

        assert self.obs_mode in ['LIFU', 'mIFU']
    
        # Get the data from the XML file

        self.data_dict_list = get_data_per_field_of_one_xml(self.xml_file)

        # Set the aspect to be used in the figures

        self.figaspect = 1 / ((1 + 5 ** 0.5) / 2)

        # Set the coordenates of the field

        self.field_ra, self.field_dec = get_coord_of_first_field(xml_file)

        # Set the fibre size, the list with the central spaxels and the field
        # of views

        self._set_fibre_size()
        self._set_central_spaxel_list()
        self._set_fov_list()

        # Make the plots

        self._make_plots()

    def _set_fibre_size(self):
        raise NotImplementedError

    def _set_central_spaxel_list(self):
        raise NotImplementedError

    def _set_fov_list(self):
        raise NotImplementedError

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.pdf.__exit__(exc_type, exc_val, exc_tb)

    def close(self):
        self.pdf.close()

    def get_pdf_filename(self):
        return self.pdf_filename

    def _get_fig_ax_vect(self):
        raise NotImplementedError

    def _set_aspect_on_axes(self, ax_vect):
        for ax in ax_vect:
            ax.set_aspect(1.0 / np.cos(np.deg2rad(self.field_dec)))

    def _set_ax_lims(self, ax_vect):

        data_dict = self.data_dict_list[0]

        for i, ax in enumerate(ax_vect):
            central_spaxel = self.central_spaxel_list[i]
            fov = self.fov_list[i]

            try:
                index = data_dict['ifu_spaxel'].index(central_spaxel)
                ra = data_dict['targra'][index]
                dec = data_dict['targdec'][index]
            except ValueError:
                # If the central bundle is not present, we put the limits out
                # of the range of the coordinates
                ra = 1000
                dec = self.field_dec

            xlim = [ra + (fov / 2) / np.cos(np.deg2rad(self.field_dec)),
                    ra - (fov / 2) / np.cos(np.deg2rad(self.field_dec))]
            ylim = [dec - (fov / 2), dec + (fov / 2)]

            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

    def _clean_ticks_on_axes(self, ax_vect):
        for ax in ax_vect:
            ax.set_xticks([])
            ax.set_yticks([])

    def _add_north_east(self, fig):

        if self.obs_mode == 'LIFU':
            origin_x = 0.66
            origin_y = 0.09
        elif self.obs_mode == 'mIFU':
            origin_x = 0.81
            origin_y = 0.07

        lenght = 0.05
        width = 0.001
        head_width = 3 * width
        head_length = 1.5 * head_width
        text_factor = 1.2

        north_arrow = FancyArrow(
            origin_x, origin_y, 0, lenght,
            width=width, head_width=head_width, head_length=head_length,
            transform=fig.transFigure, figure=fig)

        east_arrow = FancyArrow(
            origin_x, origin_y, -self.figaspect * lenght, 0,
            width=(width / self.figaspect),
            head_width=(head_width / self.figaspect),
            head_length=(head_length * self.figaspect),
            transform=fig.transFigure, figure=fig)

        fig.lines.extend([north_arrow, east_arrow])

        fig.text(origin_x, origin_y + text_factor * lenght,
                 'N', ha='right')
        fig.text(origin_x - text_factor * self.figaspect * lenght, origin_y,
                 'E', ha='right')

    def _get_setted_fig_ax_vect(self):

        fig, ax_vect = self._get_fig_ax_vect()
        
        self._set_aspect_on_axes(ax_vect)
        self._set_ax_lims(ax_vect)
        self._clean_ticks_on_axes(ax_vect)
        self._add_north_east(fig)

        return fig, ax_vect

    def _set_title_on_fig(self, fig, suffix_title='', x=0.425, y=0.91):
        title = os.path.basename(self.xml_file) + suffix_title
        fig.text(x, y, title, ha='center', va='bottom')

    def _set_title_on_axes(self, ax_vect):
        raise NotImplementedError

    def _set_fig_legend(self, fig, text_list, title=True,
                        x=0.82, y=0.87, ha='left', va='top'):

        for i, element in enumerate(text_list):

            if (i == 0) or (title is False):
                fontsize = 10
            else:
                fontsize = 8

            if type(element) == str:
                text_line = element
                color = 'k'
            else:
                text_line, color = element
            text = i * '\n' + text_line

            fig.text(x, y, text, color=color, ha=ha, va=va, fontsize=fontsize)

    def _set_fig_colorbar(self, fig, mappable, left=0.85, height=0.8,
                          fraction=1):

        aux_ax = fig.add_axes([left, 0.07, 0.05, height])

        fig.colorbar(mappable, ax=aux_ax, fraction=fraction)

        aux_ax.remove()

    def _plot_weight_map(self):

        fig, ax_vect = self._get_setted_fig_ax_vect()

        self._set_title_on_fig(fig, suffix_title=' - Weight map')
        self._set_title_on_axes(ax_vect)

        if self.obs_mode == 'LIFU':
            # pix_size_arcsec = 0.50
            pix_size_arcsec = 0.10
        elif self.obs_mode == 'mIFU':
            # pix_size_arcsec = 0.25
            pix_size_arcsec = 0.05

        cmap = cm.gray_r
        vmin = 0
        vmax = len(self.data_dict_list)

        for i, ax in enumerate(ax_vect):

            num_samples = int(self.fov_list[i] * 3600 / pix_size_arcsec)

            ra_left, ra_right = ax_vect[i].get_xlim()
            ra_left_edge_vect, ra_step = np.linspace(
                ra_left, ra_right, num_samples, endpoint=False, retstep=True)
            ra_vect = ra_left_edge_vect + ra_step / 2

            dec_left, dec_right = ax_vect[i].get_ylim()
            dec_left_edge_vect, dec_step = np.linspace(
                dec_left, dec_right, num_samples, endpoint=False, retstep=True)
            dec_vect = dec_left_edge_vect + dec_step / 2

            ra_array, dec_array = np.meshgrid(ra_vect, dec_vect)

            array = np.zeros(ra_array.shape)

            for j, data_dict in enumerate(self.data_dict_list):

                for k, (ra, dec) in enumerate(zip(data_dict['targra'],
                                                  data_dict['targdec'])):

                    if ((ra >= ra_right) and (ra <= ra_left) and
                        (dec >= dec_left) and (dec <= dec_right)):

                        dist_array = np.hypot(
                            np.cos(np.deg2rad(self.field_dec)) * (ra_array - ra),
                            dec_array - dec)

                        mask = (dist_array <= self.fibre_size / 2)

                        array[mask] += 1

            extent = [ra_left, ra_right, dec_left, dec_right]

            ax_vect[i].imshow(array, cmap=cmap, vmin=vmin, vmax=vmax,
                              extent=extent, origin='lower')

        mappable = cm.ScalarMappable(norm=colors.Normalize(vmin=vmin,
                                                           vmax=vmax),
                                     cmap=cmap)
        self._set_fig_colorbar(fig, mappable)

        self.pdf.savefig(fig)
        plt.close(fig)

    def _get_color_list(self):

        color_list = [
            'red', 'lime', 'blue', 'gold', 'peru',
            'magenta', 'olive', 'cyan', 'indigo', 'silver',
            colors.TABLEAU_COLORS['tab:red'],
              colors.TABLEAU_COLORS['tab:green'],
              colors.TABLEAU_COLORS['tab:blue'],
              colors.TABLEAU_COLORS['tab:orange'],
              colors.TABLEAU_COLORS['tab:brown'],
            colors.TABLEAU_COLORS['tab:pink'],
              colors.TABLEAU_COLORS['tab:olive'],
              colors.TABLEAU_COLORS['tab:cyan'],
              colors.TABLEAU_COLORS['tab:purple'],
              colors.TABLEAU_COLORS['tab:gray']]

        return color_list

    def _assign_color_pointings(self, i, legend=False):

        color_list = self._get_color_list()

        if legend is False:
            facecolor = color_list[i]
            edgecolor = None
            result = (facecolor, edgecolor)
        else:
            legend_text_list = [
                ('Field #{}'.format(j + 1),
                 color_list[j % (len(color_list) - 1)])
                for j in range(i)]

            result = legend_text_list

        return result

    def _assign_color_attrib(self, type_attrib, attrib=None, set_of_values=None,
                             legend=False):

        color_list = self._get_color_list()

        if type_attrib in ['cname', 'ifu_spaxel', 'sim_template']:

            regex_dict = OrderedDict()

            if legend is False:

                if attrib is None:
                    result = ('w', 'r')
                elif type(attrib) is not str:
                    assert np.isnan(attrib)
                    result = ('w', 'k')
                else:
                    result = ('w', 'b')

                    if type_attrib == 'cname':
                        if re.match('^WVE_[0-9]{8}[+-][0-9]{7}$', attrib):
                            result = ('g', None)
                    elif type_attrib == 'ifu_spaxel':
                        if len(attrib) == 3:
                            if attrib[0] == 'S':
                                result = ('lime', 'gray')
                            else:
                                result = ('blue', 'gray')
                        elif len(attrib) == 6:
                            for i in range(20):
                                if attrib[:3] == 'm{:02d}'.format(i + 1):
                                    result = (color_list[i], None)
                                    break
                    elif type_attrib == 'sim_template':
                        result = ('g', None)

            else:

                legend_text_list = [type_attrib, '']

                if type_attrib == 'cname':
                    legend_text_list.append(('WVE_[0-9]+[+-][0-9]+', 'g'))
                elif type_attrib == 'ifu_spaxel':
                    legend_text_list.append(('S??', 'lime'))
                    legend_text_list.append(('???', 'blue'))
                    for i in range(20):
                        legend_text_list.append(
                            ('m{:02d}???'.format(i + 1), color_list[i]))
                elif type_attrib == 'sim_template':
                    legend_text_list.append(('Not empty', 'g'))

                legend_text_list.append(('-Unexpected non-empty-', 'b'))
                legend_text_list.append(('-Empty-', 'k'))
                legend_text_list.append(('-Missing-', 'r'))

                result = legend_text_list

            return result

        color_dict = OrderedDict()
        
        if type_attrib == 'targuse':
            color_dict['T'] = ('blue', None)
            color_dict['S'] = ('green', None)
            color_dict['C'] = ('red', None)
            color_dict['R'] = ('orange', None)
        elif type_attrib == 'targclass':
            color_dict[''] = ('w', 'k')

            targclass_list = [
                'UNKNOWN', 'SKY', 'GALAXY', 'NEBULA', 'QSO', 'STAR', 'STAR_BHB',
                'STAR_CEP', 'STAR_EM', 'STAR_EMP', 'STAR_FGK', 'STAR_IB',
                'STAR_MLT', 'STAR_MLUM', 'STAR_OB', 'STAR_BA', 'STAR_RRL',
                'STAR_VAR', 'STAR_WD', 'STAR_YSO']

            for i, targclass in enumerate(targclass_list):
                color_dict[targclass] = (color_list[i], None)
        else:

            if '' in set_of_values:
                color_dict[''] = ('w', 'k')

            if None in set_of_values:
                color_dict[None] = ('w', 'r')

            copy_set_of_values = set_of_values[:]
            if '' in copy_set_of_values:
                copy_set_of_values.remove('')
            if None in copy_set_of_values:
                copy_set_of_values.remove(None)

            for i, value in enumerate(copy_set_of_values):
                color_dict[value] = (color_list[i % len(color_list)], None)

        if legend is False:
            facecolor, edgecolor = color_dict[attrib]
            result = (facecolor, edgecolor)
        else:
            legend_text_list = [type_attrib, '']

            if '' in color_dict.keys():
                color_dict['-Empty-'] = color_dict[''][1]
                color_dict.pop('')

            if None in color_dict.keys():
                color_dict['-Missing-'] = color_dict[None][1]
                color_dict.pop(None)

            for key in color_dict.keys():
                legend_text_list.append((key, color_dict[key][0]))

            result = legend_text_list

        return result

    def _plot_on_ax_vect(self, ax_vect, contents='pointings_color', field=None):

        width = self.fibre_size
        height = self.fibre_size * np.cos(np.deg2rad(self.field_dec))

        if contents == 'pointings_color':
            set_of_values = None
            legend_text_list = self._assign_color_pointings(
                len(self.data_dict_list), legend=True)
            mappable = None
            alpha = 1 / len(self.data_dict_list)
        elif contents == 'pointings_bw':
            set_of_values = None
            legend_text_list = None
            mappable = None
            alpha = 1 / len(self.data_dict_list)
            facecolor = 'k'
            edgecolor = None
        elif contents in ['targuse', 'targclass', 'cname', 'ifu_spaxel',
                          'sim_template']:
            set_of_values = None
            legend_text_list = self._assign_color_attrib(contents, legend=True)
            mappable = None
            alpha = 1
        elif contents in ['targcat', 'targid', 'targname', 'targprog',
                          'targsrvy', 'sim_filterid']:
            values_array = np.concatenate(
                [data_dict[contents] for data_dict in self.data_dict_list])
            set_of_values = list(set(values_array))

            if None in set_of_values:
                set_of_values.remove(None)
                set_of_values.sort()
                set_of_values.append(None)
            else:
                set_of_values.sort()

            legend_text_list = self._assign_color_attrib(
                contents, set_of_values=set_of_values, legend=True)
            mappable = None
            alpha = 1
        else:
            set_of_values = None
            legend_text_list = None

            values_array = np.concatenate(
                [data_dict[contents] for data_dict in self.data_dict_list])

            clean_values_array = np.array(
                [value for value in values_array if value is not None])

            if np.all(np.isnan(clean_values_array)):
                vmin = -1
                vmax =  1
            else:
                vmin = np.nanmin(clean_values_array)
                vmax = np.nanmax(clean_values_array)

                if vmin == vmax:
                    vmin -= 0.5
                    vmax += 0.5

            mappable = cm.ScalarMappable(norm=colors.Normalize(vmin=vmin,
                                                               vmax=vmax),
                                         cmap=cm.rainbow_r)
            alpha = 1

        for i, ax in enumerate(ax_vect):

            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

            for j, data_dict in enumerate(self.data_dict_list):

                # Skip this field if a specific one has been requested

                if field is not None:
                    if j != field:
                        continue

                # Set the color if the contents are pointings
                
                if contents == 'pointings_color':
                    facecolor, edgecolor = self._assign_color_pointings(j)

                # Plot each fibre in the field of view of the axis

                for k, (ra, dec) in enumerate(zip(data_dict['targra'],
                                                  data_dict['targdec'])):

                    if ((ra >= xlim[1]) and (ra <= xlim[0]) and
                        (dec >= ylim[0]) and (dec <= ylim[1])):

                        if contents in ['targuse', 'targclass', 'targcat',
                                        'targid', 'targname', 'targprog',
                                        'targsrvy', 'cname', 'ifu_spaxel',
                                        'sim_filterid', 'sim_template']:
                            attrib = data_dict[contents][k]
                            facecolor, edgecolor = self._assign_color_attrib(
                                contents, attrib=attrib,
                                set_of_values=set_of_values)
                        elif not contents.startswith('pointings'):
                            attrib = data_dict[contents][k]

                            if attrib is None:
                                facecolor = 'w'
                                edgecolor = 'r'
                            elif np.isnan(attrib):
                                facecolor = 'w'
                                edgecolor = 'k'
                            else:
                                facecolor = mappable.to_rgba(attrib)
                                edgecolor = None

                        ax.add_artist(
                            Ellipse((ra, dec), width=width, height=height,
                                    facecolor=facecolor, edgecolor=edgecolor,
                                    alpha=alpha))

        return legend_text_list, mappable

    def _plot_pointings(self, color=True):

        fig, ax_vect = self._get_setted_fig_ax_vect()

        self._set_title_on_fig(fig)
        self._set_title_on_axes(ax_vect)

        if color is True:
            contents = 'pointings_color'
        else:
            contents = 'pointings_bw'

        legend_text_list, mappable = self._plot_on_ax_vect(ax_vect,
                                                           contents=contents)

        if legend_text_list is not None:
            self._set_fig_legend(fig, legend_text_list, title=False)

        self.pdf.savefig(fig)
        plt.close(fig)

    def _plot_attrib(self, type_attrib, split=True):

        if split is True:

            for i in range(len(self.data_dict_list)):

                fig, ax_vect = self._get_setted_fig_ax_vect()

                suffix_title = ' - {} - Field #{}'.format(type_attrib, i + 1)

                self._set_title_on_fig(fig, suffix_title=suffix_title)
                self._set_title_on_axes(ax_vect)

                legend_text_list, mappable = self._plot_on_ax_vect(
                    ax_vect, field=i, contents=type_attrib)

                if legend_text_list is not None:
                    self._set_fig_legend(fig, legend_text_list)

                if mappable is not None:
                    self._set_fig_colorbar(fig, mappable)

                self.pdf.savefig(fig)
                plt.close(fig)

        else:

            fig, ax_vect = self._get_setted_fig_ax_vect()

            suffix_title = ' - {}'.format(type_attrib)

            self._set_title_on_fig(fig, suffix_title=suffix_title)
            self._set_title_on_axes(ax_vect)

            legend_text_list, mappable = self._plot_on_ax_vect(
                ax_vect, contents=type_attrib)

            if legend_text_list is not None:
                self._set_fig_legend(fig, legend_text_list)

            if mappable is not None:
                self._set_fig_colorbar(fig, mappable)

            self.pdf.savefig(fig)
            plt.close(fig)

    def _make_plots(self):

        logging.info('starting plots for file {}'.format(
            self.pdf_filename))

        logging.info('\tplotting pointings of {} with color'.format(
            self.xml_file))
        self._plot_pointings(color=True)

        logging.info('\tplotting pointings of {} in black and white'.format(
            self.xml_file))
        self._plot_pointings(color=False)

        if self.verbose is True:
            logging.info('\tplotting weight map of {}'.format(self.xml_file))
            self._plot_weight_map()

        mag_list = ['mag_g', 'mag_r', 'mag_i',
                    'mag_gg', 'mag_bp', 'mag_rp']

        err_mag_list = ['e' + mag for mag in mag_list]

        non_verbose_attrib_list = ['targuse', 'targclass'] + mag_list

        sim_attrib_list = [
             'sim_filterid', 'sim_fwhm', 'sim_mag', 'sim_redshift',
             'sim_template', 'sim_velocity']

        verbose_attrib_list = (err_mag_list +
            ['targparal', 'targpmra', 'targpmdec', 'targepoch',
             'targcat', 'targsrvy', 'targprog', 'targname', 'targid',
             'ifu_pa', 'targprio', 'cname', 'ifu_spaxel',
             'automatic', 'configid', 'fibreid', 'targx', 'targy'])
        
        plot_attrib_list = non_verbose_attrib_list
        
        if self.verbose is True:
            plot_attrib_list.extend(verbose_attrib_list)
        
        if self.sim is True:
            plot_attrib_list.extend(sim_attrib_list)

        for attrib in plot_attrib_list:
            logging.info('\tplotting {} of {}'.format(attrib, self.xml_file))
            self._plot_attrib(attrib, split=False)
            self._plot_attrib(attrib, split=True)


class _LIFUPlot(_IFUPlot):

    def _set_fibre_size(self):
        self.fibre_size = 2.6 / 3600

    def _set_central_spaxel_list(self):
        self.central_spaxel_list = [
            'C14', 'S46', 'S39', 'S32', 'S25', 'S53', 'S04', 'S11', 'S18']

    def _set_fov_list(self):
        self.fov_list = [110/3600] + [25/3600 for i in range(8)]

    def _get_fig_ax_vect(self, nrows=4, ncols=6,
                         left=0.05, right=0.80, bottom=0.07, top=0.87):

        # Get the size of the fibure and create it

        figsize = plt.figaspect(self.figaspect)

        fig = plt.figure(figsize=figsize)

        # Get a grid for the axes

        gridspec = fig.add_gridspec(
            nrows, ncols, left=left, bottom=bottom, right=right, top=top)

        # Create the axes of the central bundle
        
        central_ax = fig.add_subplot(gridspec[:, 1:(ncols - 1)])

        # Create the axes of the sky bundles

        sky_ax_vect = []

        for i in range(4):
            sky_ax_vect.append(fig.add_subplot(gridspec[i, 0]))

        for i in range(4):
            sky_ax_vect.append(fig.add_subplot(gridspec[i, ncols - 1]))

        # Set some properties of the axes

        ax_vect = [central_ax] + sky_ax_vect

        return fig, ax_vect

    def _set_title_on_axes(self, ax_vect):

        ax_title_list = ['Ul', 'uL', 'dL', 'Dl', 'Ur', 'uR', 'dR', 'Dr']

        for i, ax_title in enumerate(ax_title_list):
            ax_vect[i + 1].set_title(ax_title, pad=0, fontdict={'fontsize': 10})

        data_dict = self.data_dict_list[0]

        for i, central_spaxel in enumerate(self.central_spaxel_list):

            index = data_dict['ifu_spaxel'].index(central_spaxel)
            ra = data_dict['targra'][index]
            dec = data_dict['targdec'][index]

            if i == 0:
                centre_coord = coordinates.SkyCoord(ra, dec, unit='deg')
            else:
                sky_bundle_coord = coordinates.SkyCoord(ra, dec, unit='deg')

                dist = centre_coord.separation(sky_bundle_coord)
                pos_angle = centre_coord.position_angle(sky_bundle_coord)

                text_coord = centre_coord.directional_offset_by(pos_angle,
                                                                dist / 3)

                ax_vect[0].text(text_coord.ra.deg, text_coord.dec.deg,
                                ax_title_list[i - 1])


class _MIFUPlot(_IFUPlot):

    def _set_fibre_size(self):
        self.fibre_size = 1.3 / 3600

    def _set_central_spaxel_list(self):
        self.central_spaxel_list = ['m{:02d}C04'.format(i + 1)
                                    for i in range(20)]

    def _set_fov_list(self):
        self.fov_list = [15/3600 for i in range(20)]

    def _get_fig_ax_vect(self, nrows=4, ncols=5,
                         left=0.10, right=0.75, bottom=0.07, top=0.87):

        # Get the size of the fibure and create it

        figsize = plt.figaspect(self.figaspect)

        fig = plt.figure(figsize=figsize)
    
        # Create a figure and the axes

        gridspec_kw = {'left': left, 'right': right,
                       'bottom': bottom, 'top': top}

        fig, ax_array = plt.subplots(
            figsize=figsize, nrows=nrows, ncols=ncols, gridspec_kw=gridspec_kw)

        ax_vect = ax_array.flatten()

        return fig, ax_vect

    def _set_title_on_axes(self, ax_vect):
        for i, ax in enumerate(ax_vect):
            ax.set_title('m{:02d}'.format(i + 1), pad=0,
                         fontdict={'fontsize': 10})


def plot_data_from_xml(xml_file, output_dir='output', sim=False, verbose=False):
    """
    Plot the data of a XML files.

    Parameters
    ----------
    xml_file : str
        The name of a XML file.
    output_dir : str, optional
        The directory which will contain the plot.
    sim : bool, optional
        Option for plotting simulation attributes
    verbose : bool, optional
        Option for verbose mode.

    Returns
    -------
    output_filename : str
        The filename of the output plot.
    """

    # Get rid of a matplotlib warning
    
    plt.rcParams.update({'figure.max_open_warning': 0})
    
    # Get the obs_mode of the XML
    
    obs_mode = get_obs_mode(xml_file)

    # Use a plotting class depending on the obs_mode

    if obs_mode == 'LIFU':

        with _LIFUPlot(xml_file, output_dir=output_dir, sim=sim,
                       verbose=verbose) as plot:
            output_filename = plot.get_pdf_filename()
        
    elif obs_mode == 'mIFU':

        with _MIFUPlot(xml_file, output_dir=output_dir, sim=sim,
                       verbose=verbose) as plot:
            output_filename = plot.get_pdf_filename()

    else:

        logging.error(
            'unexpected obs_mode in {}: {}'.format(xml_file, obs_mode))
    
        output_filename = None
        
    return output_filename


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Plot the data of one or more XML files.')

    parser.add_argument('xml_file', nargs='+', help='name of a XML file')

    parser.add_argument('--dir', default='output',
                        help='the directory which will contain the plots')

    parser.add_argument('--sim', action='store_true',
                        help='option for plotting simulation attributes')

    parser.add_argument('--verbose', action='store_true',
                        help='option for verbose mode')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()
    
    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    if not os.path.exists(args.dir):
        logging.info('Creating the output directory')
        os.mkdir(args.dir)
    
    for xml_file in args.xml_file:
        plot_data_from_xml(xml_file, output_dir=args.dir, sim=args.sim,
                           verbose=args.verbose)

