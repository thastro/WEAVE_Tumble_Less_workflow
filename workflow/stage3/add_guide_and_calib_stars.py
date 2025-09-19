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
import os

from workflow.utils.classes import OBXML


def add_guide_and_calib_stars(
        xml_file_list, output_dir, lifu_num_guide_stars_request=1,
        mifu_num_guide_stars_request=None, mifu_num_central_guide_stars=1,
        mifu_min_guide_cut=0.3, mifu_max_guide_cut=0.9,
        num_calib_stars_request=None, num_central_calib_stars=0,
        min_calib_cut=0.0, max_calib_cut=1.0, write_useful_tables=False,
        cnames_for_guide_stars=True, cnames_for_calib_stars=True,
        overwrite=False, cache=False, guide_mag_lim=99., calib_mag_lim=99.,
        ha_range=2.0, calib_mixres=False):
    """
    Add guide and calib stars to XML files.
    
    Parameters
    ----------
    xml_file_list : list of str
        A list of input OB XML files.
    output_dir : str
        Name of the directory which will containe the output XML files.
    lifu_num_guide_stars_request : int, optional
        Maximum number of LIFU guide stars in the output. None means no limit.
    mifu_num_guide_stars_request : int, optional
        Maximum number of mIFU guide stars in the output. None means no limit.
    mifu_num_central_guide_stars : int, optional
        Number of mIFU guide stars near to centre to be selected.
    mifu_min_guide_cut : float, optional
        Minimum cut factor to be used for the non-central mIFU guide stars.
    mifu_max_guide_cut : float, optional
        Maximum cut factor to be used for the non-central mIFU guide stars.
    num_calib_stars_request : int, optional
        Maximum number of calib stars in the output. None means no limit.
    num_central_calib_stars : int, optional
        Number of calib stars near to centre to be selected.
    min_calib_cut : float, optional
        Minimum cut factor to be used for the non-central calib stars.
    max_calib_cut : float, optional
        Maximum cut factor to be used for the non-central calib stars.
    write_useful_tables : bool, optional
        Write tables with the potentially useful stars, i.e. those guide or
        calib stars available for the position of the location in the sky of
        each OB XML file.
    cnames_for_guide_stars : bool, optional
        Add or remove CNAMEs for guide stars in the output XMLs.
    cnames_for_calib_stars : bool, optional
        Add or remove CNAMEs for calib stars in the output XMLs.
    overwrite : bool, optional
        Overwrite the output FITS file.
    cache : bool, optional
            Do not download catalogue files if they aready exist
    guide_mag_lim : float, optional
            Apply a faint limit cut in Gaia G to the available guide stars        
    calib_mag_lim : float, optional
            Apply a faint limit cut in Gaia G to the available calibration stars        
    ha_range : float, optional
            Define an hour-angle range (+/-ha) required for this observation
            (relevant for LIFU, impact on observing within rotator limits)
    calib_mixres : bool, optional
            Permit selection of LR (fainter) calibration targets under a HR field
            where no (or fewer) HR calibration targets are available.
    


    Returns
    -------
    output_file_list : list of str
        A list with the output XML files.
    """

    output_file_list = []

    xml_file_list.sort()

    for xml_file in xml_file_list:

        # Check that the input XML exists and is a file

        assert os.path.isfile(xml_file)

        # Choose the output filename depedending on the input filename

        input_basename_wo_ext = os.path.splitext(os.path.basename(xml_file))[0]

        if (input_basename_wo_ext.endswith('-t') or
            input_basename_wo_ext.endswith('-')):
            output_basename_wo_ext = input_basename_wo_ext + 'gc'
        else:
            output_basename_wo_ext = input_basename_wo_ext + '-gc'

        output_file = os.path.join(output_dir, output_basename_wo_ext + '.xml')
        guide_plot_filename = os.path.join(
            output_dir, output_basename_wo_ext + '-guide_stars.png')
        calib_plot_filename = os.path.join(
            output_dir, output_basename_wo_ext + '-calib_stars.png')
        
        if write_useful_tables is True:
            guide_useful_table_filename = os.path.join(
                output_dir, output_basename_wo_ext + '-useful_guide_stars.fits')
            calib_useful_table_filename = os.path.join(
                output_dir, output_basename_wo_ext + '-useful_calib_stars.fits')
        else:
            guide_useful_table_filename = None
            calib_useful_table_filename = None

        # Save the output filename for the result

        output_file_list.append(output_file)

        # If the output file already exists, delete it or continue with the next
        # one

        if os.path.exists(output_file):
            if overwrite == True:
                logging.info('Removing previous file: {}'.format(output_file))
                os.remove(output_file)
                
                for filename in [guide_plot_filename, calib_plot_filename,
                                 guide_useful_table_filename,
                                 calib_useful_table_filename]:
                    if (filename is not None) and os.path.exists(filename):
                        logging.info(
                            'Removing previous file: {}'.format(filename))
                        os.remove(filename)
                        
            else:
                logging.info(
                    'Skipping file {} as its output already exists: {}'.format(
                        xml_file, output_file))
                continue

        # Read the input file, add the guide and calib stars and write it to the
        # output file
        
        ob_xml = OBXML(xml_file)

        ob_xml.add_guide_and_calib_stars(
            guide_plot_filename=guide_plot_filename,
            guide_useful_table_filename=guide_useful_table_filename,
            lifu_num_guide_stars_request=lifu_num_guide_stars_request,
            mifu_num_guide_stars_request=mifu_num_guide_stars_request,
            mifu_num_central_guide_stars=mifu_num_central_guide_stars,
            mifu_min_guide_cut=mifu_min_guide_cut,
            mifu_max_guide_cut=mifu_max_guide_cut,
            calib_plot_filename=calib_plot_filename,
            calib_useful_table_filename=calib_useful_table_filename,
            num_calib_stars_request=num_calib_stars_request,
            num_central_calib_stars=num_central_calib_stars,
            min_calib_cut=min_calib_cut, max_calib_cut=max_calib_cut,
            cnames_for_guide_stars=cnames_for_guide_stars,
            cnames_for_calib_stars=cnames_for_calib_stars,
            cache=cache, guide_mag_lim=guide_mag_lim,
            calib_mag_lim=calib_mag_lim, ha_range=ha_range,
            calib_mixres=calib_mixres)

        ob_xml.write_xml(output_file)

    return output_file_list


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
             description='Add guide and calib stars to XML files')

    parser.add_argument('xml_file', nargs='+',
                        help="""an input OB XML file""")

    parser.add_argument('--outdir', dest='output_dir', default='output',
                        help="""name of the directory which will contain the
                        output XML files""")

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

    parser.add_argument('--write_useful_tables', action='store_true',
                        help="""write tables with the potentially useful stars,
                        i.e. those guide or calib stars available for the
                        position of the location in the sky of each OB XML
                        file""")

    parser.add_argument('--add_cnames_for_guide_stars',
                        dest='cnames_for_guide_stars', action='store_true',
                        help='add CNAMEs for guide stars in the output XMLs')

    parser.add_argument('--remove_cnames_for_guide_stars',
                        dest='cnames_for_guide_stars', action='store_false',
                        help='remove CNAMEs for guide stars in the output XMLs')

    parser.set_defaults(cnames_for_guide_stars=True)

    parser.add_argument('--add_cnames_for_calib_stars',
                        dest='cnames_for_calib_stars', action='store_true',
                        help='add CNAMEs for calib stars in the output XMLs')

    parser.add_argument('--remove_cnames_for_calib_stars',
                        dest='cnames_for_calib_stars', action='store_false',
                        help='remove CNAMEs for calib stars in the output XMLs')

    parser.set_defaults(cnames_for_calib_stars=True)

    parser.add_argument('--overwrite', action='store_true',
                        help='overwrite the output files')

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
    
    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    if not os.path.exists(args.output_dir):
        logging.info('Creating the output directory')
        os.mkdir(args.output_dir)
    
    if args.lifu_num_guide_stars_request > -1:
        lifu_num_guide_stars_request = int(args.lifu_num_guide_stars_request)
    else:
        lifu_num_guide_stars_request = None
    
    if args.mifu_num_guide_stars_request > -1:
        mifu_num_guide_stars_request = int(args.mifu_num_guide_stars_request)
    else:
        mifu_num_guide_stars_request = None
    
    if args.num_calib_stars_request > -1:
        num_calib_stars_request = args.num_calib_stars_request
        assert type(num_calib_stars_request) == int
    else:
        num_calib_stars_request = None

    add_guide_and_calib_stars(
        args.xml_file, args.output_dir,
        lifu_num_guide_stars_request=lifu_num_guide_stars_request,
        mifu_num_guide_stars_request=mifu_num_guide_stars_request,
        mifu_num_central_guide_stars=args.mifu_num_central_guide_stars,
        mifu_min_guide_cut=args.mifu_min_guide_cut,
        mifu_max_guide_cut=args.mifu_max_guide_cut,
        num_calib_stars_request=num_calib_stars_request,
        num_central_calib_stars=args.num_central_calib_stars,
        min_calib_cut=args.min_calib_cut,
        max_calib_cut=args.max_calib_cut,
        write_useful_tables=args.write_useful_tables,
        cnames_for_guide_stars=args.cnames_for_guide_stars,
        cnames_for_calib_stars=args.cnames_for_calib_stars,
        overwrite=args.overwrite,
        cache=args.cache,
        guide_mag_lim=args.guide_maglim,
        calib_mag_lim=args.calib_maglim,
        calib_mixres=args.calib_mixres)

