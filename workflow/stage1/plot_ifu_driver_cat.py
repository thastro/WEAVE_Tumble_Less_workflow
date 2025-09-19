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
import os.path
import logging
import subprocess

from astropy.io import fits

from workflow.utils.get_resources import get_aladin_jar
from workflow.utils.get_progtemp_info import (get_progtemp_dict,
                                              get_obsmode_from_progtemp)


def _send_command(cmd, p, option='write', encoding='utf-8'):

    assert option in ['write', 'communicate']

    if option == 'write':
        send_func = p.stdin.write
        send_func_str = 'write'
    elif option == 'communicate':
        send_func = p.communicate
        send_func_str = 'communicate'

    logging.debug('{}({})'.format(send_func_str, cmd))

    send_func('{}\n'.format(cmd).encode(encoding))


def plot_ifu_driver_cat(cat_filename, output_dir='output',
                        aladin_jar='Aladin.jar', progtemp_file=None,
                        pass_datamver=False):
    """
    Plot targets contained in an IFU driver catalogue with Aladin.

    Parameters
    ----------
    cat_filename : str
        A FITS file with an IFU driver catalogue.
    output_dir : str, optional
        The directory which will contain the plots generated with Aladin.
    aladin_jar : str, optional
        The location of the Java JAR file of Aladin.
    progtemp_file : str, optional
        A progtemp.dat file with the definition of PROGTEMP. If None, it will be
        retrieved from Internet.
    pass_datamver : bool, optional
        Continue even if DATAMVER mismatch is detected.
    """

    logging.info(
        'This plotting tool needs time. If you want to monitor its progress, '
        'see the directory which will contain its ouput images.')

    # Get the dictionary to interpret PROGTEMP and its DATAMVER value

    progtemp_datamver, progtemp_dict, forbidden_dict = get_progtemp_dict(
        filename=progtemp_file, assert_orb=True)

    # Open Aladin
    
    aladin_starting_cmd = 'java -jar {} -nogui'.format(aladin_jar)

    logging.debug('Popen({})'.format(aladin_starting_cmd))
    p = subprocess.Popen([aladin_starting_cmd], shell=True,
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)

    # Set the grid on

    _send_command('grid on', p)

    # Open the catalogue

    with fits.open(cat_filename) as hdu_list:
    
        # Check DATAMVER of the IFU driver cat and PROGTEMP file are consistent
        
        datamver = hdu_list[0].header['DATAMVER']

        if datamver != progtemp_datamver:
            logging.critical(
                'DATAMVER mismatch ({} != {}) for PROGTEMP file: '.format(
                    datamver, progtemp_datamver) +
                'Stop unless you are sure!')
            
            if pass_datamver == False:
                raise SystemExit(2)

        # Get the data from the catalogue

        data = hdu_list[1].data

        # Set the format string to build the numbers of the images
        
        num_digits = len(str(len(data)))
        num_digits_fmt = '{:0' + str(num_digits) + '}'

        # For each target in the catalogue

        for i in range(len(data)):

            # Get the needed information of the target

            ra = data['GAIA_RA'][i]
            dec = data['GAIA_DEC'][i]
            progtemp = data['PROGTEMP'][i]
            targname = data['TARGNAME'][i]
            targid = data['TARGID'][i]

            # Set a filename for the output image

            cat_basename_wo_ext = os.path.splitext(
                                      os.path.basename(cat_filename))[0]

            img_filename = (cat_basename_wo_ext + '-' +
                            num_digits_fmt.format(i + 1) + '-' + targname +
                            '-' + targid + '.png')
            img_path = os.path.join(output_dir, img_filename)

            # Create a the string with the coordinates which will be used in the
            # Aladin commands

            coord_str = '{:.10f} {:.10f}'.format(ra, dec)

            # Guess the obsmode and set some parameters acording to it

            obsmode = get_obsmode_from_progtemp(progtemp,
                                                progtemp_dict=progtemp_dict)

            if obsmode == 'LIFU':
                zoom_size_str = '10arcmin'
                get_radius_str = zoom_size_str
            elif obsmode == 'mIFU':
                zoom_size_str = '30arcsec'
                get_radius_str = zoom_size_str
            else:
                logging.warning(
                    'Skipping plotting of target {} of {}'.format(
                        i + 1, cat_filename))

                continue

            # Create a list of commands to be submitted to Aladin

            cmd_list = []

            # Reset

            cmd_list.append('reset')

            # Load the catalogue

            cmd_list.append('load {}'.format(os.path.abspath(cat_filename)))

            # Set the coordinates

            cmd_list.append(coord_str)

            # Get the image

            cmd_list.append('get hips(CDS/P/DSS2/color) {} {}'.format(coord_str,
                                                                get_radius_str))


            # Draw some circles to show the field of view of the instrument:
            # - For LIFU, draw circles to show the central bundle and the sky
            #   bundles.
            # - For mIFU, draw a circle for each bundle

            if obsmode == 'LIFU':
                cmd_list.append('draw yellow circle({} 45.097arcsec)'.format(
                    coord_str))
                cmd_list.append('draw green circle({} 249.203arcsec)'.format(
                    coord_str))
                cmd_list.append('draw green circle({} 258.169arcsec)'.format(
                    coord_str))
            elif obsmode == 'mIFU':
                cmd_list.append('draw red circle({} 6.088arcsec)'.format(
                    coord_str))

            # Set the desired zoom

            cmd_list.append('zoom {}'.format(zoom_size_str))

            # Save the image

            cmd_list.append('save {}'.format(os.path.abspath(img_path)))

            # Repeat some commands from above due to a bug in Aladin

            cmd_list.append(coord_str)
            cmd_list.append('zoom {}'.format(zoom_size_str))
            cmd_list.append('save {}'.format(os.path.abspath(img_path)))

            # Join the commands for this target and send them to Aladin

            cmd = '; '.join(cmd_list)
            _send_command(cmd, p)

        # Ask Aladin to quit and wait

        _send_command('quit', p, option='communicate')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=
        """
        Plot targets contained in an IFU driver catalogue with Aladin.
        """)

    parser.add_argument('catalogue',
                        help='a FITS file with an IFU driver catalogue')

    parser.add_argument('--outdir', default='output', help=
                        """
                        the directory which will contain the plots generated
                        with Aladin
                        """)

    parser.add_argument('--aladin', default=os.path.join('aux', 'Aladin.jar'),
                        help='the location of the Java JAR file of Aladin')

    parser.add_argument('--progtemp_file', default=None,
                        help="""a progtemp.dat file with the definition of
                        PROGTEMP. If it is not provided, the latest version will
                        be retrieved from Internet""")

    parser.add_argument('--pass_datamver', action='store_true',
                        help='continue even if DATAMVER mismatch is detected')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()
    
    logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    if not os.path.exists(args.aladin):
        logging.info('Downloading the Java JAR file of Aladin')
        get_aladin_jar(file_path=args.aladin)
    
    plot_ifu_driver_cat(args.catalogue, output_dir=args.outdir,
                        aladin_jar=args.aladin,
                        progtemp_file=args.progtemp_file,
                        pass_datamver=args.pass_datamver)

