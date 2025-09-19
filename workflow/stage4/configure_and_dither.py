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

import numpy
import glob
import logging
from subprocess import Popen, PIPE
import argparse
import os
import postproc
import platform
import preproc
import time
from xml.etree.ElementTree import ElementTree as ET
import xml.etree.cElementTree as ElementTree


do_preproc = False
do_postproc = False



def check_halims(xml_filename_list, ha_range=2.0):
    # Check the output HA limits to ensure still within rotator limits

    # As this is the first we've seen of the XMLs since configure, we need to check
    # their min/max ha limits do not exceed the canonical HA range (warn if it does)

    # xml_halims = get_halims_from_xmls(xml_filename_list)


    # Checks that the XMLs have ha ranges below the canonical used for rotator
    # limit checks. If not, warn!

    ha_mins = [_get_attribute_in_simple_element_of_xml_file(
                          xml_file, 'hour_angle_limits', 'earliest')
                      for xml_file in xml_filename_list]

    ha_maxs = [_get_attribute_in_simple_element_of_xml_file(
                          xml_file, 'hour_angle_limits', 'latest')
                      for xml_file in xml_filename_list]

    obs_modes = [_get_attribute_in_simple_element_of_xml_file(
                          xml_file, 'observation', 'obs_mode')
                      for xml_file in xml_filename_list]

    failure = False
    
    for ha_min,ha_max,obs_mode in zip(ha_mins,ha_maxs,obs_modes):
        if obs_mode != 'LIFU':
            # Don't care - Configure should have selected viable HA ranges
            continue
        xml_fail = False
        if ha_min < -1.0 * ha_range:
            failure = True
            xml_fail = True
        if ha_max > ha_range:
            failure = True
            xml_fail = True

        if xml_fail:
            logging.warning("XML file has HA limits [{:1.3f},{:1.3f}] that exceeds the canonical limits [{:1.3f},{:1.3f}]".format(
                ha_min,ha_max,-1.0*ha_range,ha_range))
            

    if failure:
        logging.warning("")
        logging.warning("One or more of the provided XML files have HA ranges outside of those used to")
        logging.warning("ensure pointings (Dec, PA, HA range) are within WEAVE rotator limits.")
        logging.warning("This means it is *possible* that these XMLs cannot now be observed: you may")
        logging.warning("proceed, but be aware that WASP might reject these submissions for the reason above")


def get_config_path(version=None):
    cfg_path = '/Applications/configure.app/Contents/MacOS/'
    assert platform.system() == 'Darwin', "Until futher notice, Linux systems *must* explicitly specify the configure install location with the --cfg_path"
        
    

    if version:
        cfg_path = cfg_path.replace('configure.app','configure.app.{}'.format(version))
        assert os.path.isdir(cfg_path), "Requested version ('{}') has no install directory at {}".format(version,cfg_path)


    assert os.path.isdir(cfg_path), "Configure is not installed at {} : please specify an install loc with --cfg_path".format(cfg_path)        
    return cfg_path 

def check_lint():
    cmd = "xmllint --version"
    cmd = cmd.split()
    proc = Popen(cmd,stdout=PIPE,stderr=PIPE)
    comm = proc.communicate()
    out = comm[0].decode('utf8')
    err = comm[1].decode('utf8')

    # it puts the output into stderr for some reason?
    # check both...
    assert ('compiled with:' in out) or ('compiled with:' in err), "Please install xmllint or specify --nolint to skip"

    return True

        
def run_configure(xmls, config_path=None, version=None,
                  dither=True, clobber=True, lint=True,
                  src_cfg=None, nogui=False,set_pa=None,
                  set_pressure=None, cfg_override=False,
                  wipe_on_fail=False):
    if (not config_path) or (version):
        config_path = get_config_path(version=version)

    logging.info("Using configure at: {}".format(config_path))

    if lint:
        check_lint()
    
    failed = []
    passed = []

    logging.warning("STDERR checking in place until diagnostic outputs (3.4.5) removed - revert once Configure writes messaging to STDOUT")

    
    for xml in xmls:
        field_outputs = []
        if not 'ifu_' in xml:
            logging.info("Skipping {}".format(xml))
            continue
        assert '-tgc.xml' in xml, "Expecting filename containing -tgc.xml: please follow workflow tutorial"

        try:
            ha_file = xml.replace('-tgc.xml','.ha')
            ha = float(open(ha_file,'r').read().split('\n')[0])
        except:
            ha = 0.0
        try:
            epoch_file = xml.replace('-tgc.xml','.epoch')
            epoch = float(open(epoch_file,'r').read().split('\n')[0])
        except:
            epoch = 2022.0

        cfg_output = xml.replace('-tgc.xml','-tgcs-tmp-orig.xml')
        if not dither:
            # the cfg output needs to be ready for the next stage
            cfg_output = xml.replace('-tgc.xml','-tgcs-orig.xml')
            if not lint:
                cfg_output = xml.replace('-tgc.xml','-tgcs.xml')

        if os.path.isfile(cfg_output):
            if not clobber:
                logging.error("Will not clobber existing file {}".format(cfg_output))
                failed.append(cfg_output)
                continue
            logging.warning("Removing existing file {}".format(cfg_output))
            os.unlink(cfg_output)            
            
        ifu_mode = 'lifu_configure'
        input_switch = 'i'
        gui_mode = ''
        config_file = ''
        if '-mifu_' in xml:
            ifu_mode = 'configure'
            input_switch = 'f'
            if nogui:
                gui_mode = ' --gui 0'
        if cfg_override:
            assert os.path.isfile(cfg_override), "Provided configure XML override file {} does not exist".format(cfg_override)
            config_file = ' --config_file {}'.format(cfg_override)

                
        if (version) and (do_preproc):
            # this is not the latest release - write the RA_d and Dec_d into the field from C14 
            # (will skip if this is not LIFU)
            preproc.horiz_fix(xml)

        if set_pa != None:
            preproc.set_pa(xml,set_pa)
            
        if (src_cfg) and ('-mifu_' in xml):
            if os.path.isdir(src_cfg):
                # this is a directory - find out if there is an analogous file to the one being worked on...
                src_file = os.path.basename(xml).replace('tgc.xml','tgcs-tmp-orig.xml')
                src_file = os.path.join(src_cfg,src_file)
                if os.path.isfile(src_file):
                    preproc.auto_mifu(xml,src_file)
                    if gui_mode == '':
                        #set no-gui for the user...
                        gui_mode = ' --gui 0'
                else:
                    logging.warning("No viable source XML for bundle/guide transplant for {}".format(os.path.basename(xml)))
            else:
                preproc.auto_mifu(xml,src_cfg)
                if gui_mode == '':
                    #set no-gui for the user...
                    gui_mode = ' --gui 0'
                

        pressure_str = ''
        if set_pressure:
            pressure_str = ' --pressure {:.1f}'.format(set_pressure)
            
        cmd = "{}{}{}{} --epoch {:.3f} --ha {:.3f}{} -{} {} -o {}".format(
            os.path.join(config_path,''),
            ifu_mode,gui_mode,config_file,
            epoch,ha,pressure_str,input_switch,
            xml,cfg_output)
            
        logging.info(cmd)
        cmd = cmd.split()
        
        if ifu_mode != 'lifu_configure':
            user_choice = None
            user_done = False
            while not user_done:

                if nogui:
                    user_done = True
                elif user_choice != '3':
                    logging.info("")
                    logging.info("--------- Interactive configure mode ---------")
                    logging.info("\tOnce you press enter, the Configure UI will launch for file:")
                    logging.info("\t"+xml)
                    logging.info("\tClick 'Open' in the UI to open this field.")
                    logging.info("\tYou must now assign mIFU bundles for this field, including guide stars")
                    logging.info("\tand required calibration stars. Please consult Configure documenation for")
                    logging.info("\tguidance on this.")

                    logging.info("")
                    logging.info("\tThe output file for this has been defined as:")
                    logging.info("\t"+cfg_output)
                    logging.info("\tWhen you have completed configuration of this field, click 'Save'.")
                    logging.info("\tThe output filename provided in the save dialogue is correct - please do")
                    logging.info("\tnot change it, otherwise downstream IFU workflow operations will fail.")
                    logging.info("\tOnce the file is saved, click 'Quit'.") 
                    logging.info("\tIf you wish to re-do this field, you can quit and choose to try again.")
                    input("Press Enter to continue...")
                    
                elif user_choice == '3':


                    # get a new input/output XML (the function will write it)
                    xml, cfg_output = filter_targets(xml,field_outputs)
                    
                    # now reset the command 
                    cmd = "{}{} --epoch {:.3f} --ha {:.3f} -{} {} -o {}".format(
                        os.path.join(config_path,''),
                        ifu_mode,
                        epoch,ha,input_switch,
                        xml,cfg_output)
                    logging.info(cmd)
                    cmd = cmd.split()
                    
                    logging.info("")
                    logging.info("--------- Interactive configure mode (filtered field) ---------")
                    logging.info("\tOnce you press enter, the Configure UI will launch for file:")
                    logging.info("\t"+xml)
                    logging.info("\tThis file has filtered all previous (T)argets that have been")
                    logging.info("\tassigned bundles under this requested field. You can continue the")
                    logging.info("\telimination process after this field has been configured and saved")
                    logging.info("")
                    logging.info("\tAfter the GUI launches, click 'Open' to load and display this field.")
                    logging.info("\tYou must then assign mIFU bundles for this field, including guide stars")
                    logging.info("\tand required calibration stars. Please consult Configure documenation for")
                    logging.info("\tguidance on this.")
                    input("Press Enter to continue...")                    


                    
                proc = Popen(cmd,stdout=PIPE,stderr=PIPE)
                comm = proc.communicate()
                logging.debug(comm)

                if (comm[1].decode('utf8') != ''):
                    is_warning = False
                    for _msg in comm[1].decode('utf8').split('\n'):
                        if 'Warning' in _msg:
                            is_warning = True

                    if not is_warning:
                        # failed.append(cfg_output)
                        for _msg in comm[1].decode('utf8').split('\n'):
                            if ('guide targets lost' in _msg):
                                logging.warning(_msg)
                                continue
                            elif _msg != '':
                                if not cfg_output in failed:
                                    failed.append(cfg_output)
                                    logging.error(_msg)
                                    continue
                    
                    else:
                        for _msg in comm[1].decode('utf8').split('\n'):
                            logging.warning(_msg)
                        
                field_outputs.append(cfg_output)

                
                if not nogui:
                    logging.info("")
                    logging.info("Do you want to: continue to the next field (1)")
                    logging.info("(1) Continue to the next field")
                    logging.info("(2) Redo this field")
                    logging.info("(3) Select more targets in this field with already-allocated targets filtered out")
                    logging.info("Please press 1, 2 or 3 then hit ENTER")
                    user_choice = None
                    while user_choice not in ["1","2","3"]:
                        user_choice = input("Enter option: ")
                    user_done = user_choice == "1"
                    
        else:
            proc = Popen(cmd,stdout=PIPE,stderr=PIPE)
            comm = proc.communicate()
            logging.debug(comm)

            if comm[1].decode('utf8') != '':
                is_warning = False
                for _msg in comm[1].decode('utf8').split('\n'):
                    if 'Warning' in _msg:
                        is_warning = True

                if not is_warning:
                    failed.append(cfg_output)
                    for _msg in comm[1].decode('utf8').split('\n'):
                        logging.error(_msg)
                    continue
                else:
                    for _msg in comm[1].decode('utf8').split('\n'):
                        logging.warning(_msg)

                
            field_outputs.append(cfg_output)

            
        if do_postproc:
            # add any configure postprocs here
            for cfg_output in field_outputs:
                postproc.filter_spares(cfg_output)
                postproc.populate_sky(cfg_output)
                postproc.guide_pa_unset(cfg_output)
                postproc.reset_random_cnames(cfg_output)

                
        if lint:
            _field_outputs = []
            for cfg_output in field_outputs:
                lint_output = cfg_output
                if not dither:
                    # need to drop the -orig and write out as if it was from dither stage
                    lint_output = cfg_output.replace('-tgcs-orig.xml','-tgcs.xml')
                cmd = "xmllint --c14n {} | xmllint --format --encode utf-8 --output {} -".format(cfg_output,lint_output)
                _field_outputs.append(lint_output)
                logging.debug(cmd)
                proc = Popen(cmd,shell=True)
            field_outputs = _field_outputs
            
        for cfg_output in field_outputs:
            if os.path.isfile(cfg_output):
                passed.append(cfg_output)
            else:
                time.sleep(0.2)
                if os.path.isfile(cfg_output):
                    passed.append(cfg_output)
                else:
                    failed.append(cfg_output)
                continue

    if len(failed):
        fail_state = wait_failures(failed)
        if not fail_state:
            passed = passed + failed
            failed = []
            
    if len(failed):
        logging.error("Successfully configured {} XMLs, but {} failed - bailing".format(len(passed),len(failed)))
        if wipe_on_fail:
            logging.warning("Removing all output products generated so far (--fail_wipe):")
            for _file in passed:
                logging.warning("\t{}".format(_file))
                os.unlink(_file)

        return []
    else:
        if len(passed) == 1:
            logging.info("Successfully configured {}".format(passed[0]))
        else:
            logging.info("Successfully configured {} XMLs".format(len(passed)))
    return passed

def filter_targets(input_xml, configured_xmls,outname=None,
                   guides=False, skies=False, calibs=False):
    """For a given input XML file that you wish to allocate mIFU fibres onto,
    and one or more configured XMLs, produce a new input XML that excludes all
    targets that have already been allocated a bundle
 
    if guides / sky /calibs == True: also filter out these target types
    (not advised)
    
    Need to return:
    - filename of input file (must exist - write it in this function)
    - proposed output fileame (should not exist!)


    """
    assert len(configured_xmls)
    assert os.path.isfile(input_xml)

    # get the source XML file loaded - ensure it's just one field
    # (otherwise things get complicated)

    lookup_attributes = {}
    lookup_attributes['G'] = 'cname'
    lookup_attributes['T'] = 'targid'
    lookup_attributes['C'] = 'cname'

    filter_targets = ['T']
    if guides:
        filter_targets.append('G')
    if skies:
        filter_targets.append('S')
    if calibs:
        filter_targets.append('C')
    
    stree = ElementTree.parse(input_xml)
    sroot = stree.getroot()
    sfield = sroot.findall('.//field')
    assert len(sfield) == 1
    sfield = sfield[0]
    stargets = sfield.findall('.//target')

    filter_lookup = {}

    ob_names = []
    # Now build the filter list
    for xml in configured_xmls:
        ftree = ElementTree.parse(xml)
        froot = ftree.getroot()
        ffield = froot.find('.//field')
        filter_targs = ffield.findall(""".//target""")
        fobservaion = froot.find('observation')
        fobname = fobservaion.get('name')
        ob_names.append(fobname)
        for target in filter_targs:
            if target.get('ifu_spaxel'):
                if not 'C04' in target.get('ifu_spaxel'):
                    continue
            if target.get('fibreid') and (target.get('fibreid') not in ['','%%%']):
                if target.get('targuse') in filter_targets:
                    # this is in the allowed types that can be filtered
                    filter_lookup[target.get(lookup_attributes[target.get('targuse')])] = 1

    logging.info("Identified {} possible targets to remove".format(len(filter_lookup)))
                
    # Now filter any matches in the input file
    nfiltered = 0
    for target in stargets:
        try:
            match = filter_lookup[target.get(lookup_attributes[target.get('targuse')])]
            sfield.remove(target)
            nfiltered += 1
        except KeyError:
            # is it nevertheless allocated a bundle (eg, a calib bundle)?
            # if so, park it so it can be chosen again, or not...
            if target.get('fibreid'):
                del target.attrib['fibreid']
            if target.get('configid'):
                del target.attrib['configid']

    observation = sroot.find('observation')
    obname = observation.get('name')
    if not outname:
        # choose a destination output name for this filter
        # take the last of the outputs, iterate until the file doesn't exist
        _outname = configured_xmls[-1]
        outname,index = construct_outname(_outname)
        index = index.split('-')[1]
        obname += "/{}".format(index)
    else:
        # don't know how to define the index in this case, so will have to hope
        # this is OK, but can check against what's come before and iterate until
        # it's not in the list
        index = '/sub1'
        _obname = obname+index
        while _obname in ob_names:
            index_no = int(_obname.split('/sub')[1])
            index_np += 1
            _obname = _obname.replace(index,'/sub{}'.format(index_np))
        obname = _obname

    # update the OB name in the XML:
    observation.set('name',obname)
        
    # set the input name (and write the tree to it)

    xml = outname.replace('-tgcs-tmp-orig.xml','-tgc.xml')
    xml_tree = ET(sroot)
    file = open(xml,'w')
    logging.info("Writing filtered configure input XML to {}".format(xml))
    xml_tree.write(file,encoding='unicode')
    file.close()

    return xml, outname
    
    
def construct_outname(outname):
    while(os.path.isfile(outname)):
          suffix = '-tgcs-tmp-orig.xml'
          prefix = outname.split(suffix)[0]
          _index = prefix.split('-mifu_')[1]
          index = generate_index(_index)
          # now put it all back together
          outname = "{}{}".format(prefix.replace(_index,index),suffix)
    return outname,index
        
def generate_index(index):
    if 'sub' in index:
        ii = index.split('sub')[1]
        index = index.replace("sub{}".format(ii),"sub{}".format(str(int(ii)+1)))
    else:
        # add the first sub
        index += '-sub1'
    return index
        
        
def wait_failures(xmls,wait=2,retries=10):
    logging.warning("{} XML files did not report as completed immediately - waiting on file system to catch up...".format(len(xmls)))
    for i in range(retries):
        outmode = logging.debug
        outmode("Retry {}".format(i+1))
        if not False in [os.path.isfile(xml) for xml in xmls]:
            break
        else:
            time.sleep(wait)
    if False in [os.path.isfile(xml) for xml in xmls]:
        return True
    return False

        
def run_dither(xmls, config_path=None, version=None,
               clobber=True, lint=True,
               cfg_override=None,
               wipe_on_fail=False):
    if (not config_path) or (version):
        config_path = get_config_path(version=version)

    time.sleep(1.0)
    logging.info("Using dither at: {}".format(config_path))

    if lint:
        check_lint()

    if not len(xmls):
        logging.error("No input to dither provided")
        return []
        
    failed = []
    passed = []

    for xml in xmls:
        if not 'ifu_' in xml:
            logging.info("Skipping {}".format(xml))
            continue
        assert '-tgcs-tmp-orig.xml' in xml, "Expecting filename containing -tgcs-tmp-orig.xml: please follow workflow tutorial"

        dth_output = xml.replace('-tgcs-tmp-orig.xml','-tgcs-orig.xml')
        if not lint:
            dth_output = xml.replace('-tgcs-tmp-orig.xml','-tgcs.xml')
        
        if os.path.isfile(dth_output):
            if not clobber:
                logging.error("Will not clobber existing file {}".format(dth_output))
                failed.append(dth_output)
                continue
            logging.warning("Removing existing file {}".format(dth_output))
            os.unlink(dth_output)

        config_file = ''
        if cfg_override:
            assert os.path.isfile(cfg_override), "Provided configure XML override file {} does not exist".format(cfg_override)
            config_file = ' --config_file {}'.format(cfg_override)
            
        cmd = "{}dither{} --input {} --output {}".format(
            os.path.join(config_path,''),
            config_file,
            xml,dth_output)
            
        logging.info(cmd)

        cmd = cmd.split()
        proc = Popen(cmd,stdout=PIPE,stderr=PIPE,shell=False)
        comm = proc.communicate()
        logging.debug(comm)

        if comm[1].decode('utf8') != '':
            failed.append(dth_output)
            for _msg in comm[1].decode('utf8').split('\n'):
                logging.error(_msg)
            continue

        if do_postproc:
            # add any dither postprocs here
            postproc.add_offsets(dth_output)
            postproc.populate_random_photometry(dth_output)
            postproc.reset_targpos(dth_output,xml)

            
        output = dth_output
        if lint:
            lint_output = dth_output.replace('-tgcs-orig.xml','-tgcs.xml')
            cmd = "xmllint --c14n {} | xmllint --format --encode utf-8 --output {} -".format(dth_output,lint_output)
            logging.debug(cmd)
            proc = Popen(cmd,shell=True)
            output = lint_output
            
        if os.path.isfile(output):
            passed.append(output)
        else:
            time.sleep(0.2)
            if os.path.isfile(output):
                passed.append(output)
            else:
                failed.append(output)
            continue

    if len(failed):
        fail_state = wait_failures(failed)
        if not fail_state:
            passed = passed + failed
            failed = []
            
    if len(failed):
        logging.error("Successfully dithered {} XMLs, but {} failed - bailing".format(len(passed),len(failed)))
        if wipe_on_fail:
            logging.warning("Removing all output products generated so far (--fail_wipe):")
            for _file in passed:
                logging.warning("\t{}".format(_file))
                os.unlink(_file)
        return []
    else:
        if len(passed) == 1:
            logging.info("Successfully dithered {}".format(passed[0]))
        else:
            logging.info("Successfully dithered {} XMLs".format(len(passed)))
        
    return passed


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run configure and dither within the workflow')

    parser.add_argument('xml_files', nargs='+',
                        help="""an input OB XML file""")

    parser.add_argument('--cfg_path', default=None,
                        help='Override the default path for configure/dither')

    parser.add_argument('--cfg_ver', default=None,
                        help='Define a specific Configure version to run (must be installed). Will override cfg_path')

    parser.add_argument('--cfg_override', default=None,
                        help='Define an alternative configure.xml file (at your own risk!)')

    parser.add_argument('--src_cfg', default=None,
                        help='Provide a source mIFU XML to transplant the bundle/guide allocations from')

    parser.add_argument('--pa', default=None, type=float,
                        help="""Set preferred mIFU PA""")

    parser.add_argument('--pressure', default=775.0, type=float,
                        help="""Set preferred atmospheric pressure""")
    
    parser.add_argument('--outdir', dest='output_dir', default='output',
                        help="""name for the directory which will contain the
                        combo catalogue""")

    parser.add_argument('--nodither', action='store_true',
                        help="""Do not run dither at all""")

    parser.add_argument('--noclobber', action='store_true',
                        help="""Do not overwrite existing files""")

    parser.add_argument('--fail_wipe', action='store_true',
                        help="""If *any* of the configure/dither operations fail, wipe *all* outputs (useful for pipelines)""")
    
    parser.add_argument('--nogui', action='store_true',
                        help="""Do not launch the GUI for mIFU (req. ---src_cfg)""")
    
    parser.add_argument('--nolint', dest='nolint', action='store_true',
                        help='Do not perform a lint')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')

    args = parser.parse_args()
    if args.nogui:
        assert args.src_cfg, "If no GUI is specified, an input --src_cfg file is required"

    if args.cfg_path == '':
        args.cfg_path = None
        
    logging.basicConfig(level=getattr(logging, args.log_level.upper()))
    configured_xmls = run_configure(args.xml_files,
                                    config_path=args.cfg_path,
                                    version=args.cfg_ver,
                                    dither=not args.nodither,
                                    clobber=not args.noclobber,
                                    lint=not args.nolint,
                                    src_cfg=args.src_cfg,
                                    nogui=args.nogui,
                                    set_pa=args.pa,
                                    set_pressure=args.pressure,
                                    cfg_override=args.cfg_override,
                                    wipe_on_fail=args.fail_wipe)
    if not args.nodither:
        dithered_xmls = run_dither(configured_xmls,
                                   config_path=args.cfg_path,
                                   version=args.cfg_ver,
                                   clobber=not args.noclobber,
                                   lint=not args.nolint,
                                   cfg_override=args.cfg_override,
                                   wipe_on_fail=args.fail_wipe)

    import time
    time.sleep(1.0)
    logging.info("Configure and dither complete")
