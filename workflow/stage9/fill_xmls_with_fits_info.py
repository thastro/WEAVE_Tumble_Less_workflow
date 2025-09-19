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
import numpy as _np
import ast

from astropy.io import fits
from astropy.table import Table
from workflow.utils.classes import OBXML

import time

import xml.dom.minidom as minidom
import xml.etree.ElementTree as ET

from collections import defaultdict
from collections import Counter

import traceback
import sys

from concurrent.futures import ProcessPoolExecutor

allowed_modes = ['L', 'O']


def get_group_dict(filename):
    """ A function to extract the linking parameters from the CSV file
    
    Args: filename Str: the name of the file
    
    Returns: group_data Dict: a dictionary with all of the linking.csv values
    """
    try:
        groups = Table.read(filename, comment='#')
    except Exception as e:
        logging.error('Input linking CSV file {} could not be read:'.format(filename))
        exc_type, exc_obj, tb = sys.exc_info()
        traceback.print_tb(tb)
        tb_info = traceback.extract_tb(tb)
        filename, line, func, text = tb_info[-1]
        logging.error('An exception occurred on line {}'.format(line))
        logging.error(e)
        raise SystemExit(1)

        
    group_data = {}
	
    for name in groups.colnames:
        try:
            group_data[name] = groups[name].tolist()
            if None in group_data[name]:
                group_data[name] = list(_np.array(group_data[name], dtype=str))
        except SystemExit: 
            logging.error('Information about {} not found'.format(name))


    return group_data

def check_linking_csv(group_data):
    """A function to check the validity of the linking CSV. This includes identifying non-integer priority values 
    and incorrect LINKING_MODE values.

    Args:
        group_data (dict): the data of the linking csv file

    Returns:
        bool: whether the csv file has passed the checks
    """
    
    allowed_modes = ['L', 'O']
    expected_elements = 5
    
    for key in group_data:
        if not group_data[key]:
            AssertionError("No value found for key")
    # Check if all keys have the same length of values
    for key, values in group_data.items():
        if not all(value is not None and value != '' for value in values):
            raise AssertionError(str(key)+" Key does not have the same number of values as requested")
    
    for key in group_data:
        if key=='KEY':
            for value in group_data[key]:
                indices = find_key_index(value, group_data)
                if len(indices)>1:
                    logging.warning(f'KEY {value} in linking.csv contains duplicate elements')
                    #return False
                split_elements = value.split(':')
                if not len(split_elements) == expected_elements:
                    raise AssertionError(f'KEY Structure {value} incorrect')
                
        elif key=='LINKING_MODE':
            for mode in group_data[key]:
                #logging.error('LINKING_MODE must be either L or O')
                if not mode in allowed_modes:
                    raise AssertionError('LINKING_MODE must be either L or O')
                #return False
                        
        elif key=='SEED_VALUE':
            for value in group_data[key]:
                #logging.error('Linking priority must be an integer number')
                if not type(value)==int:
                    raise AssertionError('Linking priority must be an integer number')
                #return False

    #counting how many times a key occurs in the CSV file
    full_keys = ["{}:{}".format(_key,_mode) for _key,_mode in zip(group_data['KEY'],group_data['LINKING_MODE'])]
    if len(full_keys) != len(_np.unique(full_keys)):
        key_counts = Counter(key for key in full_keys)
        for key in key_counts.keys():
            if key_counts[key] > 1:
                logging.error("KEY:LINKING_MODE '{}' appears {} times in CSV file".format(key,key_counts[key]))
        raise AssertionError('Duplication of linking keys under the same LINKING_MODE is not allowed in the CSV linking file')
                


def find_key_index(key, group_data):
    """A function to identify the index needed to extract the linking information for a single key.

    Args:
        key (Str): The key string being searched for
        group_data (dict): The dictionary containing the whole information from the linking csv file

    Returns:
        list: The indicies where the key string has been found. If unique len should equal 1.
    """
    indices = [index for index, value in enumerate(group_data['KEY']) if value == key]
    
    if len(indices)==0:
        logging.error('No matching key found in group_data')
        
    return indices

def find_key_info_data(key, group_data)-> dict:
    """A function to generate a dictionary with only the linking information which corresponds to a unique key

    Args:
        key (Str): The key string being searched for
        group_data (dict): The dictionary containing the whole information from the linking csv file

    Returns:
        dict: A dictionary containing the linking information which corresponds to the specific key
    """
    
    key_index = find_key_index(key, group_data)
    ###DO we want an assertion here or just stick with grabbing the first one?
    
    link_data = defaultdict(list)
    try:
        for jkey in key_index:
            for k in group_data:
                link_data[k].append(group_data[k][jkey])
            
    except:
        logging.error('No matching key found')
    
    return link_data

def get_unique_seed_value_by_group_id(data, target_group_id, mode):
    """A function to return the SEED_VALUE for a specific KEY entry within the group_alloc extracted
       from the linking.csv file
    

    Args:
        data Dict: the information extracted from the linking.csv file
        target_group_id Str: the GROUP_ID being searched for
        mode Str: The linking mode being used

    Raises:
        AssertionError: _description_

    Returns:
        _type_: _description_
    """
    seed_values = []

    # Iterate through the indices of GROUP_ID
    for index, group_id in enumerate(data['GROUP_ID']):
        if group_id == target_group_id and data['LINKING_MODE'][index] == mode:
            seed_values.append(data['SEED_VALUE'][index])

    # Check if all seed values are the same
    if seed_values and all(seed == seed_values[0] for seed in seed_values):
        logging.debug('seed_values are '+str(seed_values))
        return seed_values[0]
    else:
        logging.debug('seed_values are '+str(seed_values))
        logging.warning("SEED_VALUE values for the given GROUP_ID may not be the same.")
        try:
            return seed_values[0]
        except:
            raise AssertionError("SEED_VALUE values not returned")

    
def _clear_linking_attributes(obselement)->None:
    """A function to set all of the linking group attributes to an empty string

    Args:
        obselement (_type_): the observation element of the xml file
    """
    
    ###Need to check that I don't have to extract the element before I can set it
    blank_element = ''
    
    old_str_linkgroup = str(obselement.getAttribute('linkedgroup'))
    obselement.setAttribute('linkedgroup', blank_element)
    
    #old_str_linkgroup_prio = str(obselement.getAttribute('linkedgroup_priority'))
    obselement.setAttribute('linkedgroup_priority', blank_element)
    
    #old_str_groupname = str(obselement.getAttribute('obsgroup'))
    obselement.setAttribute('obsgroup', blank_element)
    
    #old_str_group_valid = str(obselement.getAttribute('obsgroup_validity'))
    obselement.setAttribute('obsgroup_validity', blank_element)
    
def _extract_obs_key_elements(obselement):
    obs_key = str(obselement.getAttribute('progtemp')) +':'+ str(obselement.getAttribute('obstemp'))
    
    return obs_key


def _modify_observation_element(obselement, data_dict):
    """A function to modify the linking attributes of the xml file based on the dictionary values it is given.

    Args:
        obselement <class 'xml.dom.minidom.Element'>: the observation element of the xml file to be modified
        data_dict Dict: the data dictionary that contains the values to populate the xml file with

    Raises:
        AssertionError: values for a key should be a single element not a list
    """
    
    for key in data_dict:
        value = data_dict.get(key)
        if isinstance(value, list):
            raise AssertionError("Value for key is a list, but should be a single element")
  
    if data_dict['LINKING_MODE'][0]=='L':
            
        old_str_groupname = str(obselement.getAttribute('linkedgroup'))
        new_str_groupname = str(data_dict['GROUP_ID'])
        obselement.setAttribute('linkedgroup', new_str_groupname)
        
        old_str_group_prio = str(obselement.getAttribute('linkedgroup_priority'))
        new_str_group_prio = str(data_dict['SEED_VALUE'])
        obselement.setAttribute('linkedgroup_priority', new_str_group_prio)
            
                
            
    if data_dict['LINKING_MODE'][0]=='O':
        old_str_groupname = str(obselement.getAttribute('obsgroup'))
        new_str_groupname = str(data_dict['GROUP_ID'])
        obselement.setAttribute('obsgroup', new_str_groupname)
        
        old_str_group_valid = str(obselement.getAttribute('obsgroup_validity'))
        new_str_group_valid = str(data_dict['SEED_VALUE'])
        obselement.setAttribute('obsgroup_validity', new_str_group_valid)
    

def find_best_group_modes(key_list, group_data, O_group_id_counter, L_group_id_counter):
    """A function to find the most dominant group from a given list of keys, after being given a set of data and empty counters for each LINKING_MODE.

    Args:
        key_list (_type_): A list of keys which matched with targets in the xml file
        group_data (_type_): the dictionary extracted from linking.csv
        O_group_id_counter (_type_): A dictionary with every 'O' mode GROUP_ID where the counters are set to zero
        L_group_id_counter (_type_): A dictionary with every 'L' mode GROUP_ID where the counters are set to zero

    Raises:
        AssertionError: asserts that the mode of the key must match the data being iterated through

    Returns:
        best_L_group Str: the most dominant group with a linkmode of 'L' for the xml file
        best_O_group Str: the most dominant group with a linkmode of 'O' for the xml file
         
    """
    for index, item in enumerate(key_list):
        mode = ':'.join(item.split(':')[-1])
        KEY = ':'.join(item.split(':')[:-1])
        data = dict(find_key_info_data(KEY, group_data))
        if data:
            if mode == 'L':
                # Filter based on LINKING_MODE=='O'
                filtered_values = {key: [value[i] for i in range(len(value)) if data['LINKING_MODE'][i] == 'L'] for key, value in data.items()}
                reference_key = 'KEY'
                reference_length = len(filtered_values[reference_key])
                for i in range(reference_length):
                    if not filtered_values['LINKING_MODE'][i]=='L':
                        raise AssertionError('Extracted data does not match key mode')
                    data_key = filtered_values['KEY'][i]
                    group_id = filtered_values['GROUP_ID'][i]
                    if group_id in L_group_id_counter:
                        L_group_id_counter[group_id] +=1
            if mode=='O':
                # Filter based on LINKING_MODE=='O'
                filtered_values = {key: [value[i] for i in range(len(value)) if data['LINKING_MODE'][i] == 'O'] for key, value in data.items()}
                reference_key = 'KEY'
                reference_length = len(filtered_values[reference_key])
                for i in range(reference_length):
                    if not filtered_values['LINKING_MODE'][i]=='O':
                        raise AssertionError('Extracted data does not match key mode')
                    data_key = filtered_values['KEY'][i]
                    group_id = filtered_values['GROUP_ID'][i]
                    if group_id in O_group_id_counter:
                        O_group_id_counter[group_id] +=1
                        
    logging.info('L_group_id_counter '+str(L_group_id_counter))
    logging.info('O_group_id_counter '+str(O_group_id_counter))
    #Finding the best one after counters modified
    # Initialize variables to track the best value and corresponding key
    best_value_L = float('-inf')  # Start with negative infinity to handle negative values
    best_group_L = None
    
    best_value_O = float('-inf')
    best_group_O = None

    # Iterate through the dictionary and find the best value
    for group, num in L_group_id_counter.items():
        if best_group_L is None or (num > best_value_L and best_group_L is not None):
            best_value_L = num
            best_group_L = group
        elif num == best_value_L:
            best_group_L = max(best_group_L, group, key=lambda x: get_unique_seed_value_by_group_id(group_data, x, 'L'))
            
    if best_value_L == 0:
        best_group_L = None
    
    for group, num in O_group_id_counter.items():
        if best_group_O is None or (num > best_value_O and best_group_O is not None):
            best_value_O = num
            best_group_O = group
        elif num == best_value_O:
            best_group_O = min(best_group_O, group, key=lambda x: get_unique_seed_value_by_group_id(group_data, x, 'O'))
    
    if best_value_O == 0:
        best_group_O = None
                
    logging.debug("Best value LINKING_MODE L: {}".format(best_value_L))
    logging.debug("GROUP_ID with the best value (best SEED_VALUE in case of a tie): {}".format(best_group_L) if best_group_L is not None else "None")
    
    logging.debug("Best value LINKING_MODE O: {}".format(best_value_O))
    logging.debug("GROUP_ID with the best value (best SEED_VALUE in case of a tie): {}".format(best_group_O) if best_group_O is not None else "None")
    return best_group_L, best_group_O


def find_entries_by_group_id(data, target_group_id, mode):
    """A function to find all entries from linking.csv for a specific GROUP_ID value

    Args:
        data Dict: the dictionary of linking.csv values
        target_group_id Str: the GROUP_ID string value
        mode Str: The linking mode to search under, either 'L' or 'O'

    Returns:
        macthing_entries list: a list containing the matching dictionary entries
    """
    # Initialize a list to store matching entries
    matching_entries = []

    # Iterate through the data and check for matches
    for index, group_id in enumerate(zip(data['GROUP_ID'], data['LINKING_MODE'])):
        if group_id[0] == target_group_id and group_id[1]==mode:
            matching_entries.append({
                'KEY': data['KEY'][index],
                'LINKING_MODE': data['LINKING_MODE'][index],
                'GROUP_ID': group_id[0],
                'SEED_VALUE': data['SEED_VALUE'][index]
            })

    return matching_entries

def check_seed_values(matching_entries):
    """A function to check the seed values for a list of entries

    Args:
        macthing_entries list: a list containing the matching dictionary entries

    Raises:
        AssertionError: cannot have multiple SEED_VALUES for the same GROUP_ID
    """
    # Get the SEED_VALUE of the first entry
    first_seed_value = matching_entries[0]['SEED_VALUE']

    # Check if all matching entries have the same SEED_VALUE
    for entry in matching_entries[1:]:
        if entry['SEED_VALUE'] != first_seed_value:
            logging.error("SEED_VALUE mismatch detected!")
            raise AssertionError("SEED_VALUE mismatch detected!")


def open_FITS_file(fits_cat, add_sim):
    """A function to open the catalogue FITS file generated by the WASP with the CNAMES included

    Args:
        fits_cat (Str): The catalogue FITS file output from the WASP
        add_sim (Bool): a flag which states whether simulation data is included in the FITS header

    Returns:
        fits_data (astropy.io.fits.hdu.hdulist.HDUList): the FITS data extracted using astropy.io.fits
    """
    
    with fits.open(fits_cat) as hdu_list:
        
        fits_data = hdu_list[1].data
        
        # Let's check that there are CNAMES in here (a good indication that it
        # has been WASP-validated!)
        if not _np.all(_np.char.startswith(fits_data['CNAME'],'WVE_')):
            raise AssertionError("Not all CNAMEs have been populated in provided FITS catalogue: has this been validated with WASP?")

        
        if add_sim == True:
            try:
                sim_data = hdu_list[2].data
            except Exception as e:
                logging.error('Unable to open sim data from FITS catalogue, please check the correct FITS file is being used')
                logging.error('Exception occurred: '+e)
            
                exc_type, exc_obj, tb = sys.exc_info()
                traceback.print_tb(tb)
                tb_info = traceback.extract_tb(tb)
                filename, line, func, text = tb_info[-1]
                logging.error('An exception occurred on line {}'.format(line))
                logging.error(e)
                raise SystemExit(1)
        else:
            sim_data = None
        
        return fits_data, sim_data

def generate_output_file_list(xml_files, output_dir):
    """A function to generate the output file list for a series of XML files

    Args:
        xml_files (list): the list of XML files
        output_dir (Str): the directory the output files will be located in

    Returns:
        list: a list of combined XML file names with the corresponding output directory
    """
    return [os.path.join(output_dir, '{}.xml'.format(os.path.basename(xml_file).rsplit('-', 1)[0])) for xml_file in xml_files]

def generate_output_file(xml_file, output_dir):
    """A function to generate the output file path for a specific XML file

    Args:
        xml_files (str): the XML file
        output_dir (Str): the directory the output files will be located in

    Returns:
        str: the combined XML file name with the corresponding output directory
    """
    
    return os.path.join(output_dir, '{}.xml'.format(os.path.basename(xml_file).rsplit('-', 1)[0]))

def outfile_overwrite(output_path, overwrite):
    """
        A function to check if a output file already exists and overwrite it on request
    Args:
        output_path (str): a string stating the output file path
    """
    if os.path.exists(output_path):
            
            if overwrite == True:
                logging.info('Removing previous file: {}'.format(output_path))
                os.remove(output_path)
            else:
                logging.info(
                    'Skipping file because it already exists: {}'.format(
                        output_path))

    else:
        logging.debug('path does not exist')
    return


######functions here to organise the ledger printing
def combine_seed_map(group_data):
    """A function to format the group_data so that it can be used in the final summary printing.

    Args:
        group_data Dict: the dictionary of linking.csv values

    Returns:
        combined_key_seed_map: a dictionary of the formatted entries
        combined_key_group_id_map: a dictionary of the formatted entry ids
    """
    try:
        combined_key_seed_map = {}
        combined_key_group_id_map = {}
        # Iterate through the keys and create the combined keys
        for key, seed_value, linking_mode, group_id in zip(group_data['KEY'], group_data['SEED_VALUE'], group_data['LINKING_MODE'], group_data['GROUP_ID']):
            combined_key = f'{key}:{linking_mode}'
            combined_key_seed_map[combined_key] = seed_value
            combined_key_group_id_map[combined_key] = group_id
            
    except Exception as e:
        logging.error('Unable to recombine group_data for printing formatting')
        logging.error('Exception occurred: '+ e)
            
        exc_type, exc_obj, tb = sys.exc_info()
        traceback.print_tb(tb)
        tb_info = traceback.extract_tb(tb)
        filename, line, func, text = tb_info[-1]
        logging.error('An exception occurred on line {}'.format(line))
        logging.error(e)
        raise SystemExit(1)
    
    return combined_key_seed_map, combined_key_group_id_map



def print_ledger_key_files(data, key_seed_map):
    """A function to print the KEYS, GROUP_IDS, and SEED_VALUES for a given xml file based on the final linking ledger.

    Args:
        data Dict: the ledger dictionary produced by the xml file processing
        key_seed_map Dict: the dictionary of linking.csv values
    """
    unique_key_files = {}
    combined_key_seed_map, combined_key_group_id_map = combine_seed_map(key_seed_map)
        
    # Iterate through the data
    for filename, key_value_dict in data.items():
        for key, counter in key_value_dict.items():
            if counter > 0:
                if key not in unique_key_files:
                    unique_key_files[key] = []
                unique_key_files[key].append(filename)


    # What keys from the CSV had no counterparts in the provided XMLs?
    csv_keys = _np.unique(key_seed_map['KEY'])
    unused = []
    for partial_key in csv_keys:
        used = False
        for mode in allowed_modes:
            key = "{}:{}".format(partial_key,mode)
            if key in unique_key_files:
                used = True
        if not used:
            unused.append(partial_key)

    unused = _np.unique(unused)
    if len(unused):
        logging.info("\tUser-defined linking keys that had NO match within the supplied XML files:")
        for key in unused:
            logging.info("\t\t{}".format(key))
        logging.info("")
            
    # Print the key linkage results
    for key, filenames in unique_key_files.items():
        seed_value = combined_key_seed_map[key]
        group_id = combined_key_group_id_map[key]
        logging.info("\tKey '{}' with a GROUP_ID of {} and a SEED_VALUE of {} is found in the following files:".format(key,group_id,seed_value))

        for filename in filenames:
            logging.info('\t\t'+filename)
        logging.info("")

    
def print_ledger_group_files(data, key_seed_map):
    """A function to print the XML files assigned a GROUP_IDS and SEED_VALUES based on the final linking ledger.

    Args:
        data Dict: the ledger dictionary produced by the xml file processing
        key_seed_map Dict: the dictionary of linking.csv values
    """
    unique_key_files = {}
    combined_key_seed_map, combined_key_group_id_map = combine_seed_map(key_seed_map)

    group_membership = {}
    group_seedlookup = {}
    
    # Iterate through the data
    for filename, key_value_dict in data.items():
        for key, counter in key_value_dict.items():
            if counter > 0:
                if key not in unique_key_files:
                    unique_key_files[key] = []
                unique_key_files[key].append(filename)

    # Print the results
    for key, filenames in unique_key_files.items():
        seed_value = combined_key_seed_map[key]
        group_id = combined_key_group_id_map[key]

        try:
            group_membership[group_id] += filenames
        except:
            group_membership[group_id] = filenames
            group_seedlookup[group_id] = seed_value

    for grp in _np.unique(key_seed_map['GROUP_ID']):
        if grp in group_membership:
            logging.info("\tGroup '{}' (seed: {}) has been connected to the following OB files:".format(grp,group_seedlookup[grp]))
            for fname in group_membership[grp]:
                logging.info("\t\t{}".format(fname))
        else:
            logging.info("\tGroup '{}' has NO assignments to any provided OB files".format(grp))
        logging.info("")



        

def process_xml_file(xml_file, fits_data, sim_data, group_alloc, clean_non_allocated, output_dir, overwrite, strict, ledger=None, multiprocess=False):
    """A function to process the XML file and populate it with the appropriate CNAMES generated by the FITS catalogue file.
        If linking is enabled, it will tag all files with the most dominant GROUP_ID from the list of matching target keys.
        If strict is enabled, multiple files cannot be tagged with more than one group element per LINKING_MODE.

    Args:
        xml_file (Str): the XML filename
        fits_data astropy.io.fits.hdu.hdulist.HDUList: the FITS data extracted using astropy.io.fits
        sim_data astropy.io.fits.hdu.hdulist.HDUList: the simulated FITS data extracted using astropy.io.fits
        group_alloc tuple: Sets any relational linkage between the XML files generated. If empty, no groups are made"
        clean_non_allocated Bool: remove non-allocated targets in the output files
        ob_xml_output list: a list to append the combined XML file names with the corresponding output directory
        multiprocess Bool: Indicates this function is running under multiprocessing - returns both ledger and output_path
        strict Bool: Prevents multiple groups being tagged to the same xml_file.
    """
    
    #setting start time
    iter_start_time = time.time()
    
    # Create an OB from the XML file
        
    ob_xml = OBXML(xml_file)
    
    ledger[xml_file] = {}
    
    #grabbing the obs linking key info if linking enabled and clear the attributes in both cases
    if group_alloc:
        for element in ob_xml.dom.getElementsByTagName('observation'):
            obs_mode = str(element.getAttribute('obs_mode'))
            _clear_linking_attributes(element)
            obs_key = _extract_obs_key_elements(element) + ':'
            
            logging.debug("obs_key "+obs_key)
        
        #checking the apply dither option
        for element in ob_xml.dom.getElementsByTagName('dithering'):
            ifu_dither = int(element.getAttribute('apply_dither'))
        
    else:
        obs_mode = False
        for element in ob_xml.dom.getElementsByTagName('observation'):
            _clear_linking_attributes(element)
    
    # Populate the OB XML with the FITS data
    
    targ_keys = ob_xml.populate_targets_with_fits_data(fits_data, sim_data=sim_data,
                                            clean_non_allocated=clean_non_allocated, linking=obs_mode)
    

    logging.debug("targ_keys "+str(targ_keys))
    
    # Populate the OB XML with grouping data
    #ob_xml.populate_targets_with_group_data(group_alloc)
    
    
    ## Probably worth trying to abstract this into it's own function
    if group_alloc and isinstance(targ_keys, list) and len(targ_keys)>0:
        if multiprocess:
            # logging.info("")
            logging.info("Linking analysis on {}".format(os.path.basename(xml_file)))

        try:
            #skipping process if no target keys extracted from the file
            if not len(targ_keys)>=1:
                raise AssertionError('Cannot process an empty target list')
            
            #Extract keys from group_data dictionary
            #combining LINKING_MODE with the rest of the CSV KEY
            csv_keys = [f"{key}:{linkmode}" for key, linkmode in zip(group_alloc['KEY'], group_alloc['LINKING_MODE'])] 
            # To cut down key back again -> ':'.join(key.split(':')[:-1])
            
            #Setting up dictionary structures to track the GROUP_IDS that occur within this file
            # Iterate through all instances of GROUP_ID and categorize based on LINKING_MODE
            # Initialize dictionaries
            L_group_dict = {}
            O_group_dict = {}
            
            for group_id, linking_mode in zip(group_alloc['GROUP_ID'], group_alloc['LINKING_MODE']):
                if linking_mode == 'L':
                    L_group_dict[group_id] = 0  # Add to L_group_dict
                elif linking_mode == 'O':
                    O_group_dict[group_id] = 0  # Add to O_group_dict

            logging.debug("L_group_dict:" +str(L_group_dict))
            logging.debug("O_group_dict:" +str(O_group_dict))
            
            
            #counting how many times a key occurs in the CSV file
            key_counts = Counter(key for key in csv_keys)
            
            # #Adding assertion that a key can't appear more than once for strict mode
            # ##May want to move this into check_csv function
            # if not all(count <= 1 for count in key_counts.values()):
            #     raise AssertionError('Duplicate linking key is not allowed in this CSV file')
            
            logging.debug("CSV key counts "+str(key_counts))
            
            # Setting the CSV keys inside the ledger dict for this file:
            ledger[xml_file] = {value: 0 for value in set(csv_keys)} #using a set here to make sure the keys aren't repeated in ledger

            
            #combining full key together for cross matching from obs and target elements
            targ_keys = [str(obs_key)+str(item)+str(ifu_dither) for index, item in enumerate(targ_keys)]
            #Listing all possible target keys by considering either mode
            full_targ_keys = [item+':'+ m for index, item in enumerate(targ_keys) for m in allowed_modes]

            logging.info("Available keys for linking:")
            for key in targ_keys:
                logging.info("\t{}".format(key))
                
            
            #Adding to the ledger counter if any target keys match those in the CSV file
            for key in ledger[xml_file].keys():
                matches = [obj for obj in full_targ_keys if key in obj]
                ledger[xml_file][key] += len(matches)
                
                logging.debug("key_prefix "+ key+  ' '+str(len(matches)))

            # Check if all counters are zero
            all_zero_counters = all(counter == 0 for counter in ledger[xml_file].values())
            
            #grabbing the list of keys that matched in this file:
            keys_occurred = [key for key, value in ledger[xml_file].items() if value > 0]
            logging.debug('Matching keys: '+str(keys_occurred))
            

            if not all_zero_counters:
                #Finding the most grouping in the list
                dominant_L_group, dominant_O_group = find_best_group_modes(keys_occurred, group_alloc, O_group_dict, L_group_dict)
                logging.debug("Debug info: Dominant GROUP_IDS "+'L '+str(dominant_L_group)+' and O '+str(dominant_O_group))
                
                #Putting in strict condition that means multiple dominant groups cannot be found in the same file.
                #Or multiple groups at all in a single file.
                if strict:
                    #considering L and O case seperately
                    #L case:
                    #First strict case: cannot have multiple groups going to the same file:
                    count_positive_values = sum(1 for value in L_group_dict.values() if value > 0)
                    if not count_positive_values <= 1:
                        raise AssertionError("More than one grouping applies to this file, not allowed in STRICT mode.")

                    if 0:
                        #Second case: considering the strict condition where we cannot have multiple dominant groups
                        L_max_value = max(L_group_dict.values())
                        # Count the occurrences of the maximum value
                        count_max_value = sum(1 for value in L_group_dict.values() if value == L_max_value)
                        if not count_max_value == 1 and not L_max_value == 0:
                            logging.debug('L_group_dict ->'+str(L_group_dict))
                            raise AssertionError("Multiple dominant groups occur in this file: Not allowed in STRICT mode.")
                    
                    #O case
                    count_positive_values = sum(1 for value in O_group_dict.values() if value > 0)
                    if not count_positive_values <= 1:
                        raise AssertionError("More than one grouping applies to this file, not allowed in STRICT mode.")
                    
                    if 0:
                        O_max_value = max(O_group_dict.values())
                        # Count the occurrences of the maximum value
                        count_max_value = sum(1 for value in O_group_dict.values() if value == O_max_value)
                        if not count_max_value == 1 and not O_max_value == 0:
                            logging.debug('O_group_dict ->'+str(O_group_dict))
                            raise AssertionError("Multiple dominant groups occur in this file: Not allowed in STRICT mode.")
                
                #extracting the information for the dominant groups
                L_entries = find_entries_by_group_id(group_alloc, dominant_L_group, mode='L')
                O_entries = find_entries_by_group_id(group_alloc, dominant_O_group, mode='O')
                
                #Reducing entries to only those possible for the file:
                L_entries = [item for index, item in enumerate(L_entries) if item['KEY'] in targ_keys]
                O_entries = [item for index, item in enumerate(O_entries) if item['KEY'] in targ_keys]
                
                #checking if multiple SEED_VALUES occur for the dominant group, asserts they must be the same
                #checking if all values are consistent across the dominant group entries, except for the KEY value
                if L_entries:
                    try:
                        check_seed_values(L_entries)
                        reference_values = {k: L_entries[0][k] for k in ['LINKING_MODE', 'GROUP_ID', 'SEED_VALUE']}
                        all_same_values = all(all(d[k] == reference_values[k] for k in ['LINKING_MODE', 'GROUP_ID', 'SEED_VALUE']) for d in L_entries[1:])
                        if not all_same_values:
                           raise AssertionError("Inconsistent values across dominant group entries: "+L_entries[0])
                    except:
                        raise AssertionError("For LINKING_MODE L entries: Mismatched entries detected!")
                if O_entries:
                    try:
                        check_seed_values(O_entries)
                        reference_values = {k: O_entries[0][k] for k in ['LINKING_MODE', 'GROUP_ID', 'SEED_VALUE']}
                        all_same_values = all(all(d[k] == reference_values[k] for k in ['LINKING_MODE', 'GROUP_ID', 'SEED_VALUE']) for d in O_entries[1:])
                        if not all_same_values:
                            raise AssertionError("Inconsistent values across dominant group entries: "+L_entries[0])
                    except:
                        raise AssertionError("For LINKING_MODE O entries: SEED_VALUE mismatch detected!")
                
                #combining all entries
                combined_entries = L_entries + O_entries
                logging.debug('combined_entries:'+str(combined_entries))

                #start tagging the dominant information
                linked_modes = []
                    
                for entry in combined_entries:    
                    if entry['LINKING_MODE'] in linked_modes:
                        continue
                    logging.debug("Debug info: entry"+ str(entry))
                            
                    for element in ob_xml.dom.getElementsByTagName('observation'):
                        logging.info('MODIFYING file {} WITH VALUES {}'.format(xml_file,entry))
                        _modify_observation_element(element, entry)
                        linked_modes.append(entry['LINKING_MODE'])
                        
            elif all_zero_counters:
                logging.info('No matching keys found in this file')


        except AssertionError as e:
            exc_type, exc_obj, tb = sys.exc_info()
            traceback.print_tb(tb)
            tb_info = traceback.extract_tb(tb)
            filename, line, func, text = tb_info[-1]
            logging.error('An exception occurred on line {}'.format(line))
            logging.error(e)
            logging.info('Writing current file state')
            output_path = generate_output_file(xml_file, output_dir)
    
            outfile_overwrite(output_path, overwrite)
    
            ob_xml.write_xml(output_path, replace_tabs=True)
            if not multiprocess:
                raise AssertionError("Assertion error occurred. "+str(e)+" Stopping the loop.")
            return None, None, None
    
    output_path = generate_output_file(xml_file, output_dir)
    
    outfile_overwrite(output_path, overwrite)
    
    ob_xml.write_xml(output_path, replace_tabs=True)
    
    
    #Printing the timings
    iter_end_time = time.time()
    iter_execution_time = iter_end_time - iter_start_time
    logging.info("Execution time: {:.2f} seconds".format(iter_execution_time)) 
    
    if multiprocess:
        logging.info("")
        return output_path, ledger, targ_keys
    
    return output_path, targ_keys


def run_linear_processing(xml_files, fits_data, sim_data, group_alloc, clean_non_allocated, output_dir, overwrite, strict, ledger):
    """A function to run the processing of the xml_files linearly. This inserts the CNAME's and addes additional linkage.

    Args:
        xml_file (Str): the XML filename
        fits_data astropy.io.fits.hdu.hdulist.HDUList: the FITS data extracted using astropy.io.fits
        sim_data astropy.io.fits.hdu.hdulist.HDUList: the simulated FITS data extracted using astropy.io.fits
        group_alloc tuple: Sets any relational linkage between the XML files generated. If empty, no groups are made"
        clean_non_allocated Bool: remove non-allocated targets in the output files
        output_dir Str: the directory to store the output files
        strict Bool: Prevents multiple groups being tagged to the same xml_file.
        ledger Dict: The ledger where the processing information is tracked.


    Returns:
        output_paths list: A list of the processed xml file output_paths
    """
    results = []
    output_paths = []
    available_keys = {}
    for xml_file in xml_files:
        output_path, available_key = process_xml_file(xml_file, fits_data, sim_data, group_alloc, clean_non_allocated,
                                                       output_dir, overwrite, strict, ledger)
        available_keys[xml_file] = available_key
        output_paths.append(output_path)
        
    # output_paths = [result.get() for result in results]
    
    return output_paths, available_keys

def run_multiprocessing(xml_files, fits_data, sim_data, group_alloc, clean_non_allocated, output_dir, overwrite, strict, ledger, groups, ncpu=2):
    """A function to run the processing of the xml_files over multiple computer cores. This inserts the CNAME's and addes additional linkage.

    Args:
        xml_file (Str): the XML filename
        fits_data astropy.io.fits.hdu.hdulist.HDUList: the FITS data extracted using astropy.io.fits
        sim_data astropy.io.fits.hdu.hdulist.HDUList: the simulated FITS data extracted using astropy.io.fits
        group_alloc tuple: Sets any relational linkage between the XML files generated. If empty, no groups are made"
        clean_non_allocated Bool: remove non-allocated targets in the output files
        output_dir Str: the directory to store the output files
        strict Bool: Prevents multiple groups being tagged to the same xml_file.
        ledger Dict: The ledger where the processing information is tracked.


    Returns:
        output_paths list: A list of the processed xml file output_paths
        _ledger Dict: a dictionary containing the processed linking information for each file
    """

    results = []
    ledgers = []
    available_keys = {}

    jobs = []
    
    with ProcessPoolExecutor(max_workers=ncpu) as executor:
        for xml_file in xml_files:
            job = executor.submit(process_xml_file, xml_file, fits_data, sim_data, group_alloc, clean_non_allocated, output_dir, overwrite, strict, ledger, True)
            jobs.append(job)

    for job, xml_file in zip(jobs, xml_files):
        output_path, ledger, available_key = job.result()
        available_keys[xml_file] = available_key
        if not output_path:
            if len(group_alloc):
                logging.error("Multiprocess returned errors: not all OB links were necessarily implemented - aborting. CONTINUE AT YOUR OWN RISK")
            else:
                logging.error("Multiprocess returned errors - aborting")
            raise SystemExit(1)

        results.append(output_path)
        ledgers.append(ledger)



    logging.debug('RESULTS: {}'.format(results))
    logging.debug('RESULTS: {}'.format(ledgers))
    logging.debug('RESULTS: {}'.format(available_keys))

    # Get the output paths from the results
    output_paths = [result for result in results] #[result.get() for result in results]

    return output_paths, ledgers, available_keys






def fill_xmls_with_fits_info(fits_cat, xml_files, output_dir,
                             group_data=None, clean_non_allocated=True,
                             add_sim=False, overwrite=False):
    """
    Fill a set of XML files with the information from a FITS catalogue.
    
    Parameters
    ----------
    fits_cat : str
        Name of a FITS file containing a FITS catalogue.
    xml_files : list of str
        A list of strings with the filenames of the XML files.
    output_dir : str
        Name of the directory which will contain the output XML files.
    group_data : str, optional
        Location of the CSV file containing grouping data for these XML files.
    clean_non_allocated: bool, optional
        Remove non-allocated targets in the output files.
    add_sim : bool, optional
        Add simulation elements to the targets.
    overwrite : bool, optional
        Overwrite the output XML files.

    Returns
    -------
    output_file_list : list of str
        A list with the output XML files.
    """
        
    # Open the combo FITS catalogue and read the needed data

    with fits.open(fits_cat) as hdu_list:
        
        fits_data = hdu_list[1].data
        
        if add_sim == True:
            sim_data = hdu_list[2].data
        else:
            sim_data = None

    # If it exists, open the group data file and read the data

    group_alloc = {}
    if group_data:
        group_alloc = get_group_dict(group_data)
        


            
    # Let's check that there are CNAMES in here (a good indication that it
    # has been WASP-validated!)
    assert _np.all(_np.char.startswith(fits_data['CNAME'],'WVE_')), "Not all CNAMEs have been populated in provided FITS catalogue: has this been validated with WASP?"

            
    # Create a list for saving the names of the output files
    
    output_file_list = []

    # For each XML file
    logging.info("Updating a total of {} XML files using CNAME'd catalogue {}".format(len(xml_files),fits_cat))
    index = 0

    for xml_file in xml_files:
        iter_start_time = time.time()
        # Create an OB from the XML file
        
        ob_xml = OBXML(xml_file)
        
        # Choose an output filename from the input filename
        
        input_basename = os.path.basename(xml_file)

        output_basename_wo_ext = '-'.join(input_basename.split('-')[:-1])
        
        output_basename = '{}.xml'.format(output_basename_wo_ext)
        
        output_path = os.path.join(output_dir, output_basename)
        
        # Save the name of the output file
        
        output_file_list.append(output_path)

        index += 1
        logging.info('\nPopulating targets for file {} [{}/{}]'.format(
            input_basename,index,len(xml_files)))
        
        # If the output file exists, remove or skip it

        if os.path.exists(output_path):
            
            if overwrite == True:
                logging.info('Removing previous file: {}'.format(output_path))
                os.remove(output_path)
            else:
                logging.info(
                    'Skipping file because it already exists: {}'.format(
                        output_path))
                continue

        # Populate the OB XML with the FITS data

        ob_xml.populate_targets_with_fits_data(
            fits_data, sim_data=sim_data,
            clean_non_allocated=clean_non_allocated)
        
        # Populate the OB XML with grouping data

        logging.info('Setting groups for file {}'.format(output_basename))

        ob_xml.populate_targets_with_group_data(
            group_alloc)

        
        # Write the OB XML to a file
        
        ob_xml.write_xml(output_path, replace_tabs=True)
        
        #checking the timings
        iter_end_time = time.time()
        iter_execution_time = iter_end_time - iter_start_time
        logging.info("Iteration Execution time: {:.2f} seconds".format(iter_execution_time))
    return output_file_list



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Fill XMLs with FITS info')

    parser.add_argument('fits_cat',
                        help='a FITS file containing a combo FITS catalogue')

    parser.add_argument('xml_file', nargs='+', help='name of a XML file')

    parser.add_argument('--groups', dest='group_data', default=None,
                        help="""Name of the grouping data file that sets any
                        relational linkage between the XML files generated. If 
                        empty, no groups are made""")
    
    parser.add_argument('--outdir', dest='output_dir', default='output',
                        help="""name of the directory which will contain the
                        output XML files""")

    parser.add_argument('--clean_non_allocated', action='store_true',
                        help='remove non-allocated targets in the output files')
    parser.add_argument('--keep_non_allocated', dest='clean_non_allocated',
                        action='store_false',
                        help='keep non-allocated targets in the output files')
    parser.set_defaults(clean_non_allocated=True)

    parser.add_argument('--add_sim', action='store_true',
                        help='add simulation elements to the targets')

    parser.add_argument('--overwrite', action='store_true',
                        help='overwrite the output files')

    parser.add_argument('--log_level', default='info',
                        choices=['debug', 'info', 'warning', 'error'],
                        help='the level for the logging messages')
    
    parser.add_argument('--strict', action='store_true', help="""Enforce restrictions on the group linking 
                        modifications.""")
    
    parser.add_argument('--ncpu', type=int, default=1, help="""Set to >1 to enable multiproccesing to run for stage 9""") 


    args = parser.parse_args()

    if args.ncpu > 1:
        logging_format = '%(asctime)s - %(levelname)s (@%(process)d) - %(message)s'
        logging.basicConfig(level=getattr(logging, args.log_level.upper()), format=logging_format)
    else:
        logging.basicConfig(level=getattr(logging, args.log_level.upper()))

    if not os.path.exists(args.output_dir):
        logging.info('Creating the output directory')
        os.mkdir(args.output_dir)

    start_time = time.time()
    
    fits_data, sim_data = open_FITS_file(args.fits_cat, args.add_sim)
    
    
    group_alloc = {}
    if args.group_data:
        group_alloc = get_group_dict(args.group_data)
        check_linking_csv(group_alloc)

    #Old printing for multiprocessing, it doesn't really work so switching to in the loop method
    
    if args.ncpu > 1:
        logging.info("Using multiprocessing with {} CPUs".format(args.ncpu))
        
        ledger = {}
        _, shared_ledger, available_keys = run_multiprocessing(args.xml_file, fits_data, sim_data, group_alloc, 
                                                               args.clean_non_allocated, args.output_dir,
                                                               args.overwrite, args.strict, ledger, args.group_data,
                                                               ncpu=args.ncpu)

        if args.group_data:
            list_as_string = str(shared_ledger)
            list_of_dicts = ast.literal_eval(list_as_string)

            # Merge the dictionaries into a single dictionary
            merged_dict = {}
            for item in list_of_dicts:
                merged_dict.update(item)

            ledger = merged_dict
                
    else:
        ledger = {}
        _, available_keys = run_linear_processing(args.xml_file, fits_data, sim_data, group_alloc, 
                                                  clean_non_allocated=args.clean_non_allocated, 
                                                  output_dir=args.output_dir, overwrite=args.overwrite,
                                                  strict=args.strict, ledger=ledger)
        
    
        
    if args.group_data:
        logging.info("")
        logging.info("###################")        
        logging.info('SUMMARY of LINKING:')
        logging.info("###################")
        for fname, data in ledger.items():
            logging.info("")
            logging.info("\tFilename: {}".format(fname))

            logging.info("\t\tKeys within this file:")
            for key in available_keys[fname]:
                logging.info("\t\t\t{}".format(key))

            # for key in ledger[fname]['']
            
            if not sum(data.values()):
                logging.info("\t\tNo matching keys in this OB")
            else:
                logging.info("\t\tKey usage information:")
                for key, value in data.items():
                    if value==0:
                        logging.debug("\t\t\tKEY  {} -> NOT MATCHED to file".format(key))
                    else:
                        logging.info("\t\t\tKEY  {} -> MATCHED to file".format(key))                 

        logging.info("")
        logging.info("##################")        
        logging.info('KEY USAGE SUMMARY:')
        logging.info("##################")
        logging.info("")
        print_ledger_key_files(ledger, group_alloc)


        logging.info("")
        logging.info("####################")        
        logging.info('GROUP USAGE SUMMARY:')
        logging.info("####################")
        logging.info("")
        print_ledger_group_files(ledger, group_alloc)
        
        
    end_time = time.time()
    #calculating and printing the total function run time
    execution_time = end_time - start_time

    logging.info("Execution time: {:.2f} seconds".format(execution_time))

