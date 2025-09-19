import xml.etree.ElementTree as ET
import os
import argparse
import pandas as pd

calib_df = pd.DataFrame({
            'arm': ['red', 'blue', 'red', 'blue', 'red', 'blue'],
            'cal_lamp': ['QTH', 'QTH', 'ThArCr', 'ThArCr', 'ThArCr', 'ThArCr'],
            'cal_lamp_filter_A': ['WCALU_clear1', 'WCALU_clear1', 'WCALU_clear1',
                                  'WCALU_clear1', 'WCALU_clear1', 'WCALU_clear1'],
            'cal_lamp_filter_B': ['EMPTY', 'EMPTY', 'EMPTY', 'EMPTY', 'EMPTY', 'EMPTY'],
            'exp_time': [3, 9, 10, 120, 10, 120],
            'ff_ilu_1': [0, 0, 0, 0, 0, 0],
            'ff_ilu_2': [0, 0, 0, 0, 0, 0],
            'ff_ilu_3': [0, 0, 0, 0, 0, 0],
            'order': [1, 1, 2, 2, 2, 2],
            'speed': ['fast', 'fast', 'fast', 'fast', 'fast', 'fast'],
            'type': ['fibre_lampflat', 'fibre_lampflat', 'arc', 'arc', 'arc', 'arc'],
})


def indent(elem, level=0):
    indent_size = "    "
    i = "\n" + level * indent_size
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + indent_size
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def save_xml(xml_tree, filename):
    indent(xml_tree.getroot())
    pretty_xml = ET.tostring(xml_tree.getroot(), encoding="unicode")
    print(pretty_xml)
    with open(filename, 'w') as file:
        file.write(pretty_xml)


def split_xml_by_field(xml_root, output_dir, original_filename):
    os.makedirs(output_dir, exist_ok=True)

    fields = xml_root.findall('.//fields/field')
    exposures = xml_root.findall('.//programme/exposures/exposure')
    observation = xml_root.findall('.//observation')
    configure = xml_root.findall('.//observation/configure')
    focal_plane_map = xml_root.findall('.//observation/configure/focal_plane_map')
    telescope = xml_root.findall('.//observation/configure/telescope')
    hour_angle_limits = xml_root.findall('.//observation/configure/hour_angle_limits')
    rotator_axis = xml_root.findall('.//observation/configure/focal_plane_map/rotator_axis')
    optical_axis = xml_root.findall('.//observation/configure/focal_plane_map/optical_axis')
    distortion_coefficients = xml_root.findall('.//observation/configure/focal_plane_map/' +
                                               'distortion_coefficients')
    rotator_centre = xml_root.findall('.//observation/configure/focal_plane_map/' +
                                      'rotator_centre')
    conditions = xml_root.findall('.//observation/configure/conditions')
    surveys = xml_root.findall('.//observation/surveys/survey')

    pa = observation[0].attrib.get('pa')
    name = observation[0].attrib.get('name')
    program = observation[0].attrib.get('program')

    dithering_attrributes = {'apply_dither': '0', 'dither_version': '%%%'}

    configure_attributes = configure[0].attrib
    telescope_attributes = telescope[0].attrib
    hour_angle_limits_attributes = hour_angle_limits[0].attrib

    focal_plane_map_attributes = focal_plane_map[0].attrib
    rotator_axis_attributes = rotator_axis[0].attrib
    optical_axis_attributes = optical_axis[0].attrib
    distortion_coefficients_attributes = distortion_coefficients[0].attrib
    rotator_centre_attributes = rotator_centre[0].attrib

    conditions_attributes = conditions[0].attrib
    #max_calibration = configure[0].attrib.get('max_calibration')
    #maximum_gate_angle = configure[0].attrib.get('maximum_gate_angle')
    #plate = configure[0].attrib.get('plate')


    #print(len(targets))
    if not fields:
        return "No <field> elements found within <fields> in the provided XML."

    file_counter = 1
    for exposure, field in zip(exposures, fields):
        new_root = ET.Element(xml_root.tag, xml_root.attrib)
        programme_element = new_root.find('.//programme')
        programme_element = ET.SubElement(new_root, 'programme')
        exposures_element = ET.SubElement(programme_element, 'exposures')

        obs_mode = field.attrib.get('obs_mode')
        ra_d = field.attrib.get('RA_d')
        dec_d = field.attrib.get('Dec_d')
        src_id = field.attrib.get('targid')
        exposures_element.append(exposure)
        exposure.set('order', '1')

        """
        slice_calib_df_0 = calib_df.iloc[0:4]
        if file_counter == 1:
            for idx, row in slice_calib_df_0.iterrows():
                calib_exposure = ET.SubElement(exposures_element, 'exposure')
                calib_exposure.set('arm', row['arm'])
                calib_exposure.set('cal_lamp', row['cal_lamp'])
                calib_exposure.set('cal_lamp_filter_A', row['cal_lamp_filter_A'])
                calib_exposure.set('cal_lamp_filter_B', row['cal_lamp_filter_B'])
                calib_exposure.set('ff_ilu_1',str(row['ff_ilu_1'])) 
                calib_exposure.set('ff_ilu_2', str(row['ff_ilu_2'])) 
                calib_exposure.set('ff_ilu_3', str(row['ff_ilu_3'])) 
                calib_exposure.set('order', str(row['order'])) 
                calib_exposure.set('speed', str(row['speed'])) 
                calib_exposure.set('type', str(row['type'])) 
                calib_exposure.set('exp_time', str(row['exp_time'])) 

        exposures_element.append(exposure)
        if file_counter == 1:
            exposure.set('order', '3')
        elif file_counter >= 2:
            exposure.set('order', '1')

        slice_calib_df_1 = calib_df.iloc[4:]
        if file_counter > 2:
            for idx, row in slice_calib_df_1.iterrows():
                calib_exposure = ET.SubElement(exposures_element, 'exposure')
                calib_exposure.set('arm', row['arm'])
                calib_exposure.set('cal_lamp', row['cal_lamp'])
                calib_exposure.set('cal_lamp_filter_A', row['cal_lamp_filter_A'])
                calib_exposure.set('cal_lamp_filter_B', row['cal_lamp_filter_B'])
                calib_exposure.set('ff_ilu_1',str(row['ff_ilu_1'])) 
                calib_exposure.set('ff_ilu_2', str(row['ff_ilu_2'])) 
                calib_exposure.set('ff_ilu_3', str(row['ff_ilu_3'])) 
                calib_exposure.set('order', str(row['order'])) 
                calib_exposure.set('speed', str(row['speed'])) 
                calib_exposure.set('type', str(row['type']))
                calib_exposure.set('exp_time', str(row['exp_time'])) 
        """

        observation_element = ET.SubElement(new_root, 'observation')
        observation_element.set('name', name)
        observation_element.set('pa', pa)
        observation_element.set('program', program)
        observation_element.set('obs_mode', str(obs_mode))
        observation_element.set('RA_d', str(ra_d))
        observation_element.set('Dec_d', str(dec_d))

        configure_element = ET.SubElement(observation_element, 'configure',
                                          attrib=configure_attributes)
        dithering_element = ET.SubElement(configure_element, 'dithering',
                                          attrib=dithering_attrributes)
        telescope_element = ET.SubElement(configure_element, 'telescope',
                                          attrib=telescope_attributes)
        focal_plane_element = ET.SubElement(configure_element, 'focal_plane_map',
                                            attrib=focal_plane_map_attributes)
        rotator_axis_element = ET.SubElement(focal_plane_element, 'rotator_axis',
                                             attrib=rotator_axis_attributes)
        optical_axis_element = ET.SubElement(focal_plane_element, 'optical_axis',
                                             attrib=optical_axis_attributes)
        distortion_coefficients_element = ET.SubElement(focal_plane_element,
                                                        'distortion_coefficients',
                                                attrib=distortion_coefficients_attributes)
        rotator_centre_element = ET.SubElement(focal_plane_element, 'rotator_centre',
                                               attrib=rotator_centre_attributes)
        hour_angle_limits_element  = ET.SubElement(configure_element, 'hour_angle_limits',
                                                   attrib=hour_angle_limits_attributes)

        conditions_element = ET.SubElement(configure_element, 'conditions',
                                           attrib=conditions_attributes)

        surveys_element = ET.SubElement(observation_element, 'surveys')
        for survey in surveys:
            survey_attributes = survey.attrib
            survey_element = ET.SubElement(surveys_element, 'survey', attrib=survey_attributes)

        fields_element = ET.SubElement(observation_element, 'fields')
        field_element = ET.SubElement(fields_element, 'field', attrib={'order': ''})
        field_element.set('RA_d', str(ra_d))
        field_element.set('Dec_d', str(dec_d))
        field_element.set('targid', str(src_id))

        targets = field.findall('target')

        for target in targets:
            field_element.append(target)

        new_tree = ET.ElementTree(new_root)
        output_file_path = os.path.join(output_dir, f'{original_filename}_split_{file_counter:03d}.xml')        
        save_xml(new_tree, output_file_path)

        file_counter += 1

    print(f"Generated {file_counter - 1} XML files with individual <field> elements")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split XML file into multiple field xmls.")
    parser.add_argument("-i", "--input", required=True, help="Input XML file path")
    parser.add_argument("-o", "--output", required=True, help="Output directory for the \
                        split XML files")

    args = parser.parse_args()

    xml_tree = ET.parse(args.input)
    xml_root = xml_tree.getroot()

    original_filename = os.path.splitext(os.path.basename(args.input))[0]
    split_xml_by_field(xml_root, args.output, original_filename)

    # usage: python xml_splitter.py -i /path/to/input.xml -o /path/to/output/directory/
