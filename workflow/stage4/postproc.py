from xml.etree.ElementTree import ElementTree as ET
import xml.etree.cElementTree as ElementTree
import logging
from copy import deepcopy

def filter_spares(xml):
    """Configure includes the 'spare' (as defined in the LIFU.xml status) fibres as
    <target> elements. In a typical case of what we want at WASP vs what we want in
    FIBTABLE for CPS, these cause problems downstream with workflow + WASP code, so 
    easier to remove them from the configure output.


    <fibre status="spare" pos_x="0.000" slit_entry="1" slitid="21" fibre_name="" id="21" type="Spectrograph" pos_y="-20.000"/>
    <fibre status="spare" pos_x="0.000" slit_entry="16" slitid="316" fibre_name="" id="316" type="Spectrograph" pos_y="-20.000"/>
    <fibre status="spare" pos_x="0.000" slit_entry="19" slitid="393" fibre_name="" id="393" type="Spectrograph" pos_y="-20.000"/>
    <fibre status="spare" pos_x="0.000" slit_entry="20" slitid="404" fibre_name="" id="404" type="Spectrograph" pos_y="-20.000"/>
    <fibre status="spare" pos_x="0.000" slit_entry="22" slitid="451" fibre_name="" id="451" type="Spectrograph" pos_y="-20.000"/>
    <fibre status="spare" pos_x="0.000" slit_entry="24" slitid="486" fibre_name="" id="486" type="Spectrograph" pos_y="-20.000"/>

    ...this code identifies them as <target> elements with -ve fibreID and a blank ifu_spaxel ID

    The question is, will these targets "come back" and be written into the FIBTABLE, even though they are not
    in the XML? I *think* so, as the POS output will still contain them (just they will have no counterparts in the
    XML component stored at the end of the POS output).

    """

    tree = ElementTree.parse(xml)
    root = tree.getroot()
    fields = root.findall('.//field')
    for field in fields:
        targets = field.findall('.//target')
        for target in targets:
            if (target.get('fibreid')) and int(target.get('fibreid')) < 0:
                if target.get('ifu_spaxel') == '':
                    # this is a spare fibre - remove it
                    field.remove(target) 
    xml_tree = ET(root)
    file = open(xml,'w')
    logging.info("Writing filtered XML to {}".format(xml))
    xml_tree.write(file,encoding='unicode')
    file.close()

def add_offsets(xml):
    """
    Pushing an OB through dither that has no dither requirement (ie apply_dither=0)
    *does* add some of the dither-related attributes (dither_version for example), but
    not the offsets elelemt. This function adds in this element and sets blanks attribs
    for it.

    """

    
    tree = ElementTree.parse(xml)
    root = tree.getroot()
    observation = root.find('.//observation')
    dithering = root.find('.//dithering')
    if dithering.get('apply_dither') == '0':
        offsets = ElementTree.Element("offsets")
        offsets.set("offset_step_ra","")
        offsets.set("offset_step_dec","")
        observation.insert(0,offsets)

        xml_tree = ET(root)
        file = open(xml,'w')
        logging.info("Added offsets element to {}".format(xml))
        xml_tree.write(file,encoding='unicode')
        file.close()
    

def guide_pa_unset(xml):
    """Per the BlankXMLTemplate.xml template file, guidestar target elements
    have a blank ifu_pa value - they should not inherit from the <observation pa>.
    This function just blanks out the ifu_pa for both LIFU and mIFU guidestar
    target elements.

    """
    
    tree = ElementTree.parse(xml)
    root = tree.getroot()
    observation = root.find('.//observation')
    fields = root.findall('.//field')
    for field in fields:
        guides = field.findall(".//target[@targuse='G']")
        for guide in guides:
            guide.set('ifu_pa','')
    xml_tree = ET(root)
    file = open(xml,'w')
    logging.info("Blanking out guidestar ifu_pa attributes on {}".format(xml))
    xml_tree.write(file,encoding='unicode')
    file.close()

def ifu_pa_set(xml):
    """Minor fix to ifu_pa precision.

    Applies to 3.3.11, after LIFU field configuration. 
    The precision of ifu_pa in the <target> elements is
    only 1dp, causing cross-matching issues downstream in some cases.

    Temporary fix - set all ifu_pa to the value of <observation pa>

    """
    
    tree = ElementTree.parse(xml)
    root = tree.getroot()
    observation = root.find('.//observation')
    if observation.get('obs_mode') != 'LIFU':
        return
    pa = observation.get('pa')
    fields = root.findall('.//field')
    for field in fields:
        targets = field.findall(".//target")
        for target in targets:
            if (target.get('ifu_pa') != None) and (target.get('ifu_pa') != ''):
                target.set('ifu_pa',pa)
    xml_tree = ET(root)
    file = open(xml,'w')
    logging.info("Reverting to ifu_pa atributes to original precision on {}".format(xml))
    xml_tree.write(file,encoding='unicode')
    file.close()


    
def reset_random_cnames(xml):
    """An annoying one: when populating a mIFU (LIFU? TBD) bundle with a WD central,
    configure has the surrounding (targuse='R') targets inherit the central's CNAME.
    For various reasons, WASP instead wants this set as '%%%', and will calculate a new
    CNAME for it. 

    This function just sets the <target cname> attribute to '%%%' for these (R)andom
    fibres. There's no easy fix under configure for this, beyond a "set this attribute
    to %%%" directive - not nice (ie "inherit, but now change this thing to another thing").

    Thinking of a way to kludge this downstream in the workflow.

    """
    
    tree = ElementTree.parse(xml)
    root = tree.getroot()
    observation = root.find('.//observation')
    fields = root.findall('.//field')
    for field in fields:
        randoms = field.findall(".//target[@targuse='R']")
        for random in randoms:
            if random.get('targcat') == 'WD.fits':
                random.set('cname','%%%')
    xml_tree = ET(root)
    file = open(xml,'w')
    logging.info("Reset CNAMES for targuse='R' on {}".format(xml))
    xml_tree.write(file,encoding='unicode')
    file.close()

    
def populate_random_photometry(xml):
    """
    dither does not inherit the photometry sub-element for mIFU calibration bundles (ie, both central 
    targuse 'C' and surrounding 'R' fibres) copied over to the next field(s) (into targuse 'R' fibres).

    As (at this stage) any photometry data is just duplicated from the source central fibre, this function
    just copies this photometry element into successive fields. 

    """

    tree = ElementTree.parse(xml)
    root = tree.getroot()
    observation = root.find('.//observation')
    if observation.get('obs_mode') != 'mIFU':
        return
    if len(root.findall('.//field')) == 1:
        # not dithered - nothing to do here
        return
    
    calib_bundle_lookup = {}
    field = root.find('.//field')
    calibs = field.findall(".//target[@targuse='C']")
    for calib in calibs:
        bundle = calib.get('ifu_spaxel')
        if not bundle:
            # not allocated a bundle - skip
            continue
        calib_bundle_lookup[bundle[:3]] = calib.find(".//photometry")

    fields = root.findall('.//field')[1:]
    for field in fields:
        randoms = field.findall(".//target[@targuse='R']")
        for random in randoms:
            _bundle = calib.get('ifu_spaxel')
            if not _bundle:
                continue
            try:
                _photometry = calib_bundle_lookup[_bundle[:3]]
            except:
                continue
            photclone = deepcopy(_photometry)
            random.append(photclone)
        
    xml_tree = ET(root)
    file = open(xml,'w')
    logging.info("Blanking out guidestar ifu_pa attributes on {}".format(xml))
    xml_tree.write(file,encoding='unicode')
    file.close()

    
    
def populate_sky(xml):
    """In addition to existing inherited attributes, LIFU sky fibres should
    also inherit the following attributes from the central fibre:

    'targsrvy','targprog','targcat','targname','targid'

    This function simply cycles through all the sky fibres and if their attributes are 
    blank, it will set them to the value from the originally-defined <target> used to spawn
    all fibres for LIFU.

    """

    
    tree = ElementTree.parse(xml)
    root = tree.getroot()
    observation = root.find('.//observation')
    if observation.get('obs_mode') == 'LIFU':
        fields = root.findall('.//field')
        inherit = ['targsrvy','targprog','targcat','targname','targid']

        for field in fields:
            c14 = field.find(".//target[@ifu_spaxel='C14']")
            skies = field.findall(".//target[@targuse='S']")
            for target in skies:
                if (target.get('targuse') == 'S'):
                    for attr in inherit:
                        if target.get(attr) == '':
                            target.set(attr,c14.get(attr))
                            
        xml_tree = ET(root)
        file = open(xml,'w')
        logging.info("Writing sky-fixed XML to {}".format(xml))
        xml_tree.write(file,encoding='unicode')
        file.close()
    

def reset_targpos(xml,src):
    """
    Dither post-proc function
    
    If an OB has been dithered, the fibre targra and targdec will have been 
    modified slightly, even if (in the case of dither != -3) the initial dither 
    step is 0.0. This function temporarily corrects the issue by setting the targra 
    and targdec of all fibres in the first field to that of the input XML.

    match by ifu_spaxel in the case of mIFU, configid for LIFU

    """
    
    stree = ElementTree.parse(src)
    sroot = stree.getroot()

    observation = sroot.find('.//observation')
    identifier = 'ifu_spaxel'
    if observation.get('obs_mode') == 'LIFU':
        identifier = 'configid'

    dithering = sroot.find('.//dithering')
    if dithering.get('apply_dither') in ['0','-3']:
        return

    sfield = sroot.find('.//field')
    sfield = sroot.find('.//field')
    stargets = sfield.findall('.//target')

    targ_lookup = {}
    for target in stargets:
        fib_indentifier = target.get(identifier)
        if (fib_indentifier) and (fib_indentifier != ''):
            targ_lookup[fib_indentifier] = {'targra':target.get('targra'),'targdec':target.get('targdec')}

    # now dump into target XML
    dtree = ElementTree.parse(xml)
    droot = dtree.getroot()
    dfield = droot.find('.//field')

    dfield = droot.find('.//field')
    dtargets = dfield.findall('.//target')

    for target in dtargets:
        fib_indentifier = target.get(identifier)
        if (fib_indentifier) and (fib_indentifier != ''):
            tpos = targ_lookup[fib_indentifier]
            target.set('targra',tpos['targra'])
            target.set('targdec',tpos['targdec'])
            
    xml_tree = ET(droot)
    file = open(xml,'w')
    logging.info("{} initial field positions reset on {}".format(observation.get('obs_mode'),xml))
    xml_tree.write(file,encoding='unicode')
    file.close()


        

    
