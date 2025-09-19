from xml.etree.ElementTree import ElementTree as ET
import xml.etree.cElementTree as ElementTree
import logging


def horiz_fix(xml):
    """Some earlier versions of Configure complained about targets being below
    the horizon limit. This was originally fixed temporarily by setting the field 
    RA and Dec to the C14 RA and Dec.

    This function is only called when the end user specifies a version of configure
    that is not the latest (and therefore fixed)
    
    """

    tree = ElementTree.parse(xml)
    root = tree.getroot()
    observation = root.find('.//observation')
    if observation.get('obs_mode') == 'LIFU':
        field = root.find('.//field')
        c14 = field.find(".//target[@targuse='T']")
        field.set("RA_d",c14.get("targra"))
        field.set("Dec_d",c14.get("targdec"))
        xml_tree = ET(root)
        file = open(xml,'w')
        logging.info("Fixing horizon check issue on {}".format(xml))
        xml_tree.write(file,encoding='unicode')
        file.close()


def set_pa(xml,pa):
    """Might be nice (for testing, certainly) to set the PA in advance before running
    configure. Note - this can *only* be for mIFU because changing the PA under LIFU
    means we lost the guidestar.
    
    """
    tree = ElementTree.parse(xml)
    root = tree.getroot()
    observation = root.find('.//observation')
    if observation.get('obs_mode') == 'mIFU':

        observation.set('pa',str(pa))
        xml_tree = ET(root)
        file = open(xml,'w')
        logging.info("Setting PA to user-requested value on {}".format(xml))
        xml_tree.write(file,encoding='unicode')
        file.close()


        
def auto_mifu(xml,src):
    """Following discussions with DT, this is a 'short-cut' function that uses a pre-defined
    field to assign bundles in the input XML. Configure then requires a "Re-map" under the GUI
    and the result can then be saved.
    
    This function assigns the fibreids for the target and calib bundles, the guide fibres and also
    ports the sky bundles.
    
    """
    stree = ElementTree.parse(src)
    sroot = stree.getroot()
    sfield = sroot.find('.//field')

    dtree = ElementTree.parse(xml)
    droot = dtree.getroot()
    dfield = droot.find('.//field')

    observation = droot.find('.//observation')
    if observation.get('obs_mode') == 'LIFU':
        return

    sobservation = sroot.find('.//observation')
    spa = sobservation.get('pa')
    observation.set('pa',spa)
    
    sfield = sroot.find('.//field')
    stargets = sfield.findall('.//target')
    
    bundle_lookup = {}
    guide_lookup = {}
    guides = []
    skies = {}
    central_lookup = {}
    for i in range(1,21):
        central = sfield.find(".//target[@ifu_spaxel='m%sC04']"%(str(i).zfill(2)))
        if central == None:
            continue
        try:
            central_lookup[central.get('targid')].append(central)
        except:
            central_lookup[central.get('targid')] = [central]
        targid = central.get('targid')
        fibreid = central.get('fibreid')
        if not fibreid:
            continue
        assert fibreid != None
        assert targid != None
        if central.get('targuse') == 'S':
            if central.get('fibreid') != None:
                skies[central.get('fibreid')] = central
        else:
            bundle_lookup[targid] = fibreid

    for target in sfield.findall(".//target[@targuse='G']"):
        if target.get('fibreid'):
            guide_lookup[target.get('cname')] = target.get('fibreid')

    # now dump into target XML
    dfield = droot.find('.//field')
    dtargets = dfield.findall('.//target')
    done = []
    done_lookup = {}
    
    for target in dtargets:
        targid = target.get('targid')
        try:
            _fibid = bundle_lookup[targid]
            if not targid in done:
                done.append(targid)
                done_lookup[targid] = target
        except:
            continue
        target.set('fibreid',_fibid)
    assert len(done) == len(bundle_lookup), "Of {} mIFU bundles to migrate, {} were done: is this the right field specified in --src_cfg?".format(len(done),len(bundle_lookup))
        
    for target in dfield.findall(".//target[@targuse='G']"):
        cname = target.get('cname')
        try:
            _fibid = guide_lookup[cname]
            guides.append(cname)
        except:
            continue
        target.set('fibreid',_fibid)


    for sky_id in skies:
        skytarg = skies[sky_id]
        del skytarg.attrib['configid']
        dfield.append(skytarg)

    xml_tree = ET(droot)
    file = open(xml,'w')
    logging.info("mIFU bundle allocations transferred from {} to {}".format(src,xml))
    xml_tree.write(file,encoding='unicode')
    file.close()


    
