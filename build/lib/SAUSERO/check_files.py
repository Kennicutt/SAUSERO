import json
import os
import ccdproc as ccdp

def read_config(config_path):
    with open(config_path, 'r') as file:
        config = json.load(file)
    return config


def update_config(config_path, config):
    with open(config_path, 'w') as file:
        json.dump(config, file, indent=4)


def classify_images(tab):

    existence = {
        'exist_BIAS': False,
        'exist_SKYFLAT': False,
        'exist_SCIENCE': False,
        'exist_STD': False
    }

    if 'OsirisBias' in tab['OBSMODE']:
        existence['exist_BIAS'] = True
    
    if 'OsirisSkyFlat' in tab['OBSMODE']:
        existence['exist_SKYFLAT'] = True

    if 'OsirisBroadBandImage' in tab['OBSMODE']:
        existence['exist_STD'] = (len([item for item in tab['OBJECT'] if item.startswith('STD')]) != 0)
        existence['exist_SCIENCE'] = (len([item for item in tab['OBJECT'] if not item.startswith('STD')]) != 0)

    #Special case for OsirisBroadBandImage

    if ('OsirisBroadBandImage' in tab['OBSMODE']) and ('OPEN' in tab['FILTER2']):
        existence['exist_STD'] = False
        existence['exist_SCIENCE'] = True
    
    return existence


def check_files(config_path, PRG, OB):

    conf = read_config(config_path)
    directory = conf['DIRECTORIES']['PATH_DATA'] + f"{PRG}_{OB}/raw/"

    ic = ccdp.ImageFileCollection(directory, keywords=['OBSMODE','OBJECT','FILTER2','EXPTIME'])
    image_types = classify_images(ic.summary)

    # Update config based on image types
    conf['REDUCTION']['use_BIAS'] = image_types['exist_BIAS']
    conf['REDUCTION']['use_FLAT'] = image_types['exist_SKYFLAT']
    conf['REDUCTION']['use_STD'] = (image_types['exist_STD'] and image_types['exist_SKYFLAT'])
    conf['REDUCTION']['save_std'] = (image_types['exist_STD'] and image_types['exist_SKYFLAT'])
    conf['REDUCTION']['save_sky'] = image_types['exist_SCIENCE']
    conf['PHOTOMETRY']['use_photometry'] = (image_types['exist_SKYFLAT'] and image_types['exist_STD'])

    update_config(config_path, conf)