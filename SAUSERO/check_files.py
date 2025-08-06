"""
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.

Copyright (C) 2025 Gran Telescopio Canarias <https://www.gtc.iac.es>
Fabricio Manuel Pérez Toledo <fabricio.perez@gtc.iac.es>
"""

import json, os
import pkg_resources
import ccdproc as ccdp
from pathlib import Path

def read_config(config_path):
    with open(config_path, 'r') as file:
        config = json.load(file)
    return config

def readJSON():
    """
    Reads the file containing the configuration parameters.

    Returns:
        json: Collection of configuration parameters 
    """
    if os.path.exists(Path(os.getcwd())/'configuration.json'):
        return json.load(open(Path(os.getcwd())/'configuration.json'))
    else:
        config_path = pkg_resources.resource_filename(
            'SAUSERO', 'config/configuration.json')
        return json.load(open(config_path))


def update_config(config_path, config):
    with open(config_path, 'w') as file:
        json.dump(config, file, indent=4)
    
    return config


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

    if 'OsirisBroadBandImage' in tab['OBSMODE'] and ('Clear' in tab['FILTER2'] or 'CLEAR' in tab['FILTER2']):
        existence['exist_STD'] = False
        existence['exist_SCIENCE'] = True
    
    return existence


def check_files():

    conf = readJSON()

    directory = Path(os.getcwd())/'raw'

    ic = ccdp.ImageFileCollection(directory, keywords=['GTCPRGID','GTCOBID','OBSMODE','OBJECT','FILTER2','EXPTIME'])
    image_types = classify_images(ic.summary)

    conf['PRG'] = [elem for elem in set(ic.summary['GTCPRGID'].value.data) if not 'CALIB' in elem][0]
    conf['OB'] = [elem for elem in set(ic.summary['GTCOBID'].value.data) if not 'CALIB' in elem][0]

    conf['DIRECTORIES']['PATH'] = str(Path(os.getcwd()))
    conf['DIRECTORIES']['PATH_DATA'] = str(Path(os.getcwd())/'raw')
    conf['DIRECTORIES']['PATH_OUTPUT'] = str(Path(os.getcwd())/'reduced') 

    # Update config based on image types
    if conf['REDUCTION']['use_BIAS']:
        conf['REDUCTION']['use_BIAS'] = image_types['exist_BIAS']
    
    if conf['REDUCTION']['use_FLAT']:
        conf['REDUCTION']['use_FLAT'] = image_types['exist_SKYFLAT']

    if conf['REDUCTION']['use_STD']:
        conf['REDUCTION']['use_STD'] = (image_types['exist_STD'] and image_types['exist_SKYFLAT'] and conf['REDUCTION']['use_FLAT'])

    if conf['REDUCTION']['save_std']:
        conf['REDUCTION']['save_std'] = (image_types['exist_STD'] and image_types['exist_SKYFLAT'] and conf['REDUCTION']['use_STD'])
    
    if conf['REDUCTION']['save_sky']:
        conf['REDUCTION']['save_sky'] = image_types['exist_SCIENCE']
    
    if conf['PHOTOMETRY']['use_photometry']:
        conf['PHOTOMETRY']['use_photometry'] = (image_types['exist_SKYFLAT'] and image_types['exist_STD'] and conf['REDUCTION']['use_STD'])

    return update_config(conf['DIRECTORIES']['PATH'] + '/configuration.json', conf)