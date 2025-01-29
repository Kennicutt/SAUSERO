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

import yaml, glob, os

from astropy.nddata import CCDData
from astropy.wcs import WCS
from astrometry_net_client import Session, FileUpload, Settings

from SAUSERO.Color_Codes import bcolors as bcl
from loguru import logger

def settings(PATH_TO_CONFIG_FILE):
    """ This method opens the configuration file containing the parameters 
    needed for optimal astrometrization.

    Args:
        PATH_TO_CONFIG_FILE (str): Path to configuration file.

    Returns:
        dict: Contents of the configuration file.
    """
    try:
        print('Loading settings for astrometrization...')
        with open(PATH_TO_CONFIG_FILE, 'r') as file:
            prime_service = yaml.safe_load(file)
    except FileNotFoundError:
        raise('Not found the configuration file.')

    return prime_service['OSIRIS']

def apply_astrometrynet_client(filename, conf):
    """This method sends the stacked science frame along with settings information 
    to the Astrometry.net server through its API. Then, it receives the 
    WCS estimated by the server.

    Args:
        filename (str): Name of stacked science frame
        conf (dict): Collection of setting parameters

    Returns:
        WCS: WCS information given by the server.
    """    
    logger.info("Loading the settings")
    ss = Settings()
    img = CCDData.read(filename, unit='adu')
    #ss.scale_units = 'arcsecperpix'
    ss.set_scale_estimate(conf["ASTROMETRY"]["set_scale_estimate"]["scale"],
                        conf["ASTROMETRY"]["set_scale_estimate"]["unknown"], 
                        unit=conf["ASTROMETRY"]["set_scale_estimate"]["scale_units"])
    ss.center_ra = img.header['RADEG']
    ss.center_dec = img.header['DECDEG']
    ss.radius = conf["ASTROMETRY"]["radius"]
    ss.downsample_factor = conf["ASTROMETRY"]["downsample_factor"]
    ss.use_sextractor = conf["ASTROMETRY"]["use_sextractor"]
    ss.crpix_center = conf["ASTROMETRY"]["crpix_center"]
    ss.parity = conf["ASTROMETRY"]["parity"]
    ss.allow_commercial_use = 'n'
    ss.allow_modifications = 'n'
    ss.publicly_visible = 'n'
    logger.info("Settings:")
    logger.info(f"{ss}")
    #Send the image
    s = Session(api_key=conf["ASTROMETRY"]["No_Session"])
    logger.info("API connection is ready")
    upl = FileUpload(filename, session=s, settings=ss)
    logger.info("Frame has been uploaded")
    submission = upl.submit()
    logger.info("Waiting for an answer from API...")
    submission.until_done()
    job = submission.jobs[0]
    job.until_done()
    if job.success():
        wcs = job.wcs_file()
        logger.info("WCS received from API")
    logger.info(job.info())
    
    if wcs != None:
        return wcs

def modify_WCS(best_wcs, PATH_TO_FILE):
    """This method updates the original WCS with the WCS estimated by the Astrometry.net server.

    Args:
        best_wcs (str): WCS estimated by the Astrometry.net server
        PATH_TO_FILE (str): Path to science frame

    Returns:
        CCDData: Science frame with its WCS updated.
    """
    frame = CCDData.read(PATH_TO_FILE, unit='adu')
    best_wcs = WCS(best_wcs)
    new_frame =  CCDData(data=frame.data, header=frame.header, wcs=best_wcs,
                         unit='adu')
    new_frame.write(PATH_TO_FILE, overwrite=True)
    logger.info(f"The WCS for {os.path.basename(PATH_TO_FILE)} has been updated")
    return new_frame

def solving_astrometry(PRG, OB, filt, conf, sky, calib_std = False):
    """This method handles the procedure for obtaining the astrometry solution from the server.

    Args:
        PRG (str): Science program code
        OB (str): Observational block number assigned to a science program
        filt (str): Filter name

    Returns:
        WCS, list, CCDData: Several results: the new WCS, a local catalogue for 
        the FoV using the new WCS to estimate positions, and the science frame 
        with its WCS updated.
    """
    root_path = conf["DIRECTORIES"]["PATH_DATA"]

    if calib_std:
        logger.info("Astrometry calibration for STD star")
        LST_PATH_TO_FILE = glob.glob(root_path + PRG + '_' + OB + '/reduced/' + f'*STD*.fits')
    else:
        logger.info("Astrometry calibration for science target")
        LST_PATH_TO_FILE = [root_path + PRG + '_' + OB + '/reduced/' + f'{PRG}_{OB}_{filt}_stacked_{sky}.fits']

    PATH_TO_FILE = LST_PATH_TO_FILE[0]

    best_wcs = apply_astrometrynet_client(PATH_TO_FILE, conf)

    new_frame = modify_WCS(best_wcs, PATH_TO_FILE)

    return best_wcs, new_frame #pixels
