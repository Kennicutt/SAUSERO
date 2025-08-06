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

import os, sys, time, json, logging
from pathlib import Path

from astropy import units as u
from astropy import wcs
from astropy.nddata import CCDData
from astropy.io import fits
import ccdproc as ccdp
from matplotlib import pyplot as plt
import numpy as np
import yaml as py
from lacosmic.core import lacosmic
import sep

from SAUSERO.Color_Codes import bcolors as bcl
from loguru import logger

import logging, inspect
class InterceptHandler(logging.Handler):
    """Custom logging handler that redirects log records from the standard logging module 
    to the Loguru logger.

    This handler allows the use of Loguru alongside the standard logging module, ensuring 
    that logs captured by the standard logging system are properly forwarded to the Loguru logger.

    Args:
        logging (module): The logging module that provides the log records to be captured 
                           and forwarded to Loguru. Typically, this would be the standard 
                           `logging` module.
    """
    
    def emit(self, record: logging.LogRecord) -> None:
        """
        Process a log record and forward it to the Loguru logger, mapping the log level and 
        determining the caller's origin for the log message.

        This method retrieves the log level, the caller's information, and the log message 
        before passing it to Loguru for proper logging.

        Args:
            record (logging.LogRecord): The log record to be processed, which contains 
                                        information about the log level, message, and other 
                                        relevant details for the log entry.

        This method will convert the logging level to a corresponding Loguru level and 
        propagate the message along with the exception information (if any) to the Loguru logger.
        """
        # Get corresponding Loguru level if it exists.
        level: str | int
        try:
            level = logger.level(record.levelname).name
        except ValueError:
            level = record.levelno

        # Find caller from where originated the logged message.
        frame, depth = inspect.currentframe(), 0
        while frame and (depth == 0 or frame.f_code.co_filename == logging.__file__):
            frame = frame.f_back
            depth += 1

        if level != "DEBUG":
            logger.opt(depth=depth, exception=record.exc_info).log(level, record.getMessage())

logging.basicConfig(handlers=[InterceptHandler()], level=logging.INFO, force=False)

class Reduction:
    """The goal is to perform the cleaning procedure for science and photometric calibration frames. 
    First, bias frames are averaged to create a master bias. This master bias is then used to subtract 
    the offset level from the sky flat. The sky flats are subsequently combined and divided by their 
    median to normalize them. The final product is a master flat, which will be used to standardize 
    the sensor's response. Finally, the science and photometric calibration frames are corrected by 
    subtracting the master bias and dividing by the master flat, resulting in the final reduced frames.
    """

    def __init__(self, main_path, path_mask = None):
        """Object initialization to carry out the reduction.

        Args:
            gtcprgid (str): Observation program code.
            gtcobid (str): Observation block code.
            path_mask (str, optional): Path to BPM.
        """

        self.PATH = Path(main_path)
        
        if os.path.exists(self.PATH):
            logger.info("Path to raw data exists")
        else:
            logger.critical(f"{bcl.ERROR}Path to raw data does NOT exist{bcl.ENDC}")
            sys.exit()

        # The directory that will contain the intermediate and processed images is defined
        self.PATH_RESULTS = self.PATH.parent/'reduced'

        #self.PATH_RESULTS.mkdir(parents=True, exist_ok=True)
        if os.path.exists(self.PATH_RESULTS):
            logger.info("Directory to reduced files exists")
        else:
            logger.critical(f"Directory to reduced files doesn't exists")
            self.PATH_RESULTS.mkdir(parents=True, exist_ok=True)
            logger.info(f"{bcl.OKGREEN}Directory to reduced files has been created{bcl.ENDC}")

        # Define the path to mask file
        self.path_mask = Path(path_mask)
        if os.path.exists(self.path_mask):
            logger.info("Path to mask file exists")
        else:
            logger.critical(f"{bcl.ERROR}Path to mask file does NOT exist{bcl.ENDC}")
            sys.exit()
        
        # The information about the frames in that directory is gathered.
        self.ic = ccdp.ImageFileCollection(self.PATH)
        if len(self.ic.summary) != 0:
            logger.info("Data collection is ready")
        else:
            logger.critical("NOT files to reduce")
            sys.exit()

        # Define dictionaries
        self.DATA_DICT={}
        self.key_dict={'flat':'OsirisSkyFlat',
                       'target': 'OsirisBroadBandImage',
                       'std': 'OsirisBroadBandImage',
                       'bias': 'OsirisBias'
        }
        self.master_dict = {} #Ready to use
        self.std_dict = {}
        self.target_dict = {}

        logger.info("Dictionaries have been created")




    @staticmethod
    def configure_mask(mask):
        """Static method to reshape the BPM array to match any frame shape.

        Args:
            mask (int): BPM HDU or CCDData Object

        Returns:
            bool: BPM array in boolean format to apply over the frames.
        """
        mask = mask.data.astype(bool)
        matrix = np.ones(mask.shape)
        matrix[mask == False] = np.nan
        return matrix[230:2026,28:2060] # TRIM SECTION




    @staticmethod
    def get_each_data(data_dict, value):
        """Reads the contents of a dictionary containing the paths to a series of 
        frames and opens them, adding them to a list.

        Args:
            data_dict (dict): A dictionary containing the paths to 
            one or more frames in a list.
            value (str): Key to access the list of paths.

        Returns:
            list: It is a list containing the images.
        """
        ccd = []
        for frame_path in data_dict[value]:
            hdul = fits.open(frame_path)
            ccd.append(hdul[0].data[230:2026,28:2060]) #TRIM SECTION
            
        logger.info(f"List of images for key {value} is ready.")
        return ccd



    @staticmethod
    def combining(lst_frames):
        """This static method combines the images in a list to obtain an averaged image. 
        The process involves creating a data cube and averaging the images.

        Args:
            lst_frames (list): List of images to be averaged.

        Returns:
            image: Create an averaged matrix from a data cube.
        """
        cube = np.dstack(lst_frames)
        cube.sort(axis=2)
        return np.nanmedian(cube[:,:,1:-1], axis=2)




    def create_cubes(self, key):
        """Create a data cube from a list of images.

        Args:
            key (str): Key of the scientific images dictionary.

        Returns:
            cube(float): Data cube.
        """
        return np.dstack(self.target_dict[key],
                        axis=2)




    def sustractMasterBias(self, value, master, data_dict):
        """Subtract the master bias from the image.

        Args:
            value (str): Key that allows access to the images in 
            the dictionary.
            master (float): MasterBias frame.
            data_dict (dict): Images dictionary. 

        Returns:
            list: List of the frames with masterbias applied.
        """
        ccd = self.get_each_data(data_dict, value)

        frames = [fr - master for fr in ccd]

        logger.info(f"Masterbias applied to key {value} set.")
        
        return frames




    def clean_target(self, value, masterbias, masterflat):
        """Applies the subtraction of the MasterBias and the division 
        by the normalized MasterFlat to the science images.

        Returns:
            list: List of cleaned science frames.
        """
        #value, masterbias, masterflat = args
        ccd = self.get_each_data(self.DATA_DICT, value)
        if masterflat is None:
            masterflat = 1.0
        frames = [(fr - masterbias)/masterflat for fr in ccd]
        logger.info(f"Applied masterbias and masterflat on the frames for {value}.")
        return frames


    def get_imagetypes(self):
        """ This method aims to identify all the types of images present in the original images 
        directory and their corresponding filters. To achieve this, it creates an empty dictionary 
        where each key is a combination of the image type and the filter used.
        """
        logger.info('Getting types of images and filters used.')
        self.filt_wheels = []
        matches = (self.ic.summary['obsmode'] != 'OsirisBias')
        matches1 = (self.ic.summary['obsmode'] != 'OsirisBias') & (self.ic.summary['filter1'] != 'OPEN')
        matches2 = (self.ic.summary['obsmode'] != 'OsirisBias') & (self.ic.summary['filter2'] != 'OPEN')
        matches3 = (self.ic.summary['obsmode'] != 'OsirisBias') & (self.ic.summary['filter3'] != 'OPEN')
        matches4 = (self.ic.summary['obsmode'] != 'OsirisBias') & (self.ic.summary['filter3'] == 'OPEN')\
        & (self.ic.summary['filter2'] == 'OPEN') & (self.ic.summary['filter1'] == 'OPEN')\
        & (self.ic.summary['filter4'] == 'OPEN')
        

        if len(list(set(self.ic.summary['filter1'][matches1]))) >= 1:
            self.filt_wheels.append('filter1')
        elif len(list(set(self.ic.summary['filter2'][matches2]))) >= 1:
            self.filt_wheels.append('filter2')
        elif len(list(set(self.ic.summary['filter3'][matches3]))) >= 1:
            self.filt_wheels.append('filter3')
        elif len(list(set(self.ic.summary['filter4'][matches4]))) >= 1:
            self.filt_wheels.append('filter4')
        else:
            raise(ValueError, "Incompatible filter setup!!!")
            sys.exit()

        for filt in self.filt_wheels:
            for value in list(set(self.ic.summary[filt][matches])):
                for type in ['flat', 'std', 'target']:
                    self.DATA_DICT[type + '+' + value] = []

        self.DATA_DICT["bias"] = []


    def load_BPM(self):
        """
        This method opens the FITS file containing the BPM.
        """
        bpm = CCDData.read(self.path_mask, unit=u.dimensionless_unscaled,
                            hdu=0, format='fits', ignore_missing_simple=True)
        self.MASK = self.configure_mask(bpm)
        logger.info("BPM is ready.")




    def load_results(self):
        """It generates a table with specific content in the results directory, specifically for 
        those scientific images that have already been cleaned.
        """
        self.ic_r = ccdp.ImageFileCollection(self.PATH_RESULTS, keywords='*',
                                             glob_include='red*')
        logger.info("Table with several reduced science frames is ready.")




    def sort_down_drawer(self):
        """This method is responsible for storing the images by their type and filter in the previously created dictionary.
        """
        for filt in self.filt_wheels:
            for elem in list(self.DATA_DICT.keys()):
                if elem != 'bias':
                    key, value = elem.split('+')
                    types_targets = set(self.ic.summary['object'][self.ic.summary['obsmode'] == 'OsirisBroadBandImage'])
                    if key == 'flat':
                        try:
                            tmp_dict = {"obsmode":self.key_dict[key], filt: value}
                            logger.info("Including flat frames.")
                        except Exception as e:
                            logger.warning(f"There are no flat frames in the directory: {self.PATH}")
                            self.DATA_DICT[elem] = []
                            continue
                    elif key == 'std':
                        try:
                            target_type = [data for data in types_targets if 'STD' in data][0]
                            tmp_dict = {"obsmode":self.key_dict[key], filt: value, 'object':target_type}
                            logger.info("Including photometric calibration frames.")
                        except Exception as e:
                            logger.warning(f"There are no photometric calibration frames in the directory: {self.PATH}")
                            self.DATA_DICT[elem] = []
                            continue
                    elif key == 'target':
                        try:
                            target_type = [data for data in types_targets if not 'STD' in data][0]
                            tmp_dict = {"obsmode":self.key_dict[key], filt: value, 'object':target_type}
                            logger.info("Including science frames.")
                        except Exception as e:
                            logger.warning(f"There are no science frames in the directory: {self.PATH}")
                            self.DATA_DICT[elem] = []
                            continue
                    
                    self.DATA_DICT[elem] = self.ic.files_filtered(**tmp_dict, include_path=True)
                else:
                    try:
                        self.DATA_DICT[elem] = self.ic.files_filtered(obsmode='OsirisBias', include_path=True)
                        logger.info("Including bias frames.")
                    except Exception as e:
                        logger.warning(f"There are no bias frames in the directory: {self.PATH}")
                        self.DATA_DICT[elem] = []
                        continue



    def do_masterbias(self):
        """
        This method creates the master bias frame.
        """
        bias = self.get_each_data(self.DATA_DICT, "bias")

        self.masterbias = self.combining(bias) * self.MASK
        self.master_dict['bias'] = self.masterbias

        logger.info(f"{bcl.OKGREEN}Masterbias has been created{bcl.ENDC}")



    def do_masterflat(self):
        """
        This method creates the master flat frame for each filter.
        """
        lst_flat = [elem for elem in list(self.DATA_DICT.keys()) if 'flat' in elem]
        for filt in lst_flat:
            if 'flat+OPEN' in filt:
                continue
            flat = self.sustractMasterBias(filt, self.masterbias, self.DATA_DICT)
            combflat = self.combining(flat)
            median = np.nanmedian(combflat)
            masterflat = combflat/median
            self.master_dict[filt] = masterflat
            logger.info(f"{bcl.OKGREEN}Masterflat has been created for {filt} filter{bcl.ENDC}")


    def get_std(self, no_CRs=False, contrast_arg = 1.5, cr_threshold_arg = 5.,
                neighbor_threshold_arg = 5., apply_flat=False):
        """
        This method processes the photometric calibration frames.
        """
        lst_std = [elem for elem in list(self.DATA_DICT.keys()) if 'std' in elem]
        logger.info("Processing photometric calibration frames.")
        for elem in lst_std:
            if 'std+OPEN' in elem:
                continue

            key, value = elem.split('+')
            
            if apply_flat == False:
                self.master_dict['flat+' + value] = None
        
            std = self.clean_target(elem, self.master_dict['bias'],self.master_dict['flat+' + value])
            lst_sd = []
            for sd in std:
                if no_CRs:
                    logger.info(f"Removing CRs to photometric calibration frame for {value}.")
                    no_mask = np.nan_to_num(sd, nan=np.nanmedian(sd))
                    lst_sd.append(lacosmic(no_mask, contrast=contrast_arg,
                                                    cr_threshold=cr_threshold_arg,
                                                    neighbor_threshold=neighbor_threshold_arg,
                                                    effective_gain=1.9,
                                                    readnoise=4.3)[0])
                else:
                    logger.info(f"NOT treatment for CRs applied to photometric calibration frame for {value}.")
                    lst_sd.append(np.nan_to_num(sd, nan=np.nanmedian(sd)))
            self.std_dict[elem] = lst_sd



    def get_target(self, no_CRs=False, contrast_arg = 1.5, cr_threshold_arg = 5.,
                neighbor_threshold_arg = 5., apply_flat=False):
        """
        This method cleans the science frames.
        """
        lst_target = [elem for elem in list(self.DATA_DICT.keys()) if 'target' in elem]
        logger.info("Processing science frames.")
        for elem in lst_target:
            key, value = elem.split('+')
            
            if (apply_flat == False) or (value == 'OPEN'):
                self.master_dict['flat+' + value] = None

            target= self.clean_target(elem, self.master_dict['bias'],self.master_dict['flat+' + value])
            lst_tg = []
            for tg in target:
                if no_CRs:
                    logger.info(f"Removing CRs to science frames for {value}.")
                    no_mask = np.nan_to_num(tg, nan=np.nanmedian(tg))
                    lst_tg.append(lacosmic(no_mask, contrast=contrast_arg,
                                                    cr_threshold=cr_threshold_arg,
                                                    neighbor_threshold=neighbor_threshold_arg,
                                                    effective_gain=1.9,
                                                    readnoise=4.3)[0])
                else:
                    logger.info(f"NOT treatment for CRs applied to science frames for {value}.")
                    lst_tg.append(np.nan_to_num(tg, nan=np.nanmedian(tg)))
            self.target_dict[elem] = lst_tg



    def remove_fringing(self):
        """
        This method performs a special cleaning when using Sloan_z to remove the interference pattern.
        """
        lst_results = [elem for elem in list(self.DATA_DICT.keys())]
        if 'target+Sloan_z' in lst_results:
            logger.info("Removing the fringe on Sloan z filter.")
            fringe= self.target_dict['target+Sloan_z']
            combfringe = self.combining(fringe)
            median = np.nanmedian(combfringe)
            masterfringe = combfringe/median
            fr_free = [elem/masterfringe for elem in fringe]
            self.target_dict['fringe+Sloan_z'] = fr_free
            logger.info("Sci frames with free fringe.")
            try:
                fringe_std = self.std_dict['std+Sloan_z']
                fr_free_std = [elem/masterfringe for elem in fringe_std]
                self.std_dict['fringe+Sloan_z'] = fr_free_std
                logger.info("STD frame/s with free fringe.")
            except Exception as e:
                logger.warning(f"There are no photometric calibration frames with Sloan z filter: {e}")
                self.std_dict['fringe+Sloan_z'] = []
        else:
            logger.warning(f"There are no science frames with Sloan z filter to remove the fringe")
            self.target_dict['fringe+Sloan_z'] = []
            self.std_dict['fringe+Sloan_z'] = []



        
    
    def sustract_sky(self):
        """
        This method subtracts the contribution of the sky background.
        """
        lst_target_keys = [elem for elem in list(self.target_dict.keys())]
        logger.info(f"Substracting sky background.")
        for elem in lst_target_keys:
            key, value = elem.split('+')
            lst_frames = self.target_dict[elem]
            cube = np.dstack(lst_frames)
            cube.sort(axis=2)
            im_avg = np.median(cube[:,:,:], axis=2)
            if key == 'target':
                logger.info(f"Creating sky background simulated for {value}.")
            elif key == 'fringe':
                logger.info(f"Creating sky background simulated for {value} without fringe.")
            else:
                logger.error("No defined option for key (target or fringe).")
            self.bkg = sep.Background(im_avg)
            no_sky = []
            for fr in lst_frames:
                no_sky.append(fr-self.bkg)
            self.target_dict['sky+' + value] = no_sky
            if key == 'std': #This is a special case for the STD stars.
                logger.info(f"List of photometric calibration frames without sky for {value} created.")
            elif key == 'target':
                logger.info(f"List of science frames without sky for {value} created.")
            elif key == 'fringe':
                logger.info(f"List of science frames without sky and without fringe for {value} created.")
            else:
                logger.error("No defined option for key (target or fringe).")


    def save_target(self, fringing=False, std=False, sky=False, not_sky=False):
        """This method saves the images generated during the cleaning process and 
        adds information to the header to assist in future processes.

        Args:
            fringing (bool, optional): Indicates whether the interference 
            pattern correction has been performed. Defaults to False.
            std (bool, optional): Indicates whether the image(s) to be
            saved are photometric calibration star images. Defaults to False.
            sky (bool, optional): Indicates whether the science image(s) to 
            be saved contain sky or not. Defaults to False.
        """
        if not std:
            logger.info("Saving science reduced frames.")
            lst_results = [elem for elem in list(self.DATA_DICT.keys()) if 'target' in elem]
        else:
            logger.info("Saving reduced photometric calibration frames.")
            lst_results = [elem for elem in list(self.DATA_DICT.keys()) if 'std' in elem]

        if 'std+OPEN' in lst_results:
            lst_results.remove('std+OPEN')

        for key in lst_results:
            fnames = self.DATA_DICT[key]

            if key == 'target+Sloan_z' and fringing:
                FRINGING = 'NO'
                sky_status = 'SKY'
                status='REDUCED'
                imagetype='SCIENCE'
                filt = key.split('+')[1]
                target = self.target_dict['fringe+Sloan_z']
                logger.info("Science frames without fringing.")
            elif std and not fringing:
                FRINGING = 'YES'
                sky_status = 'SKY'
                status='REDUCED'
                imagetype='STD'
                filt = key.split('+')[1]
                target= self.std_dict[key]
                logger.info("Photometric calibration frames.")
            elif (key == 'std+Sloan_z' and fringing and std):
                FRINGING = 'NO'
                sky_status = 'SKY'
                status='REDUCED'
                imagetype='STD'
                filt = key.split('+')[1]
                target = self.std_dict['fringe+Sloan_z']
                logger.info("Photometric calibration frames without fringing.")
            elif sky and not fringing:
                FRINGING = 'YES'
                sky_status = 'SKY'
                status='REDUCED'
                imagetype='SCIENCE'
                filt = key.split('+')[1]
                target= self.target_dict[key]
                logger.info("Science frames WITH sky")
            elif not_sky and not fringing:
                FRINGING = 'YES'
                sky_status = 'NOSKY'
                status='REDUCED'
                imagetype='SCIENCE'
                filt = key.split('+')[1]
                target= self.target_dict['sky+'+ filt]
                logger.info("Science frames WITHOUT sky")
            else:
                #logger.warning(f"{bcl.WARNING}No defined option for saving the frames (if fringing is true, ignore this message){bcl.ENDC}")
                #print("")
                #print(f'{bcl.HEADER}{lst_results}{bcl.ENDC}')
                #print(f'{bcl.HEADER}Key: {key}{bcl.ENDC}')
                #print(f'{bcl.HEADER}Fringing: {fringing}{bcl.ENDC}')
                #print(f'{bcl.HEADER}Sky: {sky}{bcl.ENDC}')
                #print(f'{bcl.HEADER}Not sky: {not_sky}{bcl.ENDC}')
                #print(f'{bcl.HEADER}STD: {std}{bcl.ENDC}')
                continue

            for i in range(len(fnames)):
                t = time.gmtime()
                time_string = time.strftime("%Y-%m-%dT%H:%M:%S", t)
                hd = fits.open(fnames[i])[0].header
                hd_wcs = wcs.WCS(fits.open(fnames[i])[0].header).to_header()
                hd['FRINGE'] = FRINGING
                hd['imgtype'] = imagetype
                hd['STATUS'] = status
                hd['SSKY'] = sky_status
                hd['BPMNAME'] = 'BPM_5sig' #Tipo de BPM aplicada.
                hd['rdate'] = time_string #Fecha de la reducción
                hd['filtro'] = filt #Simplifica la búsqueda por filtro.
                primary_hdu = fits.PrimaryHDU(target[i], header=hd+hd_wcs)
                hdul = fits.HDUList([primary_hdu])

                filename = os.path.basename(fnames[i])
                raw_name , __ = filename.split('.')

                logger.info(f"{bcl.OKGREEN}Storing the frame: ADP_{raw_name}_{imagetype}_{sky_status}_{filt} for {hd['FILTER2']}{bcl.ENDC}")
                hdul.writeto(str(self.PATH_RESULTS / f'ADP_{raw_name}_{imagetype}_{sky_status}_{filt}.fits'),
                            overwrite=True)