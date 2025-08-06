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

__author__="Fabricio M. Pérez-Toledo"
__version__ = "1.2.0"
__license__ = "GPL v3.0"

from SAUSERO.check_files import *
from SAUSERO.reduction_osirisplus import *
from SAUSERO.aligning_osirisplus import *
from SAUSERO.astrometry_osirisplus import *
from SAUSERO.photometry_osirisplus import *

from astropy import units as u

import argparse, time, os, shutil, re
import os, json, warnings
import pkg_resources
from pathlib import Path

from SAUSERO.Color_Codes import bcolors as bcl
from loguru import logger

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


############## Predefined functions #############

def create_config_file_home():
    """
    This function creates a copy of the configuration file in .config/sausero/ for easier accessibility.
    """
    config_path = pkg_resources.resource_filename(
    'SAUSERO', 'config/configuration.json')
    shutil.copy(config_path,Path(os.getcwd())/'configuration.json')
    print(f"{bcl.OKGREEN}Configuration file created successfully in the current directory.{bcl.ENDC}")
    sys.exit()

def readJSON():
    """
    Reads the file containing the configuration parameters.

    Returns:
        json: Collection of configuration parameters 
    """
    return json.load(open(Path(os.getcwd())/'configuration.json'))

    

def Results(PATH, ZP, eZP, MASK, filt, ext_info = extinction_dict, conf = None):
    """ This function adds relevant information to the science image that will be 
    delivered to the PI.This information includes the instrumental magnitude for 
    the used filter and its error, the units used, the extinction value applied and 
    its error, dividing the image by the exposure time, and applying the mask.

    Args:
        PATH (string): Directory containing the preliminary results.
        ZP (float): The estimated instrumental magnitude for the filter.
        eZP (float):  The error of the instrumental magnitude.
        MASK (bool): Mask of bad pixels.
        filt (string): Filter used for the image acquisition.
    """
    ic = ccdp.ImageFileCollection(PATH, keywords='*', glob_include='*ast*')
    for sky in ['SKY', 'NOSKY']:
        if conf['REDUCTION']['save_not_sky'] or sky == 'SKY':
            fname = ic.files_filtered(include_path=True, filtro=filt, ssky=sky)[0]
            frame = CCDData.read(fname, unit='adu')
            logger.info(f'Load science image: {os.path.basename(fname)}')
            hd = frame.header
            hd['ZP'] = (ZP, 'ZeroPoint estimation')
            hd['eZP'] = (eZP, 'Error ZeroPoint estimation')
            logger.info('ZeroPoint information added to header')
            extinction, e_extinction = ext_info[filt]
            hd['EXT']=(extinction, 'Filter extinction')
            hd['eEXT']=(e_extinction, 'Error filter extinction')
            logger.info('Extinction information added to header')
            frame.header = hd
            frame.unit = u.adu/u.second
            hd['PHOTOMETRY'] = (True, 'Photometry applied')
            logger.info('Change units: ADUs to ADUs/second')

            frame.data = (frame.data / hd['EXPTIME'])
            frame.write(PATH / f"{hd['GTCPRGID']}_{hd['GTCOBID']}_{filt}_pho_{sky}.fits",
                        overwrite=True)
            logger.info(f"{bcl.OKGREEN}Frame generated: {hd['GTCPRGID']}_{hd['GTCOBID']}_{filt}_pho_{sky}.fits{bcl.ENDC}")
        else:
            logger.warning(f'{bcl.WARNING}The photometry is not going to be executed for NOSKY{bcl.ENDC}')
    

def run():
    """
    This function 
    """
    # Parse configuration
    parser = argparse.ArgumentParser(
                         prog = 'OsirisDRP',
                         description = 'This software reduces observations taken with OSIRIS\
                            in BBI mode. It can process any filter configuration and is suitable\
                            for observations affected by fringing (Sloan_z).')

    parser.add_argument('-e', '--execute', help='Execute the configuration file in the current directory.',
                        action='store_true')
    
    parser.add_argument('-c', '--create_config', help='Create a configuration file in the current directory.',
                        action='store_true')

    args = parser.parse_args()


    print(f"{bcl.OKBLUE}***********************************************************************{bcl.ENDC}")
    print(f"{bcl.OKBLUE}************************* WELCOME TO SAUSERO **************************{bcl.ENDC}")
    print(f"{bcl.OKBLUE}***********************************************************************{bcl.ENDC}")
    print("\n")
    print(f"{bcl.BOLD}---------------------- LICENSE ----------------------{bcl.ENDC}")
    print("\n")
    print(f"This program is free software: you can redistribute it and/or modify\n\
it under the terms of the GNU General Public License as published by\n\
the Free Software Foundation, either version 3 of the License, or\n\
(at your option) any later version.\n\n\
This program is distributed in the hope that it will be useful,\n\
but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n\
GNU General Public License for more details.\n\n\
You should have received a copy of the GNU General Public License\n\
along with this program. If not, see <https://www.gnu.org/licenses/>.")
    print("\n")
    print(f"{bcl.BOLD}************************ IMPORTANT INFORMATION ************************{bcl.ENDC}")
    print("\n")
    print(f"This software is designed to reduce Broad Band Imaging observations obtained with OSIRIS+.\n\
For proper use, you need to modify the configuration file, which can be found\n\
in the directory where this software is installed. Additionally, you need to create\n\
an account on Astrometry.net. Once you have the code that allows you to use the API,\n\
you need to fill in the correct variable.")
    print(f"\n")

    # Check if the configuration file exists (2025-08-04)
    if args.create_config:
        print(f"{bcl.OKGREEN}Creating the configuration file in the current directory.{bcl.ENDC}")
        print(f"{bcl.WARNING}You can edit it before you execute the reduction.{bcl.ENDC}")
        create_config_file_home()
        sys.exit()

    ########## Checking files (2025-01-22) ##########
    conf = check_files()

    PRG = conf['PRG']
    OB = conf['OB']

    hora_local = time.localtime()
    print(conf)
    logger.add(Path(conf['DIRECTORIES']['PATH'])/f"sausero_{time.strftime('%Y-%m-%d_%H:%M:%S', hora_local)}.log", format="{time} {level} {message} ({module}:{line})", level="INFO",
            filter=lambda record: 'astropy' not in record["name"])
        
    logger.info(f'{bcl.OKGREEN}Log file created{bcl.ENDC}')
    logger.info(f'{bcl.OKGREEN}Configuration has been updated successfully{bcl.ENDC}')
    logger.info(f'{bcl.OKGREEN}Read the configuration file successfully{bcl.ENDC}')
    
    #Reduction Recipe. This recipe is responsible for cleaning the images by subtracting 
    #the masterbias and dividing by the normalized masterflat.
    #Subsequently, the cleaned images are saved.
    logger.info(f'{bcl.OKBLUE}---------- Starting the reduction for {PRG}-{OB} ----------{bcl.ENDC}')
    
    bpm_path = pkg_resources.resource_filename('SAUSERO', 'BPM/BPM_OSIRIS_PLUS.fits')
    o = Reduction(main_path=conf['DIRECTORIES']['PATH_DATA'],
                path_mask=bpm_path)
    o.get_imagetypes()
    o.load_BPM()
    o.sort_down_drawer()

    if conf['REDUCTION']['use_BIAS']:
        o.do_masterbias()

    else:
        logger.warning(f'{bcl.WARNING}The masterbias is not going to be created{bcl.ENDC}')

    if conf['REDUCTION']['use_FLAT']:
        o.do_masterflat()

    else:
        logger.warning(f'{bcl.WARNING}The masterflat is not going to be created{bcl.ENDC}')
    
    if conf['REDUCTION']['use_STD']:
        o.get_std(no_CRs=conf['REDUCTION']['no_CRs'], contrast_arg = conf['REDUCTION']['contrast'],
            cr_threshold_arg = conf['REDUCTION']['cr_threshold'],
            neighbor_threshold_arg = conf['REDUCTION']['neighbor_threshold'], apply_flat=conf['REDUCTION']['use_FLAT'])

    else:
        logger.warning(f'{bcl.WARNING}The STD star is not going to be reduced{bcl.ENDC}')
    
    o.get_target(no_CRs=conf['REDUCTION']['no_CRs'], contrast_arg = conf['REDUCTION']['contrast'],
            cr_threshold_arg = conf['REDUCTION']['cr_threshold'],
            neighbor_threshold_arg = conf['REDUCTION']['neighbor_threshold'], apply_flat=conf['REDUCTION']['use_FLAT'])
    
    
    if conf['REDUCTION']['save_fringing']:
        o.remove_fringing()
        logger.info(f'{bcl.OKGREEN}Fringing correction applied successfully{bcl.ENDC}')
    else:
        logger.warning(f'{bcl.WARNING}The fringing correction is not going to be executed{bcl.ENDC}')

    if conf['REDUCTION']['save_not_sky']:    
        o.sustract_sky()
        logger.info(f'{bcl.OKGREEN}The sky subtraction has been applied successfully{bcl.ENDC}')
    else:
        logger.warning(f'{bcl.WARNING}The sky substraction is not going to be executed{bcl.ENDC}')
    
    o.save_target(std=conf['REDUCTION']['save_std'])
    o.save_target(std=conf['REDUCTION']['save_std'], fringing=conf['REDUCTION']['save_fringing'])
    o.save_target(sky=conf['REDUCTION']['save_sky'])
    o.save_target(fringing=conf['REDUCTION']['save_fringing'])    
    o.save_target(not_sky=conf['REDUCTION']['save_not_sky'])
    logger.info(f'{bcl.OKBLUE}-------------- End of the reduction successfully --------------{bcl.ENDC}')
    print(2*"\n")

#-----------------------------------------------------------------------------------------------------------------------------------    
    #Aligned Recipe. The cleaned science images are aligned based on the filter used in each case. 
    #Then, they are saved as aligned images.
    if conf['ALIGNING']['use_aligning']:
        logger.info(f"{bcl.OKBLUE}---------- Starting the alignment ----------{bcl.ENDC}")
        al = OsirisAlign(conf)
        for filt in list(set(al.ic.summary['filtro'])):
            for sky in ['SKY', 'NOSKY']:
                if conf['REDUCTION']['save_not_sky'] or sky == 'SKY':
                    logger.info(f'{bcl.OKCYAN}++++++++++ Aligment for {filt} & {sky} ++++++++++{bcl.ENDC}')
                    align = al.aligning(filt, sky=sky)
                    lst = al.load_frames(filt, sky=sky)
                    fr = CCDData.read(lst[0], unit='adu')
                    header = fr.header
                    header['STACKED'] = (True, 'Stacked image')
                    header['exptime'] = al.total_exptime * al.num #(al.num + 1.)
                    logger.info(f"Estimated total exposure time: {header['exptime']} sec")
                    wcs = fr.wcs
                    logger.info(f"Updating the WCS information")
                    save_fits(align, header, wcs, str(al.PATH_REDUCED / f'{PRG}_{OB}_{filt}_stacked_{sky}.fits'))
                    
                else:
                    logger.warning(f'{bcl.WARNING}Alignment is not going to be executed for NOSKY{bcl.ENDC}')

        logger.info(f'{bcl.OKBLUE}------------------- End of the alignment -------------------{bcl.ENDC}')
        print(2*"\n")
    else:
        logger.warning(f'{bcl.WARNING}The alignment is not going to be executed{bcl.ENDC}')

#-----------------------------------------------------------------------------------------------------------------------------------
     #Astrometry Recipe. The aligned image for each filter undergoes an astrometric process to 
     #accurately determine the real positions of the celestial bodies present in the scene.
    if conf['ASTROMETRY']['use_astrometry']:
        logger.info(f"{bcl.OKBLUE}---------- Start the astrometrization ----------{bcl.ENDC}")
        del filt
        ic_ast = ccdp.ImageFileCollection(al.PATH_REDUCED, keywords='*', glob_include='*stacked*', glob_exclude='*NOSKY*')
        lst_filt = list(ic_ast.summary['filtro'])
        for filt in lst_filt:
            logger.info(f'{bcl.OKCYAN}++++++++++ Astrometrization for the stacked image with {filt} filter ++++++++++{bcl.ENDC}')
            try:
                best_wcs, new_frame = solving_astrometry(PRG, OB, filt, conf, sky='SKY', calib_std=False)
                logger.info(f'{bcl.OKGREEN}New WCS for the stacked image with {filt} filter.{bcl.ENDC}')
            except:
                new_frame = CCDData.read(al.PATH_REDUCED / f'{PRG}_{OB}_{filt}_stacked_SKY.fits', unit='adu')
                best_wcs = None
                logger.error(f'{bcl.ERROR}Failed astrometrization for the stacked image with {filt} filter{bcl.ENDC}')
            for sky in ['SKY', 'NOSKY']:
                if sky == 'SKY':
                    if best_wcs is not None:
                        new_frame.header['ASTROMETRY'] = (True, 'Astrometrized image')
                        new_frame.write(al.PATH_REDUCED / f'{PRG}_{OB}_{filt}_ast_SKY.fits', overwrite=True)
                        logger.info(f'{bcl.OKGREEN}Successful astrometrization for the stacked image with {filt} and SKY{bcl.ENDC}')
                    else:
                        new_frame.header['ASTROMETRY'] = (False, 'Astrometrized image')
                        new_frame.write(al.PATH_REDUCED / f'{PRG}_{OB}_{filt}_ast_SKY.fits', overwrite=True)
                        logger.warning(f'{bcl.WARNING}Failed astrometrization for the stacked image with {filt} and SKY. Conservation of original WCS{bcl.ENDC}')
                else:
                    if conf['REDUCTION']['save_not_sky']:
                        nosky = CCDData.read(al.PATH_REDUCED / f'{PRG}_{OB}_{filt}_stacked_NOSKY.fits', unit='adu')
                        if best_wcs is not None:
                            nosky.wcs = WCS(best_wcs)
                            nosky.header['ASTROMETRY'] = (True, 'Astrometrized image')
                            nosky.write(al.PATH_REDUCED / f'{PRG}_{OB}_{filt}_ast_NOSKY.fits', overwrite=True)
                            logger.info(f'{bcl.OKGREEN}Successful astrometrization for the stacked image with {filt} and NOSKY{bcl.ENDC}')
                        else:
                            nosky.wcs = nosky.wcs
                            nosky.header['ASTROMETRY'] = (False, 'Astrometrized image')
                            nosky.write(al.PATH_REDUCED / f'{PRG}_{OB}_{filt}_ast_NOSKY.fits', overwrite=True)
                            logger.warning(f'{bcl.WARNING}Failed astrometrization for the stacked image with {filt} and NOSKY. Conservation of original WCS{bcl.ENDC}')
                    else:
                        logger.warning(f'{bcl.WARNING}The astrometry is not going to be executed for NOSKY{bcl.ENDC}')

                logger.info(f'{bcl.OKCYAN}END of the astrometrization for {filt} and {sky}{bcl.ENDC}')
                time.sleep(10)
    else:
        logger.warning(f'{bcl.WARNING}The astrometry is not going to be executed{bcl.ENDC}')

    if conf['REDUCTION']['use_STD'] and conf['ASTROMETRY']['use_astrometry']:
        logger.info(f'{bcl.OKCYAN}---------- Start astrometrization for STD star ----------{bcl.ENDC}')
        time.sleep(30)
        ic_std = ccdp.ImageFileCollection(al.PATH_REDUCED, keywords='*', glob_include='*STD*', glob_exclude='*NOSKY*')
        lst_object = list(set(ic_std.summary['object']))
        try:
            best_wcs_std, _ = solving_astrometry(PRG, OB, filt, conf, sky='SKY', calib_std=True)
            logger.info(f'{bcl.OKGREEN}New WCS for the STD star with {filt} filter.{bcl.ENDC}')
        except:
            best_wcs_std = None
            logger.error(f'{bcl.ERROR}Failed astrometrization for the STD star with {filt} filter{bcl.ENDC}')

        for path_to_std in ic_std.files_filtered(include_path=True):
            std_img = CCDData.read(path_to_std, unit='adu')
            if best_wcs_std is not None:
                std_img.wcs = WCS(best_wcs_std)
                std_img.header['ASTROMETRY'] = (True, 'Astrometrized image')
                std_img.write(path_to_std, overwrite=True)
                logger.info(f'{bcl.OKGREEN}Successful astrometrization done for the {std_img.header["OBJECT"]} with {std_img.header["FILTER2"]}{bcl.ENDC}')
            else:
                std_img.wcs = std_img.wcs
                std_img.header['ASTROMETRY'] = (False, 'Astrometrized image')
                std_img.write(path_to_std, overwrite=True)
                logger.warning(f'{bcl.WARNING}Failed astrometrization for the {std_img.header["OBJECT"]} with {std_img.header["FILTER2"]}. Conserve the original WCS{bcl.ENDC}')


        logger.info(f'{bcl.OKCYAN}END of astrometrization for STD{bcl.ENDC}')
        print(2*"\n")
    else:
        logger.warning(f'{bcl.WARNING}The astrometry for STDs are not going to be executed{bcl.ENDC}')

    logger.info(f'{bcl.OKBLUE}------------------- End of the astrometrization -------------------{bcl.ENDC}')
    print(2*"\n")

#-----------------------------------------------------------------------------------------------------------------------------------
    #Photometry Recipe. The instrumental magnitude is estimated from the calibration star. This allows us 
    # to later estimate the apparent magnitude of the celestial bodies present in the science image. 
    # This process is done for each filter used. Additionally, the Results() function allows including 
    # this information in the header of the cleaned, aligned, and astrometrically processed science image.
    if conf['PHOTOMETRY']['use_photometry']:
        logger.info(f"{bcl.OKBLUE}---------- Starting the estimation of ZeroPoint ----------{bcl.ENDC}")
        del filt
        ic_pho = ccdp.ImageFileCollection(al.PATH_REDUCED, keywords='*', glob_include='*ADP*')
        lst_filt = list(set(ic_pho.summary['filtro']))
        for filt in lst_filt:
            if filt == "OPEN":
                logger.warning(f'{bcl.WARNING}The OPEN filter is not going to be used for photometry{bcl.ENDC}')
                continue
            logger.info(f'{bcl.OKCYAN}++++++++++ Filter selected is {filt} ++++++++++{bcl.ENDC}')
            try:
                ZP, eZP = photometry(PRG, OB, ic_pho.files_filtered(imgtype='STD',filter2=filt)[0], conf)
                Results(al.PATH_REDUCED, ZP, eZP, o.MASK, filt, conf=conf)
                logger.info(f'{bcl.OKGREEN}Photometry done for {filt} filter{bcl.ENDC}')
            except:
                ZP, eZP = None, None
                Results(al.PATH_REDUCED, ZP, eZP, o.MASK, filt, conf=conf)
                logger.warning(f'{bcl.WARNING}Failed photometry for {filt} filter{bcl.ENDC}')
    
    else:
        logger.warning(f'{bcl.WARNING}The photometry is not going to be executed{bcl.ENDC}')

    logger.info(f'{bcl.OKBLUE}------------------- End of the photometry -------------------{bcl.ENDC}')
    print(2*"\n")
    # Final message
    logger.info(f'{bcl.OKBLUE}End of the reduction. The results are available in {conf["DIRECTORIES"]["PATH_OUTPUT"]}{bcl.ENDC}')
    
    for archivo in glob.glob(str(al.PATH_REDUCED/"ADP*.fits")):
        nuevo_nombre = re.sub(r"_(Sloan)_[a-zA-Z]+", "", archivo)
        os.rename(archivo, nuevo_nombre)

    try:
        for archivo in glob.glob(str(al.PATH_REDUCED/"ADP*OPEN*.fits")):
            nuevo_nombre = re.sub(r"_(OPEN)_+", "", archivo)
            os.rename(archivo, nuevo_nombre)
    except:
        logger.info(f'No OPEN filter')

    print(2*"\n")
    print(f"{bcl.OKBLUE}************************* THANK YOU FOR USING SAUSERO *************************{bcl.ENDC}")
    print(2*"\n")
    print(f"{bcl.OKBLUE}***********************************************************************{bcl.ENDC}")
    print(f"{bcl.OKBLUE}************************* END OF SAUSERO **************************{bcl.ENDC}")
    print(f"{bcl.OKBLUE}***********************************************************************{bcl.ENDC}")

#--------------------------------------------------------------------
if __name__ == '__main__':
    run()
