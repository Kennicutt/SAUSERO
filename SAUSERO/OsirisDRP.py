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

Copyright (C) 2024 Gran Telescopio Canarias <https://www.gtc.iac.es>
Fabricio Manuel Pérez Toledo <fabricio.perez@gtc.iac.es>
"""

__author__="Fabricio Manuel Pérez Toledo"
__version__ = "0.1.0"
__license__ = "GPL v3.0"

from reduction.reduction_osirisplus import *
from aligning.aligning_osirisplus import *
from astroMetry.astrometry_osirisplus import *
from photometry.photometry_osirisplus import *

from astropy import units as u

import argparse, time
import os, json, warnings

from Color_Codes import bcolors as bcl
from loguru import logger

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Parse configuration
parser = argparse.ArgumentParser(
                     prog = 'OsirisDRP',
                     description = 'This software reduces observations taken with OSIRIS\
                        in BBI mode. It can process any filter configuration and is suitable\
                        for observations affected by fringing (Sloan_z).')

parser.add_argument('-pr','--program', help='Select the GTC program.',
                     required=True, type=str)

parser.add_argument('-bl','--block', help='Select the block of the program.',
                     required=True, type=str)

args = parser.parse_args()

############## Predefined functions #############


def readJSON():
    return json.load(open("configuration.json"))

def Results(PATH, ZP, eZP, MASK, filt, ext_info = extinction_dict):
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
    ic = ccdp.ImageFileCollection(PATH, keywords='*', glob_include='ast*')
    for sky in ['SKY', 'NOSKY']:
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
        logger.info('Change units: ADUs to ADUs/second')
        
        frame.data = (frame.data / hd['EXPTIME'])
        frame.write(PATH / f"{hd['GTCPRGID']}_{hd['GTCOBID']}_{filt}_{(hd['DATE'].split('T')[0]).replace('-','')}_{sky}_BBI.fits",
                    overwrite=True)
        logger.info(f"Frame generated: {hd['GTCPRGID']}_{hd['GTCOBID']}_{filt}_{(hd['DATE'].split('T')[0]).replace('-','')}_{sky}_BBI.fits")
    

if __name__ == '__main__':

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

    PRG = args.program
    OB = args.block

    hora_local = time.localtime()
    conf = readJSON()

    logger.add(f"{conf['DIRECTORIES']['PATH_DATA']}{PRG}_{OB}/sausero_{time.strftime('%Y-%m-%d_%H:%M:%S', hora_local)}.log", format="{time} {level} {message} ({module}:{line})", level="INFO",
               filter=lambda record: 'astropy' not in record["name"])
    
    #Reduction Recipe. This recipe is responsible for cleaning the images by subtracting 
    #the masterbias and dividing by the normalized masterflat.
    #Subsequently, the cleaned images are saved.
    logger.info(f'{bcl.HEADER}---------- Starting the reduction ----------{bcl.ENDC}')
    o = Reduction(PRG, OB, main_path=conf['DIRECTORIES']['PATH_DATA'],
                path_mask="BPM/BPM_OSIRIS_PLUS.fits")
    o.get_imagetypes()
    o.load_BPM()
    o.sort_down_drawer()
    o.do_masterbias()
    o.do_masterflat()
    o.get_std(no_CRs=conf['REDUCTION']['no_CRs'], contrast_arg = conf['REDUCTION']['contrast'],
            cr_threshold_arg = conf['REDUCTION']['cr_threshold'],
            neighbor_threshold_arg = conf['REDUCTION']['neighbor_threshold'])
    o.get_target(no_CRs=conf['REDUCTION']['no_CRs'], contrast_arg = conf['REDUCTION']['contrast'],
            cr_threshold_arg = conf['REDUCTION']['cr_threshold'],
            neighbor_threshold_arg = conf['REDUCTION']['neighbor_threshold'])
    o.remove_fringing()
    o.sustract_sky()
    o.save_target(std=conf['REDUCTION']['save_std'])
    o.save_target(sky=conf['REDUCTION']['save_sky'])
    o.save_target(fringing=conf['REDUCTION']['save_fringing'])
    o.save_target()
    logger.info(f'{bcl.HEADER}¡¡¡¡¡¡¡ Finished the reduction successfully !!!!!!!{bcl.ENDC}')
    print(2*"\n")
    
    #Aligned Recipe. The cleaned science images are aligned based on the filter used in each case. 
    #Then, they are saved as aligned images.
    logger.info(f"{bcl.HEADER}---------- Starting the alignment ----------{bcl.ENDC}")
    al = OsirisAlign(PRG, OB, conf)
    for filt in list(set(al.ic.summary['filtro'])):
        for sky in ['SKY', 'NOSKY']:
            logger.info(f'{bcl.WARNING}++++++++++ Aligment for {filt} & {sky} ++++++++++{bcl.ENDC}')
            align = al.aligning(filt, sky=sky)
            lst = al.load_frames(filt, sky=sky)
            fr = CCDData.read(lst[0], unit='adu')
            header = fr.header
            header['exptime'] = al.total_exptime * (al.num + 1.)
            logger.info(f"Total exposure time estimated: {header['exptime']} sg")
            wcs = fr.wcs
            logger.info(f"Update the WCS information")
            save_fits(align, header, wcs, al.PATH_REDUCED / f'aligned_result_{filt}_{sky}.fits')
            print(f'{bcl.HEADER}¡¡¡¡¡¡¡ Alineado para {filt} & {sky} realizado exitosamente !!!!!!!{bcl.ENDC}')

    print(f'{bcl.HEADER}¡¡¡¡¡¡¡ El alineado para cada filtro finalizado !!!!!!!{bcl.ENDC}')
    print(2*"\n")

     #Astrometry Recipe. The aligned image for each filter undergoes an astrometric process to 
     #accurately determine the real positions of the celestial bodies present in the scene.
    logger.info(f"{bcl.HEADER}---------- Start the astrometrization ----------{bcl.ENDC}")
    del filt
    ic_ast = ccdp.ImageFileCollection(al.PATH_REDUCED, keywords='*', glob_include='ali*', glob_exclude='*NOSKY*')
    lst_filt = list(ic_ast.summary['filtro'])
    for filt in lst_filt:
        logger.info(f'{bcl.WARNING}++++++++++ Astrometrization for {filt} ++++++++++{bcl.ENDC}')
        best_wcs, new_frame = solving_astrometry(PRG, OB, filt, conf, sky='SKY', calib_std=False)
        for sky in ['SKY', 'NOSKY']:
            if sky == 'SKY':
                new_frame.write(al.PATH_REDUCED / f'ast_result_{filt}_{sky}.fits', overwrite=True)
            else:
                nosky = CCDData.read(al.PATH_REDUCED / f'aligned_result_{filt}_{sky}.fits', unit='adu')
                nosky.wcs = WCS(best_wcs)
                nosky.write(al.PATH_REDUCED / f'ast_result_{filt}_{sky}.fits', overwrite=True)
            print(f'{bcl.HEADER}¡¡¡¡¡¡¡ Astrometrization done successfully for {filt} !!!!!!!{bcl.ENDC}')
            time.sleep(10)
    
    logger.info(f'{bcl.HEADER}---------- Start astrometrization for STD star ----------{bcl.ENDC}')
    time.sleep(30)
    ic_std = ccdp.ImageFileCollection(al.PATH_REDUCED, keywords='*', glob_include='*std*', glob_exclude='*NOSKY*')
    best_wcs_std, new_frame_std = solving_astrometry(PRG, OB, filt, conf, sky='SKY', calib_std=True)
    for path_to_std in ic_std.files_filtered(include_path=True):
        logger.info(f'{bcl.WARNING}++++++++++ Astrometrization for: {path_to_std} ++++++++++{bcl.ENDC}')
        std_img = CCDData.read(path_to_std, unit='adu')
        std_img.wcs = WCS(best_wcs_std)
        std_img.write(path_to_std, overwrite=True)


    print(f'{bcl.HEADER}¡¡¡¡¡¡¡ Astrometrization finished !!!!!!!{bcl.ENDC}')
    print(2*"\n")


    #Photometry Recipe. The instrumental magnitude is estimated from the calibration star. This allows us 
    # to later estimate the apparent magnitude of the celestial bodies present in the science image. 
    # This process is done for each filter used. Additionally, the Results() function allows including 
    # this information in the header of the cleaned, aligned, and astrometrically processed science image.
    logger.info(f"{bcl.HEADER}---------- Starting the estimation of ZeroPoint ----------{bcl.ENDC}")
    del filt
    ic_pho = ccdp.ImageFileCollection(al.PATH_REDUCED, keywords='*', glob_include='red*')
    lst_filt = list(set(ic_pho.summary['filtro']))
    for filt in lst_filt:
        logger.info(f'{bcl.WARNING}++++++++++ Filter selected is {filt} ++++++++++{bcl.ENDC}')
        ZP, eZP = photometry(PRG, OB, ic_pho.files_filtered(imgtype='STD',filter2=filt)[0], conf)
        Results(al.PATH_REDUCED, ZP, eZP, o.MASK, filt)
    
    logger.info(f'{bcl.HEADER}¡¡¡¡¡¡¡ End of the reduction. The results are available in reduced directory. !!!!!!!{bcl.ENDC}')
