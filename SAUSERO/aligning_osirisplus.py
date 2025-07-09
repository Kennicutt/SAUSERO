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

import astroalign as aa
from astropy.nddata import CCDData
from astropy.table import Table
import ccdproc as ccdp
from pathlib import Path
import time
import matplotlib.pyplot as plt
from astropy.visualization import LogStretch,imshow_norm, ZScaleInterval

from SAUSERO.Color_Codes import bcolors as bcl
from loguru import logger


class OsirisAlign:
    """This class allows the alignment of science frames (or photometry calibration frames, if applicable). 
    Afterward, we can stack them to enhance the measured flux. This step is essential for observing 
    faint sources.
    """
    
    def __init__(self, program, block, conf):
        """We initialize the class by defining important parameters.

        Args:
            program (str): Science program code
            block (str): Observational block number assigned to a science program
        """
        self.conf = conf
        self.PATH = self.conf["DIRECTORIES"]["PATH_DATA"]
        self.program = program
        self.block = block
        self.reduced_data_directory = "reduced"
        self.PATH_REDUCED = Path(self.PATH + self.program + "_" + self.block + "/" + "reduced/")
        #self.PATH_TARGET = Path(self.PATH + self.program + "_" + self.block + "/" + "object/")
        self.ic = ccdp.ImageFileCollection(self.PATH_REDUCED, keywords='*', glob_include='ADP*')


    def load_frames(self, filt, sky):
        """This method retrieves a list of science frames for each filter.

        Args:
            filt (str): Filter name

        Returns:
            list: A list of science frames for a given filter and its path
        """
        self.tab = Table(self.ic.summary)
        sub_tab = self.tab[(self.tab['filter2']==filt) & (self.tab['ssky']==sky)]
        logger.info(f"Looking for frames that have {filt} and {sky}")
        logger.info(f"{sub_tab}")
        self.total_exptime = sub_tab['exptime'].value.data[0]
        logger.info(f"Exposure time per frame: {self.total_exptime} sec")
        return self.ic.files_filtered(imgtype="SCIENCE",
                                      filtro=filt,
                                      ssky = sky,
                                      include_path=True)



    def get_each_data(self, filt, sky):
        """This methods reads a list containing the paths to several frames 
        and opens them to add them to a new list.

        Args:
            filt (str): Filter name

        Returns:
            list: List of science frames for each filter (matrices)
        """
        ccd = []
        for frame_path in self.load_frames(filt, sky=sky):
            ccd.append(CCDData.read(frame_path, unit='adu', hdu=0).data)

        return ccd



    def aligning(self, filt, sky='SKY'): #default: 30
        """This method aligns the science frames taken with the same filter.

        Args:
            filt (str): Filter name

        Returns:
            float: Stacked image obtained by combining multiple science frames.
        """
        logger.info(f"Creating cube with frames for {filt}")
        cube = self.get_each_data(filt, sky=sky)
        REF = cube[0]

        self.num=0
        logger.info(f"Number of frames in cube is: {len(cube)}")
        for IMG in cube[1:]:
            try:
                t, __ = aa.find_transform(IMG, REF, 
                                          max_control_points=self.conf["ALIGNING"]["max_control_points"]) #default: 30
                time.sleep(1)
                REF = REF.byteswap().newbyteorder()
                IMG = IMG.byteswap().newbyteorder()
                align, footprint = aa.apply_transform(t, IMG, REF)
                REF = REF + align
                self.num += 1
                logger.info(f"Image NO: {self.num}/{len(cube)}")
            except aa.MaxIterError:
                #print(aa.MaxIterError)
                logger.error(f"{bcl.FAIL}ERROR{bcl.ENDC}: {aa.MaxIterError}")
                pass
            except ValueError:
                #print(ValueError)
                logger.error(f"{bcl.FAIL}ERROR{bcl.ENDC}: {ValueError}")
                pass
            except TypeError:
                #print(TypeError)
                logger.error(f"{bcl.FAIL}ERROR{bcl.ENDC}: {TypeError}")
                pass

        return REF


def show_picture(cube, a=1):
    """This method allows us to show the image.

    Args:
        cube (float): Stacked image.
        a (int, optional): Factor applied. Defaults: 1.
    """
    fig = plt.figure(figsize=(15,20))
    ax = fig.add_subplot(1, 1, 1)
    im, norm = imshow_norm(cube, ax, origin='lower',
                           interval=ZScaleInterval(),
                           stretch=LogStretch(a=a))
    fig.colorbar(im)
    plt.show()
    
def save_fits(image, header, wcs, fname, allow_nosky=True):
    """This method saves the stacked image.

    Args:
        image (float): Data to save.
        header (str): Header for the stacked image.
        wcs (str): WCS information for the stacked image.
        fname (str): Name for the stacked image.
    """
    header['STACKED'] = 'YES'
    ccd = CCDData(data=image, header=header, wcs=wcs, unit='adu')
    ccd.write(fname, overwrite=True)
    logger.info(f"New image has been created: {fname}")