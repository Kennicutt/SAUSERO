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
"""

import astroalign as aa
from astropy.nddata import CCDData
from astropy.table import Table
import ccdproc as ccdp
from pathlib import Path
import time
import matplotlib.pyplot as plt
from astropy.visualization import LogStretch,imshow_norm, ZScaleInterval

from Color_Codes import bcolors as bcl
from loguru import logger


class OsirisAlign:
    """Esta clase permite alinear las imágenes de ciencia (o de calibración fotométrica si procede).
    De esta forma podemos sumarlas y así obtener más flujo. Proceso vital para observar cuerpos débiles.
    """
    
    def __init__(self, program, block, conf):
        """Inicializamos la clase definiendo parámetros importantes.

        Args:
            program (str): Programa científico a alinear. Cómo incluir cambios
            block (str): Número del bloque del programa científico.
        """
        self.conf = conf
        self.PATH = self.conf["DIRECTORIES"]["PATH_DATA"]
        self.program = program
        self.block = block
        self.reduced_data_directory = "reduced"
        self.PATH_REDUCED = Path(self.PATH + self.program + "_" + self.block + "/" + "reduced/")
        self.PATH_TARGET = Path(self.PATH + self.program + "_" + self.block + "/" + "object/")
        self.ic = ccdp.ImageFileCollection(self.PATH_REDUCED, keywords='*', glob_include='red*')


    def load_frames(self, filt, sky):
        """Este método pretende obtener una lista de imágenes de ciencia para un filtro
        determinado.

        Args:
            filt (str): Filtro de interés.

        Returns:
            list: Lista de imágenes científicas para un determinado filtro y su ruta.
        """
        self.tab = Table(self.ic.summary)
        sub_tab = self.tab[(self.tab['filter2']==filt) & (self.tab['ssky']==sky)]
        logger.info(f"Looking for frames with {filt} and {sky}")
        logger.info(f"{sub_tab}")
        self.total_exptime = sub_tab['exptime'].value.data[0]
        logger.info(f"Exposure time per frame: {self.total_exptime} sg")
        return self.ic.files_filtered(imgtype="SCIENCE",
                                      filtro=filt,
                                      ssky = sky,
                                      include_path=True)



    def get_each_data(self, filt, sky):
        """Lee el contenido de una lista que contiene
        la ruta a una serie de frames y abre estos frames
        añadiéndolos a una lista nueva.

        Args:
            filt (str): Filtro de interés.

        Returns:
            list: Lista de imágenes de ciencia para un filtro (matrices).
        """
        ccd = []
        for frame_path in self.load_frames(filt, sky=sky):
            ccd.append(CCDData.read(frame_path, unit='adu', hdu=0).data)

        return ccd



    def aligning(self, filt, sky='SKY'): #default: 30
        """Este método tiene por objetivo alinear las imágenes de ciencia que han sido adquiridas
        con un mismo filtro.

        Args:
            filt (str): Filtro de interés.

        Returns:
            float: Una imagen que está compuesta por varias imágenes alineadas incrementando el flujo
            de la misma.
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
    """Permite visualizar el contenido

    Args:
        cube (float): Imagen alineada.
        a (int, optional): Factor a aplicar. Defaults to 1.
    """
    fig = plt.figure(figsize=(15,20))
    ax = fig.add_subplot(1, 1, 1)
    im, norm = imshow_norm(cube, ax, origin='lower',
                           interval=ZScaleInterval(),
                           stretch=LogStretch(a=a))
    fig.colorbar(im)
    plt.show()
    
def save_fits(image, header, wcs, fname):
    """Permite salvar la imagen generada.

    Args:
        image (float): Matriz que representa la imagen.
        header (str): Header para la imagen.
        wcs (str): WCS para la imagen.
        fname (str): Nombre a dar a la imagen.
    """
    ccd = CCDData(data=image, header=header, wcs=wcs, unit='adu')
    ccd.write(fname, overwrite=True)
    logger.info(f"New image has been created: {fname}")