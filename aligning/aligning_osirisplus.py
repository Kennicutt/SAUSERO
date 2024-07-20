import astroalign as aa
from astropy.nddata import CCDData
from astropy.table import Table
import ccdproc as ccdp
from pathlib import Path
import time
import matplotlib.pyplot as plt
from astropy.visualization import LogStretch,imshow_norm, ZScaleInterval


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
        self.total_exptime = sub_tab['exptime'].value.data[0]
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
        cube = self.get_each_data(filt, sky=sky)
        REF = cube[0]

        self.num=0
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
                print(f"Image NO: {self.num}")
            except aa.MaxIterError:
                print(aa.MaxIterError)
                pass
            except ValueError:
                print(ValueError)
                pass
            except TypeError:
                print(TypeError)
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