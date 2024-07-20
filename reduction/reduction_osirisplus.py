import os, sys, time, json
from pathlib import Path

from astropy import units as u
from astropy import wcs
from astropy.nddata import CCDData
from astropy.io import fits
import ccdproc as ccdp
from matplotlib import pyplot as plt
import numpy as np
import yaml as py
import lacosmic
import sep

class Reduction:
    """El objetivo de esta clase es llevar a cabo el proceso de limpieza de las imágenes de ciencia y
    de calibración fotométrica.
    """

    def __init__(self, gtcprgid, gtcobid, main_path, path_mask = None):
        """Inicializamos la Clase mediante la definición de una serie de parámetros a partir de los que
        se definirán otros de relevante importancia para el correcto proceso de limpieza de las imágenes.

        Args:
            gtcprgid (str): Programa de observación que se desea reducir.
            gtcobid (str): Bloque del programa definido que se desea reducir.
            path_mask (str, optional): Directorio que contiene la BPM. 
            Defaults to "/home/phase3/Fabricio-drm/SAUSERO/BPM/BPM_OSIRIS_CMS_sig5.fits".
        """
        print(f'Initializing the reduction for {gtcprgid}_{gtcobid}.')
        self.gtcprgid = gtcprgid
        self.gtcobid = gtcobid
        ROOT = main_path
        #Se define a la ruta que contiene las imágenes crudas:
        self.PATH = Path(ROOT + str(gtcprgid) + '_' + str(gtcobid) + '/' + 'raw/')
        #Se define el directorio que contendrá las imágenes intermedias y procesadas:
        self.PATH_RESULTS = Path(ROOT + str(gtcprgid) + '_' + str(gtcobid) + '/' + 'reduced/')
        self.PATH_RESULTS.mkdir(parents=True, exist_ok=True)
        self.path_mask = Path(path_mask)
        #Se recopila la información sobre los frames ubicados en dicho directorio.
        self.ic = ccdp.ImageFileCollection(self.PATH)
        self.DATA_DICT={}
        self.key_dict={'flat':'OsirisSkyFlat',
                       'target': 'OsirisBroadBandImage',
                       'std': 'OsirisBroadBandImage',
                       'bias': 'OsirisBias'
        }
        self.master_dict = {} #Listo para usar
        self.std_dict = {} #Listo para usar
        self.target_dict = {} #Listo para ciencia




    @staticmethod
    def configure_mask(mask):
        """Método estático que se emplea para adecuar la BPM
        para su uso en las imágenes.

        Args:
            mask (int): BPM

        Returns:
            bool: BPM en términos booleanos para aplicar sobre
            las imágenes.
        """
        mask = mask.data.astype(bool)
        matrix = np.ones(mask.shape)
        matrix[mask == False] = np.nan
        return matrix[230:2026,28:2060] # TRIM SECTION




    @staticmethod
    def get_each_data(data_dict, value):
        """Lee el contenido de un diccionario que contiene
        la ruta a una serie de frames y abre estos frames
        añadiéndolos a una lista.

        Args:
            data_dict (dict): Diccionario que contiene las
            rutas a uno o varios frames contenidos en una lista.
            value (str): Key para acceder a lista de las rutas. 

        Returns:
            list: Es una lista que contiene las imágenes (matrices).
        """
        ccd = []
        for frame_path in data_dict[value]:
            hdul = fits.open(frame_path)
            ccd.append(hdul[0].data[230:2026,28:2060]) #TRIM SECTION
        return ccd



    @staticmethod
    def combining(lst_frames):
        """Este método estático combinan las imágenes contenidas
        en una lista para obtener una imagen promediada. La
        estrategia consiste en crear un cubo de datos y promediarlos.

        Args:
            lst_frames (list): Lista de imágenes a promediar.

        Returns:
            float: Genera una matriz promediada a partir del cubo de
            imágenes.
        """
        cube = np.dstack(lst_frames)
        cube.sort(axis=2)
        return np.nanmedian(cube[:,:,1:-1], axis=2)




    def create_cubes(self, key):
        """Genera un cubo de datos a partir de una lista de
        imágenes (matrices).

        Args:
            key (str): Key del diccionario de imágenes científicas.

        Returns:
            cube(float): Cubo de imágenes en profundidad.
        """
        return np.dstack(self.target_dict[key],
                        axis=2)






    def sustractMasterBias(self, value, master, data_dict):
        """Sustracción del masterbias a una imagen.

        Args:
            value (str): Key que permite acceder a las imágenes
            del diccionario.
            master (float): Matriz que representa al MasterBias.
            data_dict (dict): Diccionario que alberga las imágenes. 

        Returns:
            list: Lista de frames con el MasterBias aplicado.
        """
        ccd = self.get_each_data(data_dict, value)

        frames = [fr - master for fr in ccd]
        return frames






    def clean_target(self, *args):
        """Aplica a la imágenes de ciencia la sustracción del MasterBias y
        la división por el MasterFlat normalizado.

        Returns:
            list: Lista de imágenes cienfíficas limpias.
        """
        value, masterbias, masterflat = args
        ccd = self.get_each_data(self.DATA_DICT, value)
        frames = [(fr - masterbias)/masterflat for fr in ccd]
        return frames


    def get_imagetypes(self):
        """ Este método tiene por objetivo establecer todos los tipos de imágenes que se encuentran
        presentes en el directorio de las imágenes originales y sus correspondientes filtros empleados.
        Para ello crea un diccionario vacío que cada key está constituida por el tipo de imagen y el
        filtro empleado.
        """
        print('Create all components ')
        self.filt_wheels = []
        matches = (self.ic.summary['obsmode'] != 'OsirisBias')
        matches1 = (self.ic.summary['obsmode'] != 'OsirisBias') & (self.ic.summary['filter1'] != 'OPEN')
        matches2 = (self.ic.summary['obsmode'] != 'OsirisBias') & (self.ic.summary['filter2'] != 'OPEN')
        matches3 = (self.ic.summary['obsmode'] != 'OsirisBias') & (self.ic.summary['filter3'] != 'OPEN')

        if len(list(set(self.ic.summary['filter1'][matches1]))) >= 1:
            self.filt_wheels.append('filter1')
        elif len(list(set(self.ic.summary['filter2'][matches2]))) >= 1:
            self.filt_wheels.append('filter2')
        elif len(list(set(self.ic.summary['filter3'][matches3]))) >= 1:
            self.filt_wheels.append('filter3')
        else:
            raise(ValueError, "All filters are open!!!")
            sys.exit()

        for filt in self.filt_wheels:
            for value in list(set(self.ic.summary[filt][matches])):
                for type in ['flat', 'std', 'target']:
                    self.DATA_DICT[type + '+' + value] = []

        self.DATA_DICT["bias"] = []


    def load_BPM(self):
        """
        Este método abre el FITS que contiene la BPM.
        """
        bpm = CCDData.read(self.path_mask, unit=u.dimensionless_unscaled,
                            hdu=1)
        
        #bpm = np.ones((2056,2073))

        self.MASK = self.configure_mask(bpm)




    def load_results(self):
        """Genera una tabla de contenido específico en el directorio
        de resultados. En concreto, de aquellas imágenes científicas
        que ya han sido limpiadas.
        """
        self.ic_r = ccdp.ImageFileCollection(self.PATH_RESULTS, keywords='*',
                                             glob_include='red*')




    def sort_down_drawer(self):
        """Este método se encarga de almacenar las imágenes en función de su tipo y filtro en el
        diccionario creado anteriormente.
        """
        for filt in self.filt_wheels:
            for elem in list(self.DATA_DICT.keys()):
                if elem != 'bias':
                    key, value = elem.split('+')
                    types_targets = set(self.ic.summary['object'][self.ic.summary['obsmode'] == 'OsirisBroadBandImage'])
                    if key == 'flat':
                        tmp_dict = {"obsmode":self.key_dict[key], filt: value}
                    elif key == 'std':
                        target_type = [data for data in types_targets if 'STD' in data][0]
                        tmp_dict = {"obsmode":self.key_dict[key], filt: value, 'object':target_type}
                    elif key == 'target':
                        target_type = [data for data in types_targets if not 'STD' in data][0]
                        tmp_dict = {"obsmode":self.key_dict[key], filt: value, 'object':target_type}
                    self.DATA_DICT[elem] = self.ic.files_filtered(**tmp_dict, include_path=True)
                else:
                    self.DATA_DICT[elem] = self.ic.files_filtered(obsmode='OsirisBias', include_path=True)



    def do_masterbias(self):
        """Este método permite crear el MasterBias.
        """
        bias = self.get_each_data(self.DATA_DICT, "bias")

        self.masterbias = self.combining(bias) * self.MASK
        self.master_dict['bias'] = self.masterbias



    def do_masterflat(self):
        """Este método crea el MasterFlat para cada filtro.
        """
        lst_flat = [elem for elem in list(self.DATA_DICT.keys()) if 'flat' in elem]
        for filt in lst_flat:
            flat = self.sustractMasterBias(filt, self.masterbias, self.DATA_DICT)
            combflat = self.combining(flat)
            median = np.nanmedian(combflat)
            masterflat = combflat/median
            self.master_dict[filt] = masterflat



    def get_std(self, no_CRs=False, contrast_arg = 1.5, cr_threshold_arg = 5.,
                neighbor_threshold_arg = 5. ):
        """
        Este método limpia las imágenes de calibración fotométrica.
        """
        lst_std = [elem for elem in list(self.DATA_DICT.keys()) if 'std' in elem]
        for elem in lst_std:
            key, value = elem.split('+')
            args = [elem, self.master_dict['bias'],
                    self.master_dict['flat+' + value]] 
            std = self.clean_target(*args)
            lst_sd = []
            for sd in std:
                if no_CRs:
                    no_mask = np.nan_to_num(sd, nan=np.nanmedian(sd))
                    lst_sd.append(lacosmic.lacosmic(no_mask, contrast=contrast_arg,
                                                    cr_threshold=cr_threshold_arg,
                                                    neighbor_threshold=neighbor_threshold_arg,
                                                    effective_gain=1.9,
                                                    readnoise=4.3)[0])
                else:
                    lst_sd.append(np.nan_to_num(sd, nan=np.nanmedian(sd)))
            self.std_dict[elem] = lst_sd



    def get_target(self, no_CRs=False, contrast_arg = 1.5, cr_threshold_arg = 5.,
                neighbor_threshold_arg = 5.):
        """
        Este método limpia las imágenes de ciencia.
        """
        lst_target = [elem for elem in list(self.DATA_DICT.keys()) if 'target' in elem]
        for elem in lst_target:
            key, value = elem.split('+')
            args = [elem, self.master_dict['bias'],
                    self.master_dict['flat+' + value]] 
            target= self.clean_target(*args)
            lst_tg = []
            for tg in target:
                if no_CRs:
                    no_mask = np.nan_to_num(tg, nan=np.nanmedian(tg))
                    lst_tg.append(lacosmic.lacosmic(no_mask, contrast=contrast_arg,
                                                    cr_threshold=cr_threshold_arg,
                                                    neighbor_threshold=neighbor_threshold_arg,
                                                    effective_gain=1.9,
                                                    readnoise=4.3)[0])
                else:
                    lst_tg.append(np.nan_to_num(tg, nan=np.nanmedian(tg)))
            self.target_dict[elem] = lst_tg



    def remove_fringing(self):
        """Este método permite llevar a cabo una limpieza especial en el caso del uso de Sloan_z
        para eliminar el patrón de interferencia.
        """
        lst_results = [elem for elem in list(self.DATA_DICT.keys())]
        if 'Sloan_z' in lst_results:
            fringe= self.target_dict['target+Sloan_z']
            combfringe = self.combining(fringe)
            median = np.nanmedian(combfringe)
            masterfringe = combfringe/median
            fr_free = [elem/masterfringe for elem in fringe]
            self.target_dict['fringe+Sloan_z'] = fr_free

        
    
    def sustract_sky(self):
        """En esta instancia se procede a sustraer la contribución del fondo de cielo.
        """
        lst_target_keys = [elem for elem in list(self.target_dict.keys())]
        for elem in lst_target_keys:
            key, value = elem.split('+')
            lst_frames = self.target_dict['target+' + value]
            #if len(lst_frames) != 0:
            cube = np.dstack(lst_frames)
            cube.sort(axis=2)
            im_avg = np.median(cube[:,:,:], axis=2)
            self.bkg = sep.Background(im_avg)
            no_sky = []
            for fr in lst_frames:
                no_sky.append(fr-self.bkg)
            self.target_dict['sky+' + value] = no_sky
            #else:
            #    print(f"ERROR: Target list for {value} is empty!!!")
            #    continue



    def save_target(self, fringing=False, std=False, sky=False):
        """Este método permite salvar las imágenes generadas durante el proceso de limpieza.
        Además añade información al header que permitirá ayudar en proceso futuros.

        Args:
            fringing (bool, optional): Para indicar si se ha hecho corrección del patrón de interferencia.
            Defaults to False.
            std (bool, optional): Indicar si la/s imágen/es a salvar son imágenes de la estrella de
            calibración fotométrica. Defaults to False.
        """
        if not std:
            lst_results = [elem for elem in list(self.DATA_DICT.keys()) if 'target' in elem]
        else:
            lst_results = [elem for elem in list(self.DATA_DICT.keys()) if 'std' in elem]

        for key in lst_results:
            fnames = self.DATA_DICT[key]

            if key == 'target+Sloan_z' and fringing:
                adjetive = 'fringe'
                sky_status = 'SKY'
                status='REDUCED'
                imagetype='SCIENCE'
                target = self.target_dict['fringe+Sloan_z']
            elif std:
                adjetive = 'std'
                sky_status = 'SKY'
                status='REDUCED'
                imagetype='STD'
                filt = key.split('+')[1]
                target= self.std_dict[key]
            elif sky:
                adjetive = 'target'
                sky_status = 'SKY'
                status='REDUCED'
                imagetype='SCIENCE'
                filt = key.split('+')[1]
                target= self.target_dict[key]
            else:
                adjetive = 'target'
                sky_status = 'NOSKY'
                status='REDUCED'
                imagetype='SCIENCE'
                filt = key.split('+')[1]
                target= self.target_dict['sky+'+ filt]

            for i in range(len(fnames)):
                t = time.gmtime()
                time_string = time.strftime("%Y-%m-%dT%H:%M:%S", t)
                hd = fits.open(fnames[i])[0].header
                hd_wcs = wcs.WCS(fits.open(fnames[i])[0].header).to_header()
                hd['imgtype'] = imagetype
                hd['STATUS'] = status
                hd['SSKY'] = sky_status
                hd['BPMNAME'] = 'BPM_5sig' #Tipo de BPM aplicada.
                hd['rdate'] = time_string #Fecha de la reducción
                hd['filtro'] = filt #Simplifica la búsqueda por filtro.
                primary_hdu = fits.PrimaryHDU(target[i], header=hd+hd_wcs)
                hdul = fits.HDUList([primary_hdu])

                filename = os.path.basename(fnames[i])
                print(f"Storing the frame: reduced_5sig_{adjetive}_{sky_status}_{filename} for {hd['FILTER2']}")
                hdul.writeto(str(self.PATH_RESULTS / (f'reduced_5sig_{adjetive}_{sky_status}_' + filename)),
                            overwrite=True)