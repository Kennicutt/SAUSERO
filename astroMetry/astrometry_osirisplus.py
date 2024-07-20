import yaml, glob

from astropy.nddata import CCDData
from astropy.wcs import WCS
from astrometry_net_client import Session, FileUpload, Settings


def settings(PATH_TO_CONFIG_FILE):
    """ Esta función abre el archivo que contiene los parámetros
    de configuración para una adecuada astrometrización de la imagen.

    Args:
        PATH_TO_CONFIG_FILE (str): Ruta al archivo YAML que contiene
        los valores para cada parámetro.

    Returns:
        dict: Diccionario con el contenido del archivo.
    """
    try:
        print('Loading settings for astrometrization...')
        with open(PATH_TO_CONFIG_FILE, 'r') as file:
            prime_service = yaml.safe_load(file)
    except FileNotFoundError:
        raise('Not found the configuration file.')

    return prime_service['OSIRIS']

def apply_astrometrynet_client(filename, conf):
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
    #Send the image
    s = Session(api_key=conf["ASTROMETRY"]["No_Session"])
    upl = FileUpload(filename, session=s, settings=ss)
    submission = upl.submit()
    submission.until_done()
    job = submission.jobs[0]
    job.until_done()
    if job.success():
        wcs = job.wcs_file()
    print(job.info())
    
    if wcs != None:
        return wcs

def modify_WCS(best_wcs, PATH_TO_FILE):
    """Esta función permite modificar WCS de la imagen de ciencia.

    Args:
        best_wcs (str): Nuevo WCS.
        PATH_TO_FILE (str): Ruta a la imagen de ciencia.

    Returns:
        CCDData: Imagen de ciencia con la nueva solución astrométrica.
    """
    frame = CCDData.read(PATH_TO_FILE, unit='adu')
    best_wcs = WCS(best_wcs)
    new_frame =  CCDData(data=frame.data, header=frame.header, wcs=best_wcs,
                         unit='adu')
    new_frame.write(PATH_TO_FILE, overwrite=True)
    return new_frame

def solving_astrometry(PRG, OB, filt, conf, sky, calib_std = False):
    """Esta función permite establecer en qué orden deben efectuarse los diferentes pasos
    para obtener la solución astrométrica.

    UPDATE (2024-01-15): Se añade líneas para calibrar la STD mejor.

    Args:
        PRG (str): Programa científico a observar.
        OB (str): Número del programa científico.
        filt (str): Filtro de interés.

    Returns:
        WCS, list, CCDData: Como resultado obtenemos la nueva WCS, una lista que contiene
        la posición de las estrellas en función de su nuevo WCS y la imagen de ciencia
        astrometrizada.
    """
    root_path = conf["DIRECTORIES"]["PATH_DATA"]

    if calib_std:
        LST_PATH_TO_FILE = glob.glob(root_path + PRG + '_' + OB + '/reduced/' + f'*std*.fits')
    else:
        LST_PATH_TO_FILE = [root_path + PRG + '_' + OB + '/reduced/' + f'aligned_result_{filt}_{sky}.fits']

    PATH_TO_FILE = LST_PATH_TO_FILE[0]

    best_wcs = apply_astrometrynet_client(PATH_TO_FILE, conf)

    new_frame = modify_WCS(best_wcs, PATH_TO_FILE)

    return best_wcs, new_frame #pixels
