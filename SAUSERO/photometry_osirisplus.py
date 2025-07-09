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

from astropy.nddata import CCDData
import numpy as np
import matplotlib.pyplot as plt
import sep, json, os
from astroquery.simbad import Simbad
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
from astropy.visualization import LogStretch,imshow_norm, ZScaleInterval
from matplotlib.patches import Ellipse

from SAUSERO.Color_Codes import bcolors as bcl
from loguru import logger
import pkg_resources

extinction_dict = {
    'Sloan_u': [0.45, 0.02],
    'Sloan_g': [0.15, 0.02],
    'Sloan_r': [0.07, 0.01],
    'Sloan_i': [0.04, 0.01],
    'Sloan_z': [0.03, 0.01]
}

def readJSON_STD():
    """
    Reads the file containing the STD parameters.

    Returns:
        json: Collection of configuration parameters 
    """
    if os.path.exists(f"{os.path.expanduser('~')}/sausero/photometric_standards.json"):
        return json.load(open(f"{os.path.expanduser('~')}/sausero/photometric_standards.json"))
    else:
        std_path = pkg_resources.resource_filename(
            'SAUSERO', 'config/photometric_standards.json')
        #return json.load(open("SAUSERO/config/photometric_standards.json"))
        return json.load(open(std_path))
    


def get_ugriz(values):
    """From Sloan r apparent magnitude and a small list of the color,
    this function estimates the magnitude for each Sloan band.

    Args:
        values (list): List with Sloan r magnitude and colors. The order
        must be Sloan r, u-g , g-r, r-i, i-z. 

    Returns:
        float: Five magnitudes: Sloan u, Sloan g, Sloan r, Sloan i, Sloan z. 
    """
    r, ug, gr, ri, iz = values
    g = gr + r
    u = ug + g
    i = -ri + r
    z = -iz + i
    return u, g, r, i, z

bands_dict = readJSON_STD()

def photometry(programa, bloque, filename, conf, extinction_dict=extinction_dict, bands_dict=bands_dict):
    """This method estimates the instrumental magnitude (zeropoint) for the night, depending on the filter used.

    Args:
        programa (str): Science program code
        bloque (str): Observational block number assigned to a science program
        filename (str): Name of STD star file.

    Returns:
        float: Estimation of the instrumental magnitude and its error.
    """
    
    path = conf["DIRECTORIES"]["PATH_DATA"] + f"{programa}_{bloque}/reduced/"

    frame = CCDData.read(path + filename, unit='adu')
    logger.info("STD frame has been loaded successfully")

    hd = frame.header
    W = frame.wcs

    filtro = hd['FILTER2']
    #ind = index_dict[filtro]
    t = hd['EXPTIME']
    target_name = hd['OBJECT'].split('_')[1]
    
    ra, dec = bands_dict[target_name][:2]
    logger.info("Giving the coordinates for STD")

    c = SkyCoord(ra,dec, frame=FK5, unit=(u.hourangle, u.deg), obstime="J2000")

    x,y = W.world_to_pixel(c)
    logger.info("Transformed the original coordinates to pixel on the image.")

    fig = plt.figure(figsize=(15,20))
    ax = fig.add_subplot(1, 1, 1, projection=W)
    im, _ = imshow_norm(frame.data, ax, origin='lower',
                            interval=ZScaleInterval(),
                            stretch=LogStretch(a=1))
    plt.plot(x, y, 'xr')
    fig.colorbar(im)
    fig.savefig(path +f'STD_IN_FIELD-{filtro}.png')
    logger.info("STD's FoV has been saved as a PNG file")

    # Extract Sources

    frame_data = frame.data.byteswap().newbyteorder()

    bkg = sep.Background(frame_data)
    logger.info(f"Background estimated: {bkg.globalback} +- {bkg.globalrms}")

    kernel = np.array([[1., 2., 3., 2., 1.],
                   [2., 3., 5., 3., 2.],
                   [3., 5., 8., 5., 3.],
                   [2., 3., 5., 3., 2.],
                   [1., 2., 3., 2., 1.]])

    clean_data = frame_data - bkg
    logger.info("Background subtracted")

    objects = sep.extract(clean_data, conf["PHOTOMETRY"]["threshold"], 
                          err=bkg.globalrms, filter_kernel=kernel)
    
    logger.info("Objects catalogue in FoV has been created")

    # plot background-subtracted image
    fig = plt.figure(figsize=(15,20))
    ax = fig.add_subplot(1, 1, 1, projection=W)
    im, norm = imshow_norm(clean_data, ax, origin='lower',
                            interval=ZScaleInterval(),
                            stretch=LogStretch(a=1))

    # plot an ellipse for each object
    for i in range(len(objects)):
        e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                    width=6*objects['a'][i],
                    height=6*objects['b'][i],
                    angle=objects['theta'][i] * 180. / np.pi)
        e.set_facecolor('none')
        e.set_edgecolor('red')
        ax.add_artist(e)

    fig.savefig(path +f'SOURCES_DETECTED-{filtro}.png')
    logger.info("STD's FoV has been saved adding the objects")

    distance = np.sqrt((objects['x'] - x)**2 + (objects['y']-y)**2)

    index = np.argmin(distance)

    target = objects[index]

    X = np.array([objects['x'][index]])
    Y = np.array([objects['y'][index]])
    a = np.array([objects['a'][index]])
    b = np.array([objects['b'][index]])
    theta = np.array([objects['theta'][index]])
    flux = np.array([objects['flux'][index]])

    logger.info(f"STD name: {hd['OBJECT']}")
    logger.info(f"RA: {ra}, Dec: {dec}")
    logger.info(f"Exposure time: {t} sec")
    logger.info(f"Position: {X[0]}, {Y[0]}")
    logger.info(f"Ellipse info -> a: {a[0]}, b: {b[0]} & theta: {theta[0]}")
    logger.info(f"Flux: {flux[0]} counts")

    kronrad, krflag = sep.kron_radius(clean_data, X, Y, a, b, theta, 6.0)
    flux, fluxerr, flag = sep.sum_ellipse(clean_data, X, Y, a, b, theta, 2.5*kronrad,
                                        subpix=1)
    flag |= krflag  # combine flags into 'flag'

    logger.info(f"Kron Radius: {kronrad[0]} pxs")
    logger.info(f"AUTO FLUX: {flux[0]} counts")
    logger.info(f"AUTO ERROR FLUX: {fluxerr[0]} counts")
    logger.info(f"AUTO FLAG: {flag[0]}")

    #Estimation ZeroPoint
    mags = get_ugriz(bands_dict[target_name][2:])
    
    ZP = mags[2] + 2.5*np.log10(flux[0]/t) + (extinction_dict[filtro][0] * hd['AIRMASS'])

    dm = 0
    dflux = fluxerr[0]
    dZP = dm + dflux * (t/flux[0]) * np.log10(np.e) + (extinction_dict[filtro][1] * hd['AIRMASS'])

    #print(f"ZP value: {ZP} +- {dZP}")
    logger.info(f'Estimated ZP: {ZP} +- {dZP} for {filtro}')

    return ZP, dZP