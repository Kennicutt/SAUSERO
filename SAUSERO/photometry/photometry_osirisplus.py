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

from astropy.nddata import CCDData
import numpy as np
import matplotlib.pyplot as plt
import sep
from astroquery.simbad import Simbad
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
from astropy.visualization import LogStretch,imshow_norm, ZScaleInterval
from matplotlib.patches import Ellipse

from Color_Codes import bcolors as bcl
from loguru import logger

extinction_dict = {
    'Sloan_u': [0.45, 0.02],
    'Sloan_g': [0.15, 0.02],
    'Sloan_r': [0.07, 0.01],
    'Sloan_i': [0.04, 0.01],
    'Sloan_z': [0.03, 0.01]
}

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

def photometry(programa, bloque, filename, conf):
    """This method estimates the instrumental magnitude (zeropoint) for the night, depending on the filter used.

    Args:
        programa (str): Science program code
        bloque (str): Observational block number assigned to a science program
        filename (str): Name of STD star file.

    Returns:
        float: Estimation of the instrumental magnitude and its error.
    """
    # Dict with STD star collection: Coordinates, Sloan r magnitude and colors.
    bands_dict = {
        'G158-100': ['00 33 54.60', '-12 07 58.9', 14.691, 1.101, 0.510, 0.222, 0.092],
        'SA92-282': ['00 56 46.86',	'+00 38 30.9', 12.936, 1.000, 0.136, 0.021, -0.009],
        'Feige22': ['02 30 16.62',	'+05 15 50.6', 13.024, 0.050, -0.333, -0.303, -0.273],
        'SA95-193': ['03 53 20.59' , '+00 16 34.7', 13.844, 2.489, 1.097, 0.407, 0.214],
        'Ross49': ['05 44 56.81','+09 14 32.1',11.163, 1.130, 0.467, 0.162, 0.049],
        'Hilt_566': ['06 32 09.67','+03 34 44.3',10.787, 1.125, 0.673, 0.341, 0.211],
        'Ru-149F': ['07 24 14.02','-00 31 38.2', 13.119, 2.469, 0.867, 0.317, 0.166],
        'SA100-280': ['08 53 34.9','-00 37 38.5',11.689, 1.143, 0.308, 0.084, 0.003],
        'PG0918+029D': ['09 21 34.0','+02 46 39.0',11.937, 2.227, 0.817, 0.324, 0.166],
        'SA101-316': ['09 54 52.03','-00 18 34.5',11.438, 1.152, 0.309, 0.073, 0.007],
        'G162-66': ['10 33 42.81','-11 41 38.7', 13.227, -0.183, -0.387, -0.354, -0.303],
        'Feige_34': ['10 39 36.74','+43 06 09.2',11.423, -0.509, -0.508, -0.347, -0.265],
        'PG1047+003A': ['10 50 09.0','-00 01 08.0',13.303, 1.385, 0.519, 0.212, 0.087],
        'G163-50': ['11 07 59.97','-05 09 26.0', 13.266, 0.215, -0.277, -0.272, -0.271],
        'Feige66': ['12 37 23.52','+25 03 59.9',10.747, -0.345, -0.476, -0.367, -0.316],
        'SA104-428': ['12 41 41.29','-00 26 26.1',12.330, 2.153, 0.763, 0.279, 0.147],
        'PG1323-086D': ['13 26 05.26', '-08 50 35.7', 11.928, 1.210, 0.397, 0.132, 0.032],
        'Ross_838': ['14 01 44.48','+08 55 17.5',11.327, 1.277, 0.573, 0.239, 0.111],
        'PG1528+062B': ['15 30 39.55', '+06 01 13.1', 11.828, 1.235, 0.419, 0.143, 0.036],
        'G15-24': ['15 30 41.76','+08 23 40.4', 11.277, 1.035, 0.412, 0.151, 0.052],
        'BD+33 2642': ['15 51 59.89','+32 56 54.3',10.979, -0.018, -0.332, -0.284, -0.212],
        'Ross_530': ['16 19 51.65','+22 38 20.4',11.319, 1.273, 0.558, 0.229, 0.103],
        'Wolf_629': ['16 55 25.22','-08 19 21.3',11.129, 3.013, 1.413, 1.466, 0.648],
        'SA109-381': ['17 44 12.27','-00 20 32.8',11.514, 1.477, 0.547, 0.223, 0.094],
        'Ross_711': ['18 35 19.17','+28 41 55.4',11.295, 0.837, 0.282, 0.104, 0.015],
        'SA110-232': ['18 40 52.33', '+00 01 54.8', 12.287, 1.390, 0.552, 0.237, 0.094],
        'SA111-1925': ['19 37 28.64','+00 25 02.8',12.345, 1.397, 0.200, 0.061, 0.051],
        'Wolf_1346': ['20 34 21.88','+25 03 49.8',11.753, -0.016, -0.351, -0.309, -0.291],
        'SA112-805': ['20 42 46.75','+00 16 08.1',12.174, 1.183, -0.087, -0.135, -0.090],
        'BD+28 4211': ['21 51 11.02','+28 51 50.4',10.750, -0.517, -0.511, -0.379, -0.313],
        'G93-48': ['21 52 25.37','+02 23 19.6', 12.96, 10.107, -0.308, -0.307, -0.261],
        'SA114-656': ['22 41 35.07','+01 11 09.8',12.326, 1.961, 0.756, 0.293, 0.156],
        'GD_246': ['23 12 21.62','+10 47 04.3',13.346, -0.491, -0.504, -0.378, -0.367],
        'PG2336+004B': ['23 38 38.26','+00 42 46.4', 12.312, 1.101, 0.336, 0.100, 0.014]
    }
    
    # Extintion coefficients depends on the filter used.
    extinction_dict = {
        'Sloan_u': [0.45, 0.02],
        'Sloan_g': [0.15, 0.02],
        'Sloan_r': [0.07, 0.01],
        'Sloan_i': [0.04, 0.01],
        'Sloan_z': [0.03, 0.01]
    }
    
    #index_dict = {
    #    'Sloan_u': 0,
    #    'Sloan_g': 1,
    #    'Sloan_r': 2,
    #    'Sloan_i': 3,
    #    'Sloan_z': 4
    #}
    
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
    logger.info("Give the coordinates for STD")

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
    logger.info("STD's FoV has been saved like PNG file")

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
    logger.info(f"Exposure time: {t} sg")
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