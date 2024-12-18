# SAUSERO

SAUSERO is a reduction software for Broad Band Imaging mode of OSIRIS+.

Developed by Fabricio M. Pérez-Toledo

## General Description

The science frames always present noises, cosmetic defects or heterogeneity between pixels (in other words raws)
and I have to perform several operations to correct them before the analising step. The operations on the images 
depend on the kind of observation. In this case, this software has been developed with the goal to prepare the 
science frames for photometric studies. The correction steps are:

1. Application of Bad Pixel Mask to all frames.
2. Creation of Master Bias.
3. Creation of Master Flat.
4. Application of masters to STD star and science frames.
5. Remove the cosmic rays.
6. Sky subtraction.
7. Aligning science frames.
8. Astrometrization.
9. Flux calibration.

It requires bias frames, skyflat frames, photometric standard star frames and science frames.

The generated results consist of one image per observed band. In each image, the following processes will have
been applied: bias subtraction, flat-field correction (and fringing correction in the case of Sloan z), image 
alignment and stacking, astrometric calibration, and photometric calibration (estimation of the ZP ± error). 
To remove cosmetic effects, a BPM is used, and the LACosmic algorithm is applied to handle cosmic rays.

## Requirements

Operative System: Any (it should run in conda enviroment).

Dependencies:

  astroalign>=2.4.1
  
  astrometry_net_client>=0.3.0
  
  astropy>=5.3.4
  
  astroquery>=0.4.6
  
  ccdproc>=2.4.1
  
  lacosmic>=1.1.0
  
  loguru>=0.7.2
  
  matplotlib>=3.8.0
  
  numpy>=1.25.2
  
  PyYAML>=6.0.2
  
  sep>=1.2.1

Hardware: RAM 4 GB

## Instalation

The installation is very easy. You only need to use pip as I show you below:

    conda activate <your_env>

    pip install sausero

## Usage

Using the same conda enviroment, we should command as the following example:

$ sausero -pr <your_program> -bl <your_ob>

where -pr is your GTC program indicator and -bl is the observed block number.

The first time, it will fail because it doesn't know which is the path to frames.
You have to modify the configuration file which is placed in your home directory
inside of a folder called as sausero. You need to fill some configuration parameters
as:

1. "PATH_DATA": "/path/to/your/frames/"
2. "No_Session":"astrometry api key"

The first one will be the root directory. You should have concern for the standard
structure. The directory with the frames must be as "<Your_Program>_<Your_OB>/".
Inside, it must create a raw directory where you will store the original frames.
During the execution, sausero will create other folder called as "reduced" where
the reduced frames will be saved.

The second one needs that you create an account in Astrometry.net. You will find
astrometry api key and paste it in configuration file.

After you finish to fill and save the configuration file, if you execute the code
again, it will run in this ocassion.

About the results, you will find the your frame directory a new collection of the 
frames. It will have:

A. Each science frame reduced, in other words, once the master bias and master
flat have been applied. Each frame is saved two time, one with sky and another
without sky.
B. Once aligned both cases.
C. They have after been applied the astrometrization.
D. A PNG file showing the sources in the FoV.
E. A PNG file showing the STD star.
F. Final reduced science frames with sky and without sky.

WARNING: By default, the internal configuration always avoid to share data with 
the astrometry.net community. So that's why your data are safe using
astrometry.net.

## Project Structure

    SAUSERO/
        BPM/
            BPM_OSIRIS_PLUS.fits -> BAD PIXEL MASK
        config/
            configuration.json   -> Configuration file.
        aligning_osirisplus.py   -> Aligns the science frames. 
        astrometry_osirisplus.py -> Astrometrization of the science frames.
        Color_Codes.py           -> Gives color to the comments
        OsirisDRP.py             -> Handles all the sofware and manages the frames. 
        photometry_osirisplus.py -> Carries out the photometric calibration.
        reduction_osirisplus.py  -> Carries out the clean process.

## Note about the frames

The code is designed to work with OSIRIS+ frames. They must be in FITS format.

## LICENSE

This software is under GPL v3.0 license. More information is available in the
repository.

## CONTACT

By e-mail: fabricio.perez@gtc.iac.es

Repository: https://github.com/Kennicutt/SAUSERO

SAUSERO

SAUSERO is a reduction software for the Broad Band Imaging mode of OSIRIS+.

Developed by Fabricio M. Pérez-Toledo.
General Description

SAUSERO processes raw science frames to address noise, cosmetic defects, and pixel heterogeneity, preparing them for photometric studies. These corrections are essential before any analysis can be performed. The operations applied to the images depend on the type of observation. This software has been specifically designed to reduce and prepare science frames for photometric studies.
Key Reduction Steps:

    Application of a Bad Pixel Mask (BPM) to all frames.
    Creation of the Master Bias.
    Creation of the Master Flat.
    Application of master calibration frames to both standard star and science frames.
    Removal of cosmic rays.
    Sky subtraction.
    Alignment of science frames.
    Astrometric calibration.
    Flux calibration.

Input Requirements:

The software requires the following frames as input:

    Bias frames
    Sky flat frames
    Photometric standard star frames
    Science frames

Outputs

The generated results consist of one image per observed band. For each image, the following corrections and calibrations will have been applied:

    Bias subtraction
    Flat-field correction (including fringing correction for the Sloan z band, if applicable)
    Image alignment and stacking
    Astrometric calibration
    Photometric calibration (estimation of the zero-point, ZP ± error)

To address cosmetic defects, a Bad Pixel Mask (BPM) is applied, and the LACosmic algorithm is used to handle cosmic ray removal.