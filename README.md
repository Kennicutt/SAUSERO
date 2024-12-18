# SAUSERO

__SAUSERO__ is a reduction software for Broad Band Imaging mode of OSIRIS+.

Developed by __Fabricio M. Pérez-Toledo__

## General Description

SAUSERO processes raw science frames to address noise, cosmetic defects, and pixel heterogeneity, preparing them for photometric studies. These corrections are essential before any analysis can be performed. The operations applied to the images depend on the type of observation. This software has been specifically designed to reduce and prepare science frames for photometric studies.

### Key Reduction Steps:

1. Application of a __Bad Pixel Mask (BPM)__ to all frames.
2. Creation of the __Master Bias__.
3. Creation of the __Master Flat__.
4. Application of master calibration frames to both __standard star__ and __science frames__.
5. Removal of __cosmic rays__.
6. __Sky subtraction__.
7. Alignment of __science frames__.
8. __Astrometric calibration__.
9. __Flux calibration__.

### Input Requirements:

The software requires the following frames as input:

- Bias frames
- Sky flat frames
- Photometric standard star frames
- Science frames

## Outputs

The generated results consist of one image per observed band. For each image, the following corrections and calibrations will have been applied:

- Bias subtraction
- Flat-field correction (including fringing correction for the Sloan z band, if applicable)
- Image alignment and stacking
- Astrometric calibration
- Photometric calibration (estimation of the zero-point, ZP ± error)

To address cosmetic defects, a __Bad Pixel Mask (BPM)__ is applied, and the __LACosmic algorithm__ is used to handle cosmic ray removal.

## Requirements

### Operative System
- __Any__: The software is designed to run within a __Conda environment__, ensuring compatibility across platforms.

### Dependencies
The following Python packages are required (minimum versions specified):

-  astroalign>=2.4.1
-  astrometry_net_client>=0.3.0
-  astropy>=5.3.4
-  astroquery>=0.4.6
-  ccdproc>=2.4.1
-  lacosmic>=1.1.0
-  loguru>=0.7.2
-  matplotlib>=3.8.0
-  numpy>=1.25.2
-  PyYAML>=6.0.2
-  sep>=1.2.1

## Instalation

Installing SAUSERO is straightforward. Follow these steps:

1. __Activate your Conda environment__ (or create a new one if needed):

    conda activate <your_env>

2. __Install SAUSERO__ using pip:

    pip install sausero

That's it! SAUSERO is now ready to use.

## Usage

Using the same conda enviroment, we should command as the following example:

    sausero -pr <your_program> -bl <your_ob>

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