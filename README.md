# SAUSERO

__SAUSERO__ is a reduction software for the Broad Band Imaging mode of OSIRIS+ at GTC.

Developed by __Fabricio M. Pérez-Toledo__

## General Description

**S**oftware to **AU**omatize in a **S**imple **E**nvironment the **R**eduction of **O**siris+ data (**SAUSERO**) processes OSIRIS+ raw science frames to address noise, cosmetic defects, and pixel heterogeneity, preparing them for photometric analysis. Correcting these artifacts is a critical prerequisite for reliable scientific analysis. The software applies observation-specific reduction steps, ensuring optimized treatment for different data types. Developed with a focus on simplicity and efficiency, **SAUSERO** streamlines the reduction pipeline, enabling researchers to obtain calibrated data ready for photometric studies.

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

- __Bias frames__
- __Sky flat frames__
- __Photometric standard star frames__
- __Science frames__

## Outputs

The generated results consist of one image per observed band. For each image, the following corrections and calibrations will have been applied:

- __Bias subtraction__
- __Flat-field correction__ (including fringing correction for the Sloan z band, if applicable)
- __Image alignment and stacking__
- __Astrometric calibration__
- __Photometric calibration__ (estimation of the zero-point, ZP ± error)

To address cosmetic defects, a __Bad Pixel Mask (BPM)__ is applied, and the __LACosmic algorithm__ is used to handle cosmic ray removal.

## Requirements

### Operative System
- __Any__: The software is designed to run within a __Conda environment__, ensuring compatibility across platforms.

### Dependencies
The following Python packages are required (minimum versions specified), however, they will be installed
automatically together the :

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
    sep>=1.2.1`

### Hardware Requirements
- __RAM__: Minimum 4GB (higher is recommended for large datasets).

## Installation

Installing SAUSERO is straightforward. Follow these steps:

1. __Activate your Conda environment__ (or create a new one if needed (see below)):
    ```
    conda activate <your_env>

2. __Install SAUSERO__ using `pip`:
    ```
    pip install sausero

That's it! SAUSERO is now almost ready to use ;)

### Optional: Creating a New Conda Environment

If you don’t have an existing Conda environment, you can create one specifically for SAUSERO with the following commands:

    conda create -n sausero_env python=3.9 -y
    conda activate sausero_env
    pip install sausero

## First-Time Setup

Once Conda is set up, you should run __SAUSERO__ for the first time to create the file `configuration.json` that has to be configured.

    sausero -pr <your_program> -bl <your_ob>

- `-pr`: Your GTC program indicator.
- `-bl`: The observed block number.

**ATTENTION**: The first time, the code will 'fail' because the configuration file does not know the root 
directory where the images are stored and your astrometry-api-key. To fix this, follow the instructions below.

You must edit the configuration file, which is located in your home directory inside 
a folder named `sausero/`.

You need to set the following parameters in the configuration file:

1. `PATH_DATA`: Set this to the root directory containing your frames. Example:

    ```
    "PATH_DATA": "/path/to/your/frames/"

The directory structure must follow the format `<Your_Program>_<Your_OB>/`. Inside this directory, you should have 
a `raw/` folder where the original frames are stored. During execution, __SAUSERO__ will create a new folder named 
`reduced/`, where the reduced frames will be saved.

2. `No_Session`: This is your Astrometry.net API key. Example:

    ```
    "No_Session":"astrometry-api-key"

To obtain this key, create an account on [Astrometry.net](https://nova.astrometry.net/). Copy your API key and paste it into the configuration file.

## Running SAUSERO

After updating and saving the configuration file, you can run the command again. This time, the software will execute successfully.

    sausero -pr <your_program> -bl <your_ob>

### Outputs and Results

Once the process is complete, you will find a collection of reduced frames in the `reduced/` folder inside your frame 
directory. The output includes:

A. __Reduced science frames__:
- One version with the sky included.
- One version with the sky subtracted.

B. __Aligned frames__:
- Both sky-included and sky-subtracted versions.

C. __Astrometrized frames__:
- Frames with astrometric calibration applied.

D. __Visualization PNG files__:
- A PNG showing the detected sources in the Field of View (FoV).
- A PNG showing the photometric standard star.

E. __Final reduced science frames__:
- Both sky-included and sky-subtracted versions.


### Important Notes

- By default, __SAUSERO__ ensures your data remains private when using Astrometry.net. The software's internal configuration avoids sharing any data with the Astrometry.net community, ensuring your data's security.

## Project Structure

    SAUSERO/
        BPM/
            BPM_OSIRIS_PLUS.fits -> BAD PIXEL MASK
        config/
            configuration.json   -> Configuration file.
        check_files.py           -> It determines which steps can be performed by the pipeline based on the available FITS files.
        aligning_osirisplus.py   -> Aligns the science frames. 
        astrometry_osirisplus.py -> Astrometrization of the science frames.
        Color_Codes.py           -> Gives color to the comments
        OsirisDRP.py             -> Handles all the sofware and manages the frames. 
        photometry_osirisplus.py -> Carries out the photometric calibration.
        reduction_osirisplus.py  -> Carries out the clean process.

## Note about the frames

The code is designed to work with __OSIRIS+__ frames. They must be in __FITS__ format.

## LICENSE

This software is under __GPL v3.0__ license. More information is available in the
repository.

## CONTACT

- __Email__: [fabricio.perez@gtc.iac.es](fabricio.perez@gtc.iac.es)

- __Repository__: [https://github.com/Kennicutt/SAUSERO](https://github.com/Kennicutt/SAUSERO)