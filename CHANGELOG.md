## [1.2.6] - 2026‑09‑11
- Update the version number in the pyproject.toml file and in the OsirisDRP.py file to 1.2.6.

## [1.2.5] - 2026‑09‑11
- Fixed a bug in the installation process related to pyproject.toml file.

## [1.2.4] - 2026‑09‑10
- Fixed a bug in relation to pyproject.toml file.

## [1.2.3] - 2026‑09‑10
- Fix a bug with the selection of science and STD images during the aligning. No afected the reduction, but it took exposure time from the STD images.
- Add a pyproject.toml file to the project to facilitate the installation of the package.

## [1.2.2] - 2026‑09‑09
- Change pkg_resources to importlib.resources to avoid deprecation warnings.
- Fix a bug with the exposure time in the header of the reduced files.
- Fix a bug with photometric calibration workflow when the user does not want to perform astrometric calibration or stacking.

## [1.2.1] - 2026‑01‑06
- Add a new file to describe new changes, corrections and improvements.
- An error is fixed regarding the inability to complete the reduction if it is indicated that astrometric calibration is not desired.
