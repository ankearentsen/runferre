# explore

Files you need to run LAMFER: 

- FERRE-compatible grid of synthetic spectra at LAMOST resolution, with teff, logg, feh and cfe ("gridfile")
- a folder with LAMOST spectra downloaded from the spectroscopy archive ("specfol")
- LAMOST summary .csv file downloaded from the spectroscopy archive, make sure to have added the gaia_source_id column ("data")

The code uses only standard packages (numpy, astropy, pandas, matplotlib), apart from pyraf if you want to do the optional step of correcting the radial velocities. In my experience this can be helpful for the LAMOST spectra. 

You also need to have FERRE installed (https://github.com/callendeprieto/ferre). 

**************

An example pipeline is provided in the lamfer.py file. 
