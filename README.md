# 3D_dust_temperature_map

Code for creating the 3D  interstellar medium dust temperature map from Zelko et al. 2022: https://arxiv.org/abs/2211.07667


We start from 3D dust reddening map. We can use different versions of maps available in the research community. For now, we have used the Green et al. 2019 version, and others will be added in the future.


## Installation Requirements

This project requires the following Python packages beyond the standard library:

- `configparser`
- `corner`
- `datetime`
- `healpy`
- `h5py`
- `scipy.optimize`
- `pandas`
- `ptemcee`
- `str2bool`

### Installing the Requirements

You can install all required packages using `pip`. Run the following command in your terminal:

`pip install configparser datetime healpy h5py scipy str2bool pandas ptemcee`

Please ensure that you are using Python 3.x, as these packages are not compatible with Python 2.x.

## Installation with Conda

If you prefer to manage your Python environments and dependencies with Conda, you can follow these steps to set up your environment.

### Setting Up the Conda Environment

First, if you do not have Conda installed, download and install it from Miniconda or Anaconda. Then, create a new Conda environment for this project:


`conda create -n myenv python=3.x

(replace x wtih the latest python version or the one desired)
Activate the new environment:

`conda activate myenv`

### Installing Packages

Install the required packages in your Conda environment:

`conda install -c anaconda configparser` 

`conda install -c conda-forge corner healpy h5py scipy pandas ptemcee`

You may still need to install the `str2bool` package with pip using
`pip install str2bool`
Note: Some packages might not be available in the default Anaconda channel and require installation from `conda-forge` or another channel.

### Verifying Installation

To verify that all packages are installed correctly, you can list all packages in the active Conda environment:

`conda list`

Now, you are ready to use the project with all dependencies set up correctly.

### Setting up the location settings
Copy the file `general_settings_template.json` and rename it to `general_settings.json`.
Edit the file with the corresponding directory paths for where you would like data, plots, code and paper files to be stored.


### Testing the installation

Please run the script 
`sky_analysis.py`


Note to David:
Don't read below here.

-------------------------------

## Recreating the analysis of Zelko et al. 2022


See master_script.py



# Data

### Dust emission data
##### IRAS/DIRBE 100 mum (2998) GHz data:

we use the map created by [Schlegel98, SFD](https://ui.adsabs.harvard.edu/abs/1998ApJ...500..525S/abstract). The map is downloaded through the IDL interface.

The exact emission map that was used, before and after smoothing, can be found at:
(give dataverse link)


##### Planck

### Dust reddening data


Step 1, get the 3D reddening data, at the highest resolution.
Step 2. Decide the binning for it
In the verision from Zelko2022, nside 1024 was used, and 17 distance bins.

Step 3.
Smooth the map to match the resolution of the planck map. In fact, I had to smooth that map as well, bringing both to 10 arcminutes.


For Step 1, the code is in:

Smoothing the maps:



##### loading_all_3D_data.ipyn is a notebook showing how to load the end result 3D dust temperature map, with the other emission parameters. It also loads the reddening components, and they can all be combined to generate emission maps based on the fits.


Step 3: 


# Descriptions of each module:


#### ebv_data_processing
Functions to load the pre-processed reddening data. 


#### emission_data_processing
import utils


#### constant_configuration.cfg 
-- Holds the values for the general physics constants, so that they are not multiply defined in the code