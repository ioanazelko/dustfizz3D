# dustfizz3D


Code for creating the **3D  interstellar medium dust temperature map** from Zelko et al. 2022: https://arxiv.org/abs/2211.07667


We start from 3D dust reddening map. We can use different versions of maps available in the research community. For now, we have used the Green et al. 2019 version, and others will be added in the future.

## Introduction to git and GitHub

If this is your first time using git and GitHub, please:
- [ ] Read the Chapter 1 this free online book https://git-scm.com/book/en/v2/Getting-Started-About-Version-Control
- [ ] Follow the instructions here to clone this repository https://git-scm.com/book/en/v2/Getting-Started-About-Version-Control
- [ ] To make changes to the repository later on, read Chapters 2 and 3, and other chapters of interest.
## Installation Requirements

This project requires the following Python packages beyond the standard library:

- `corner`
- `healpy`
- `h5py`
- `scipy`
- `pandas`
- `ptemcee`
- `dustmaps`

Please ensure that you are using Python 3.x, as these packages are not compatible with Python 2.x.

### Installing the Requirements
Two popular tools for managing these dependencies are `pip` and `conda`. Understanding the differences between these tools can help you choose the right one for your needs.

**pip** is the Python Packaging Authority’s recommended tool for installing packages from the Python Package Index (PyPI). It is suitable for all Python environments and is great for installing Python packages that do not require complex dependencies outside the Python ecosystem.

**conda**, on the other hand, is a package manager from the Anaconda distribution. It excels in handling packages from both Python and other languages, managing complex dependencies, and creating isolated environments. Conda is particularly useful in data science and scientific computing where dependencies often include non-Python libraries.

For projects that require simple Python libraries, `pip` is usually sufficient. However, for projects that involve data science and need consistent handling of diverse dependencies across platforms, `conda` may be a better choice.

Conda's ability to create isolated environments can prevent conflicts between package versions and allow for more stable development platforms across different projects.

While Conda has a vast repository, it doesn't cover as many packages as pip. For some less common packages, users might still need to resort to pip in addition to conda.

I recommend you use conda where possible.

Note: the `dustmaps` package will prompt you to specify a configuration file upon use, to set up the path were the dust map data should be stored. Follow the instructions given in the message.
#### Installation with Conda

If you prefer to manage your Python environments and dependencies with Conda, you can follow these steps to set up your environment.

##### Setting Up the Conda Environment

First, if you do not have Conda installed, download and install it from Miniconda or Anaconda. Then, create a new Conda environment for this project:


`conda create -n myenv python=3.x

(replace x wtih the latest python version or the one desired)
Activate the new environment:

`conda activate myenv`

##### Installing Packages

Install the required packages in your Conda environment:

`conda install -c anaconda configparser` 

`conda install -c conda-forge corner healpy h5py scipy pandas ptemcee`

Note: Some packages might not be available in the default Anaconda channel and require installation from `conda-forge` or another channel.

##### Verifying Installation

To verify that all packages are installed correctly, you can list all packages in the active Conda environment:

`conda list`

Now, you are ready to use the project with all dependencies set up correctly.

#### Installation with pip
You can install all required packages using `pip`. Run the following command in your terminal:

`pip install corner dustmaps healpy h5py scipy pandas ptemcee`


### Setting up the location settings
Copy the file `general_settings_template.json` and rename it to `general_settings.json`.
Edit the file with the corresponding directory paths for where you would like data, plots, code and paper files to be stored. Please use absolute paths, not relative paths.


## Data


Download the data folder from https://www.dropbox.com/scl/fo/qsjiv5ejktmnyr14jc11s/AAUyIkXYCCnQkOWOujUYrRg?rlkey=fl9v14efr5irha85u26n0cg1w&st=yqxh9xfw&dl=0
and save it under a name of your choice.
Set the data location in `general_settings.json` to match to the folder path you chose above.


## Testing the installation


Note: please make sure you have at least 10GB of RAM.
Please run the script 
`sky_analysis.py`



-------------------------------



# Descriptions of each module:

![[sky_analysis_level_0_simple.svg]]


![[sky_analysis_level_0_simple.png]]

#### ebv_data_processing
Functions to load the pre-processed reddening data. 

#### emission_data_processing
Processes the emission data

#### data_processing
This has functionality to smooth both the reddening and the emission maps to match a desired  resolution

#### constant_configuration.cfg 
-- Holds the values for the general physics constants, so that they are not multiply defined in the code

#### optimizer.py

Functions that run the optimization process
#### model.py

Sets up the bayesian model
#### sky_analysis.py

Sets up the higher level part of the analysis with all the customized settings.
#### sky_plots.py

Functions that help in plotting. Ioana adds to this script
#### master_script.py

Showing what commands to call to reproduce some of the Figures of Zelko et al. 2022
(work in progress, Ioana keeps adding to this script)
#### plot_utils.py

Helpful functions for the plotting functions in sky_plots.

#### utils.py

General helpful functions for the analysis
