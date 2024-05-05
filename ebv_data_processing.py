
from __future__ import division,print_function

import h5py
import numpy as np

from general_settings import Settings
settings = Settings()
DATA_LOCATION = settings.data_location

import utils

### Functions to load some of the reddening data

def load_bayestar_2019():
    """
    Load the bayestar 2019 data calculated from the medians of all the samples stored.

    Returns:
        numpy.ndarray: The loaded bayestar 2019 data with transposed axes, where the distance is the first axis and healpix is the second.
    """
    print("Loading the bayestar 2019 calculated from the medians of all the samples stored")
    filename =DATA_LOCATION +"/dust_reddening_maps/bayestar2019.hdf5"
    with h5py.File(filename, "r") as g: 
        bayestar_3D_map=g["EBV"][()]
        print("The type of the data was:",g["EBV"].dtype)
    print("Loading bayestar 2019 done.")
    return bayestar_3D_map.T

#### Distance Functions
def query_bayestar_distances(bayestar_version):
    """
    Query the Bayestar catalog for distances.
    
    Parameters:
        bayestar_version (str): The version of the Bayestar catalog to query.
    
    Returns:
        list: A list of distances from the Bayestar catalog.
    """
    bayestar = BayestarQuery(max_samples=2, version=bayestar_version)
    return bayestar.distances
def saving_distances():
    """
    Save distances obtained from querying bayestar distances.

    This function iterates over a list of bayestar versions and queries the distances using the
    `query_bayestar_distances` function. The distances are then saved to files using the `saveArray`
    function from the `utils` module.

    Parameters:
        None

    Returns:
        None
    """
    bayestar_version_list = ["bayestar2017","bayestar2019"]
    for bv in bayestar_version_list:
        distances = query_bayestar_distances(bv)
        utils.saveArray(np.array(distances), filename=DATA_LOCATION+"/dust_reddening_maps/distance_slices/"+bv+"_distances")
def load_bayestar_distances(bayestar_version='bayestar2019'):
    """
    Load distances obtained from querying bayestar distances.

    This function loads the distances from the file using the `openFits` function from the `utils` module.

    Parameters:
        bayestar_version (str): The bayestar version to load the distances for. Default is 'bayestar2019'.

    Returns:
        numpy.ndarray: The distances loaded from the file.
    """
    return utils.openFits(filename=DATA_LOCATION+"/dust_reddening_maps/distance_slices/"+bayestar_version+"_distances")

##### this function may not be in use in the code anymore, it's more for archival purposes
def loadBayestarSlices(dorebin=0):
    """
    Load the Bayestar Slices Data as were made by Albert Lee (Ioana made a newer version)

    Parameters:
    - dorebin (int): Flag indicating whether to rebin the data or not. Default is 0.

    Returns:
    - ebv: The loaded Bayestar Slices Data.
    """
    # Get Bayestar slices, as made by Albert
    print('loading the Bayestar Slices Data')
    if dorebin:
        ebv = utils.openFits(DATA_LOCATION+'/3D_dust_temperature/albert/maps_EBV_rebin-256.fits')
    else:
        ebv = utils.openFits(DATA_LOCATION+'/3D_dust_temperature/albert/maps_ebv-mean-full.fits')
        #ebv = utils.openFits('../data/maps_ebv.fits')
    return ebv