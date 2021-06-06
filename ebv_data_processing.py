

from __future__ import division,print_function

import os
import sys
import h5py
import numpy as np



DUST_3D_TEMPERATURE_MAP_DATA_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_DATA_LOCATION"]
DUST_3D_TEMPERATURE_MAP_CODE_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_CODE_LOCATION"]
DUST_3D_TEMPERATURE_MAP_PAPER_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_PAPER_LOCATION"]
sys.path.insert(0, DUST_3D_TEMPERATURE_MAP_CODE_LOCATION)
import utils





def load_bayestar_2019():
    print("Loading the bayestar 2019 calculated from the medians of all the samples stored")
    filename = DUST_3D_TEMPERATURE_MAP_DATA_LOCATION +"/../dust_maps/ioana/bayestar2019.hdf5"
    with h5py.File(filename, "r") as g: 
        bayestar_3D_map=g["EBV"][()]
    print("Loading bayestar 2019 done")
    return bayestar_3D_map.T ## we transpose it because we want the distance to be the first axis and healpix the second

def loadBayestarSlices(dorebin=0):
    # Get Bayestar slices, as made by Albert
    print('loading the Bayestar Slices Data')
    if dorebin:
        ebv = utils.openFits(DUST_3D_TEMPERATURE_MAP_DATA_LOCATION+'/albert/maps_EBV_rebin-256.fits')
    else:
        ebv = utils.openFits(DUST_3D_TEMPERATURE_MAP_DATA_LOCATION+'/albert/maps_ebv-mean-full.fits')
        #ebv = utils.openFits('../data/maps_ebv.fits')
    return ebv

#### Distance Functions
def query_bayestar_distances(bayestar_version):
    bayestar = BayestarQuery(max_samples=2, version=bayestar_version)
    return bayestar.distances
def saving_distances():
    bayestar_version_list = ["bayestar2017","bayestar2019"]
    for bv in bayestar_version_list:
        distances = query_bayestar_distances(bv)
        utils.saveArray(np.array(distances),filename=DUST_3D_TEMPERATURE_MAP_DATA_LOCATION+"/distance_slices/"+bv+"_distances")
def load_bayestar_distances(bayestar_version='bayestar2019'):
	return utils.openFits(filename=DUST_3D_TEMPERATURE_MAP_DATA_LOCATION+"/distance_slices/"+bayestar_version+"_distances")