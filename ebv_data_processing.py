

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
    filename = DUST_3D_TEMPERATURE_MAP_DATA_LOCATION +"/../dust_maps/ioana/bayestar_test"+".hdf5"
    with h5py.File(filename, "r") as g: 
        bayestar_3D_map=g["EBV"][()]
    print("Loading done")
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
