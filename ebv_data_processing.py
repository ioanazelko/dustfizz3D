

from __future__ import division,print_function

import os
import sys
import numpy as np



DUST_3D_TEMPERATURE_MAP_DATA_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_DATA_LOCATION"]
DUST_3D_TEMPERATURE_MAP_CODE_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_CODE_LOCATION"]
DUST_3D_TEMPERATURE_MAP_PAPER_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_PAPER_LOCATION"]
sys.path.insert(0, DUST_3D_TEMPERATURE_MAP_CODE_LOCATION)
import utils


def loadBayestarSlices(dorebin=0):
    # Get Bayestar slices
    print('loading the Bayestar Slices Data')
    if dorebin:
        ebv = utils.openFits(DUST_3D_TEMPERATURE_MAP_DATA_LOCATION+'/albert/maps_EBV_rebin-256.fits')
    else:
        ebv = utils.openFits(DUST_3D_TEMPERATURE_MAP_DATA_LOCATION+'/albert/maps_ebv-mean-full.fits')
        #ebv = utils.openFits('../data/maps_ebv.fits')
    return ebv
