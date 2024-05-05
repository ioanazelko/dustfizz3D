from __future__ import division,print_function

import os
import sys
import numpy as np


from general_settings import Settings
import utils
settings = Settings()
DATA_LOCATION = settings.data_location


Kcmb2MJy = np.array([483.69, 287.45, 58.04, 2.27])
def convertBands(planck):
    print('Converting from Kcmb to MJySr-1...')
    planck[0,:] = planck[0,:]*Kcmb2MJy[0]
    planck[1,:] = planck[1,:]*Kcmb2MJy[1]
    return planck


def subtractCMB(planck,cmb):
    print('Subtracting CMB from Planck w approp. conversions from Kcmb to MJySr-1...')
    planck[0,:] = planck[0,:]-cmb
    planck[1,:] = planck[1,:]-cmb
    planck[2,:] = planck[2,:]-cmb*Kcmb2MJy[2]
    planck[3,:] = planck[3,:]-cmb*Kcmb2MJy[3]
    return planck



def loadPlanckData(dorebin=0):
    ### these might all be NEST ordered - yes they are

    # Get Planck data
    print("Reading map data...")

    if dorebin:
        planck = utils.openFits(DATA_LOCATION+'/3D_dust_temperature/albert/maps_planck2015-I_rebin-2048-256.fits', hdu=0)
        cmb = utils.openFits(DATA_LOCATION+'/3D_dust_temperature/albert/maps_CMB_rebin-1024-256.fits', hdu=0)
    else:
        # planck = utils.openFits(DATA_LOCATION+'/3D_dust_temperature/albert/maps_planck2015-I_rebin-2048-1024.fits', hdu=0)
        # cmb = utils.openFits(DATA_LOCATION+'/3D_dust_temperature/albert/COM_CMB_IQU-smica_1024_R2.02_full.fits', hdu=1)
        planck = utils.openFits(DATA_LOCATION+'/emission_maps/maps_planck2015-I_rebin-2048-1024.fits', hdu=0)
        cmb = utils.openFits(DATA_LOCATION+'/emission_maps/COM_CMB_IQU-smica_1024_R2.02_full.fits', hdu=1)
        cmb = cmb['i_stokes']
        cmb = cmb.reshape([1,cmb.shape[0]])

    planck = subtractCMB(planck,cmb)
    planck = convertBands(planck)

    # Add SFD data
    if dorebin:
        i100 = utils.openFits(DATA_LOCATION+'/3D_dust_temperature/albert/maps_i100_rebin-256.fits')
    else:
        ### This is NSIDE 1024
        # i100 = utils.openFits(DATA_LOCATION+'/3D_dust_temperature/albert/maps_i100.fits')       
        i100 = utils.openFits(DATA_LOCATION+'/emission_maps/maps_i100.fits')        
 
    planck = np.append(planck,i100.reshape([1,i100.shape[0]]),axis=0)

    return planck


