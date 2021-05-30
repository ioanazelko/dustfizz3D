import os
import sys

from configparser import ConfigParser                  ### For parsing the configuration file
from dustmaps.bayestar import BayestarQuery
import h5py
import numpy as np


DUST_3D_TEMPERATURE_MAP_DATA_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_DATA_LOCATION"]
DUST_3D_TEMPERATURE_MAP_CODE_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_CODE_LOCATION"]
DUST_3D_TEMPERATURE_MAP_PAPER_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_PAPER_LOCATION"]
DUST_3D_TEMPERATURE_MAP_PLOTS_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_PLOTS_LOCATION"]

sys.path.insert(0, DUST_3D_TEMPERATURE_MAP_CODE_LOCATION)




import ebv_data_processing
import emission_data_processing
import utils


class DataProcessing():
    def __init__(self,data_options_dict):
       
        self.bayestar_version = data_options_dict["bayestar_version"]
        self.bayestar_nside = data_options_dict["bayestar_nside"]
        self.planck_version=data_options_dict["planck_version"]
        self.planck_nside = data_options_dict["planck_nside"]
        # ### Load Data
        # self.load_data()
        # self.load_smooth_planck()
        # self.load_smooth_ebv()
        #self.load_fixed_smooth_dEBV()

        ### Save bayestar distances to file
        #self.saving_distances()
    
    def load_planck(self):
        self.planck = emission_data_processing.loadPlanckData() ### (5, 12582912) , (nr of frequencies, number of emission pixels), NESTED
        return self.planck
    def load_bayestar(self):
        self.ebv = ebv_data_processing.loadBayestarSlices()### (31,12582912) , (nr of slices, number of extinction pixels)
        return self.ebv

    #### Distance Functions
    def query_bayestar_distances(self, bayestar_version):
        bayestar = BayestarQuery(max_samples=2, version=bayestar_version)
        return bayestar.distances
    def saving_distances(self):
        bayestar_version_list = ["bayestar2017","bayestar2019"]
        for bv in bayestar_version_list:
            distances = self.query_bayestar_distances(bv)
            utils.saveArray(np.array(distances),filename=DUST_3D_TEMPERATURE_MAP_DATA_LOCATION+"/distance_slices/"+bv+"_distances")
    def load_distances(self):
        self.model_dist_slices = utils.openFits(filename=DUST_3D_TEMPERATURE_MAP_DATA_LOCATION+"/distance_slices/"+self.bayestar_version+"_distances")
        return  self.model_dist_slices  
    ### Smoothing Functions
    def smooth_planck(self,final_psf=10.):
        smoothing_psf_array=utils.calculate_smoothing_psf(data_type="Planck",target_psf=final_psf)
        smooth_final_maps = np.zeros(self.planck.shape)
        for freq_index in range(self.nfreq):
        #for freq_index in range(1):

            smoothing_psf_arcmin = smoothing_psf_array[freq_index] #arcmin
            smoothing_psf = utils.get_rad_from_arcmin(smoothing_psf_arcmin) #rad

            freq_str = str(int(self.freq_array[freq_index]))
            nested_original_map = self.planck[freq_index]
            ring_original_map=hp.reorder(map_in=nested_original_map,inp="NESTED",out='RING')
            ring_smooth_map = hp.smoothing(ring_original_map, fwhm=smoothing_psf) #NEED RING INPUT!
            nested_smooth_map = hp.reorder(map_in=ring_smooth_map,inp="RING",out='NESTED')
            smooth_final_maps[freq_index] = nested_smooth_map 
            ###Plotting the maps
            hp.mollview(nested_smooth_map, title="Smoothed Planck dust emission "+freq_str+"GHz", nest=True, min=0,max=15.)
            hp.gnomview(nested_original_map, title="Planck dust emission "+freq_str+"GHz", nest=True, min=0,max= 15.,rot=self.rot,xsize=self.xsize)
            hp.gnomview(nested_smooth_map, title="Smoothed Planck dust emission "+freq_str+"GHz", nest=True,min=0, max=15.,rot=self.rot,xsize=self.xsize)           
        ### Saving the smoothed maps    
        
        filename = DUST_3D_TEMPERATURE_MAP_DATA_LOCATION +"/smoothed_maps/albert/smoothed_planck_"+str(int(self.planck_nside))+".hdf5"
        with h5py.File(filename, "w") as f:
            f.create_dataset("planck_array",data=smooth_final_maps,dtype='float64')
            f.close()
    def load_smooth_planck(self):
        filename = DUST_3D_TEMPERATURE_MAP_DATA_LOCATION +"/smoothed_maps/albert/smoothed_planck_"+str(int(self.planck_nside))+".hdf5"
        with h5py.File(filename, "r") as g:    
            self.planck_smooth=g["planck_array"][()]
            g.close()
        return self.planck_smooth

    def load_fixed_smooth_dEBV(self):
        """
        not to be used because this depends exactly on the choice of distance bins
        """
        filename = DUST_3D_TEMPERATURE_MAP_DATA_LOCATION +"/smoothed_maps/albert/smoothed_dEBV_"+str(int(self.bayestar_nside))+".hdf5"

        with h5py.File(filename, "r") as g:    
            self.fixed_dEBV_smooth=g["dEBV_array"][()]
            g.close()
        return self.fixed_dEBV_smooth

    def load_smooth_ebv(self):
        filename = DUST_3D_TEMPERATURE_MAP_DATA_LOCATION +"/smoothed_maps/albert/smoothed_ebv_"+str(int(self.bayestar_nside))+".hdf5"

        with h5py.File(filename, "r") as g:    
            self.ebv_smooth=g["ebv_array"][()]
            g.close()
        return self.ebv_smooth
    def fixed_smooth_dEBV(self,final_psf=10.):
        """
        not to be used because this depends exactly on the choice of distance bins
        """
        #nest = True for nested, False for Ring. Default is False
        smooth_psf_arcmin = utils.calculate_smoothing_psf(data_type="Bayestar",target_psf=final_psf,NSIDE=self.bayestar_nside) #arcmin
        smoothing_psf = utils.get_rad_from_arcmin(smooth_psf_arcmin) #rad
        smooth_final_maps = np.zeros(self.dEBV.shape)

        for ds_index in range(self.model_nslices):
        #for ds_index in range(1):
            nested_original_map = self.dEBV[ds_index]
            ring_original_map=hp.reorder(map_in=nested_original_map,inp="NESTED",out='RING')
            ring_smooth_map = hp.smoothing(ring_original_map, fwhm=smoothing_psf) #NEEDS RING INPUT!
            nested_smooth_map = hp.reorder(map_in=ring_smooth_map,inp="RING",out='NESTED')
            smooth_final_maps[ds_index] = nested_smooth_map 
            ###Plotting the maps
            hp.mollview(nested_smooth_map, title="Smoothed differential E(B-V) at distance slice "+str(ds_index) +\
                         " at "+'{:.2f}'.format(self.model_dist_slices[ds_index])+" kpc", nest=True, max=1)
            hp.gnomview(nested_original_map, title="Original differential E(B-V) at distance slice "+str(ds_index) +\
                         " at "+'{:.2f}'.format(self.model_dist_slices[ds_index])+" kpc", nest=True, max=1, rot=self.rot,xsize=self.xsize)
            hp.gnomview(nested_smooth_map, title="Smoothed differential E(B-V) at distance slice "+str(ds_index) +\
                         " at "+'{:.2f}'.format(self.model_dist_slices[ds_index])+" kpc", nest=True, max=1, rot=self.rot,xsize=self.xsize)
        ### Saving the smoothed maps    
            
        filename = DUST_3D_TEMPERATURE_MAP_DATA_LOCATION +"/smoothed_maps/albert/smoothed_dEBV_"+str(int(self.bayestar_nside))+".hdf5"
        with h5py.File(filename, "w") as f:
            f.create_dataset("dEBV_array",data=smooth_final_maps,dtype='float64')
            f.close()

    def smooth_ebv(self,final_psf=10.):
        #nest = True for nested, False for Ring. Default is False
        smooth_psf_arcmin = utils.calculate_smoothing_psf(data_type="Bayestar",target_psf=final_psf,NSIDE=self.bayestar_nside) #arcmin
        smoothing_psf = utils.get_rad_from_arcmin(smooth_psf_arcmin) #rad
        smooth_final_maps = np.zeros(self.ebv.shape)

        for ds_index in range(self.bayestar_ndistances):
        #for ds_index in range(1):
            nested_original_map = self.ebv[ds_index]
            ring_original_map=hp.reorder(map_in=nested_original_map,inp="NESTED",out='RING')
            ring_smooth_map = hp.smoothing(ring_original_map, fwhm=smoothing_psf) #NEEDS RING INPUT!
            nested_smooth_map = hp.reorder(map_in=ring_smooth_map,inp="RING",out='NESTED')
            smooth_final_maps[ds_index] = nested_smooth_map 
            ###Plotting the maps
            hp.mollview(nested_smooth_map, title="Smoothed integral E(B-V) at distance slice "+str(ds_index) +\
                         " at "+'{:.2f}'.format(self.bayestar_distances[ds_index])+" kpc", nest=True, max=1)
            hp.gnomview(nested_original_map, title="Original integral E(B-V) at distance slice "+str(ds_index) +\
                         " at "+'{:.2f}'.format(self.bayestar_distances[ds_index])+" kpc", nest=True, max=1, rot=self.rot,xsize=self.xsize)
            hp.gnomview(nested_smooth_map, title="Smoothed integral E(B-V) at distance slice "+str(ds_index) +\
                         " at "+'{:.2f}'.format(self.bayestar_distances[ds_index])+" kpc", nest=True, max=1, rot=self.rot,xsize=self.xsize)
        ### Saving the smoothed maps    

        filename = DUST_3D_TEMPERATURE_MAP_DATA_LOCATION +"/smoothed_maps/albert/smoothed_ebv_"+str(int(self.bayestar_nside))+".hdf5"
        with h5py.File(filename, "w") as f:
            f.create_dataset("ebv_array",data=smooth_final_maps,dtype='float64')
            f.close()