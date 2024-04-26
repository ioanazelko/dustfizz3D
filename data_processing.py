import os

import h5py
import healpy as hp
import numpy as np
import time
import logging
### to save the logged status output of the code
log_file_name = os.path.join(os.path.dirname(__file__), 'app.log') #this adds the right slash for the OS being using
logging.basicConfig(level=logging.DEBUG, filename=log_file_name, filemode='w', format='%(name)s - %(levelname)s - %(message)s')


from general_settings import Settings
settings = Settings()
DATA_LOCATION = settings.data_location

import ebv_data_processing
import emission_data_processing
import utils


class DataProcessing():
    def __init__(self,data_options_dict):
        """
        Initializes the DataProcessing class.
        This class is used to load and smooth the Planck and Bayestar data.

        Args:
        - data_options_dict (dict): A dictionary containing the data options, 
        like the Bayestar version and nside, and the Planck version and nside.

        """
        self.bayestar_version = data_options_dict["bayestar_version"]
        self.bayestar_nside = data_options_dict["bayestar_nside"]
        self.planck_version=data_options_dict["planck_version"]
        self.planck_nside = data_options_dict["planck_nside"]
        self.freq_array = utils.dust_temperature_map_frequency_array()
        self.nfreq = len(self.freq_array)

        ########## Zooming in parameters
        # this is not necessarily the same as sky area, because if sky area is the full sky, we may still want to have a way to zoom in
        zoom_in_dict=utils.get_sky_area_zoom_in_parameters("rho_ophiuchi")
        self.rot = zoom_in_dict['rot']
        self.xsize= zoom_in_dict['xsize']

        # ### Load Data
        # self.load_data()
        # self.load_smooth_planck()
        # self.load_smooth_ebv()
        #self.load_fixed_smooth_dEBV()

        ### Save bayestar distances to file
        #self.saving_distances()
    
    def load_planck(self):
        """
        Loads the Planck data.

        Returns:
        - planck (numpy.ndarray): The loaded Planck data.

        """
        self.planck = emission_data_processing.loadPlanckData() ### (5, 12582912) , (nr of frequencies, number of emission pixels), NESTED
        return self.planck
    
    def load_bayestar(self):
        """
        Loads the Bayestar data.

        Returns:
        - ebv (numpy.ndarray): The loaded Bayestar data.

        """
        if self.bayestar_version == "bayestar2019":
            ### Load the full data set, indexed by (distance index, healpix index)

            self.ebv = ebv_data_processing.load_bayestar_2019()
        else:
            ### Load the deprecated map from Albert, not in use anymore

            self.ebv = ebv_data_processing.loadBayestarSlices()### (31,12582912) , (nr of slices, number of extinction pixels)
        return self.ebv

   
    def load_distances(self):
        """
        Loads the Bayestar distances.

        Returns:
        - model_dist_slices (numpy.ndarray): The loaded Bayestar distances.

        """
        self.model_dist_slices = ebv_data_processing.load_bayestar_distances(self.bayestar_version)
        return  self.model_dist_slices  
    
    ### Smoothing Functions
    def smooth_planck(self,final_psf=10.):
        """
        Smooths the Planck intensity data.

        Args:
        - final_psf (float): The final point spread function (PSF) value for the 

        """

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
        
        filename =DATA_LOCATION +"/3D_dust_temperature/smoothed_maps/albert/smoothed_planck_"+str(int(self.planck_nside))+".hdf5"
        with h5py.File(filename, "w") as f:
            f.create_dataset("planck_array",data=smooth_final_maps,dtype='float64')
            f.close()
        logging.info("Smoothed Planck data calculated and saved.")

    def load_smooth_planck(self):
        """
        Loads the smoothed Planck data.

        Returns:
        - planck_smooth (numpy.ndarray): The loaded smoothed Planck data.

        """
        filename =DATA_LOCATION +"/3D_dust_temperature/smoothed_maps/albert/smoothed_planck_"+str(int(self.planck_nside))+".hdf5"
        with h5py.File(filename, "r") as g:    
            self.planck_smooth=g["planck_array"][()]
            g.close()
        return self.planck_smooth



    def load_smooth_ebv(self):
        """
        Loads the smoothed Bayestar data.

        Returns:
        - ebv_smooth (numpy.ndarray): The loaded smoothed Bayestar data.

        """
        if self.bayestar_version=='albert':
            filename =DATA_LOCATION +"/3D_dust_temperature/smoothed_maps/albert/smoothed_ebv_"+str(int(self.bayestar_nside))+".hdf5"
        elif self.bayestar_version=="bayestar2019":
            filename =DATA_LOCATION +"/3D_dust_temperature/smoothed_maps/ioana/bayestar19_smoothed_ebv_"+str(int(self.bayestar_nside))+".hdf5"
        else:
            raise ValueError("wrong bayestar type")
        with h5py.File(filename, "r") as g:    
            self.ebv_smooth=g["ebv_array"][()]
            g.close()
        print("loading smooth ebv done")
      
        return self.ebv_smooth



    def smooth_ebv(self,final_psf=10.):
        """
        Smooths the Bayestar data.

        Args:
        - final_psf (float): The final point spread function (PSF) value.

        """
        #nest = True for nested, False for Ring. Default is False
        self.load_bayestar()
        smooth_psf_arcmin = utils.calculate_smoothing_psf(data_type="Bayestar",target_psf=final_psf,NSIDE=self.bayestar_nside) #arcmin
        smoothing_psf = utils.get_rad_from_arcmin(smooth_psf_arcmin) #rad
        smooth_final_maps = np.zeros(self.ebv.shape)
        self.bayestar_distances=self.load_distances()
        self.bayestar_ndistances=len(self.bayestar_distances)
        
        for ds_index in range(self.bayestar_ndistances):
        #for ds_index in range(1):
            print("I am doing distance slice ", str(ds_index))
            nested_original_map = self.ebv[ds_index]
            ring_original_map=hp.reorder(map_in=nested_original_map,inp="NESTED",out='RING')
            ring_smooth_map = hp.smoothing(ring_original_map, fwhm=smoothing_psf) #NEEDS RING INPUT!
            nested_smooth_map = hp.reorder(map_in=ring_smooth_map,inp="RING",out='NESTED')
            smooth_final_maps[ds_index] = nested_smooth_map 
        ### Saving the smoothed maps    
        if self.bayestar_version=='albert':
            filename =DATA_LOCATION +"/3D_dust_temperature/smoothed_maps/albert/smoothed_ebv_"+str(int(self.bayestar_nside))+".hdf5"
        elif self.bayestar_version=="bayestar2019":
            filename =DATA_LOCATION +"/3D_dust_temperature/smoothed_maps/ioana/bayestar19_smoothed_ebv_"+str(int(self.bayestar_nside))+".hdf5"
        else:
            raise ValueError("wrong bayestar type")
        with h5py.File(filename, "w") as f:
            f.create_dataset("ebv_array",data=smooth_final_maps,compression='lzf', chunks=True)
            f.close()
        logging.info("Smoothed reddening data calculated and saved.")


def make_smooth_ebv():
    ### Making the smooth ebv maps; this takes 1h 10 min , and needs a lot of RAM, so better be done with a machine that has 
    ### more than 50GB of RAM
    start_time = time.time()

    #### Making the smooth ebv maps
    data_proc_dict = {'bayestar_version':"bayestar2019",'bayestar_nside':1024,'planck_version':"albert",'planck_nside':1024}
    data_proc= DataProcessing(data_proc_dict)

    data_proc.smooth_ebv()
    ### calculate the total time this operation took
    time_string = utils.end_time(start_time)

def make_smooth_planck():
    start_time = time.time()
    
    #### Making the smooth planck maps
    data_proc_dict = {'bayestar_version':"bayestar2019",'bayestar_nside':1024,'planck_version':"albert",'planck_nside':1024}
    data_proc= DataProcessing(data_proc_dict)
    data_proc.load_planck()
    data_proc.smooth_planck()
    ### calculate the total time this operation took
    time_string = utils.end_time(start_time)

if __name__ == "__main__":
    pass
    #### Making the smooth ebv maps
    #make_smooth_ebv()
    ### Making the smooth planck maps
    #make_smooth_planck()