#### Retired code
    #############################################################
    #### Graveyard
    #############################################################
    def fixed_smooth_dEBV(self,final_psf=10.):
        """
        not to be used because this depends exactly on the choice of distance bins; it won't even run the way the class is set up now.
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
    def load_fixed_smooth_dEBV(self):
        """
        not to be used because this depends exactly on the choice of distance bins
        """
        filename = DUST_3D_TEMPERATURE_MAP_DATA_LOCATION +"/smoothed_maps/albert/smoothed_dEBV_"+str(int(self.bayestar_nside))+".hdf5"

        with h5py.File(filename, "r") as g:    
            self.fixed_dEBV_smooth=g["dEBV_array"][()]
            g.close()
        return self.fixed_dEBV_smooth
