

from __future__ import division,print_function

import os
import sys
import healpy as hp
import matplotlib as mpl             ### For plotting options
import matplotlib.pyplot as plt
import numpy as np
import time
  

DUST_3D_TEMPERATURE_MAP_DATA_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_DATA_LOCATION"]
DUST_3D_TEMPERATURE_MAP_CODE_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_CODE_LOCATION"]
DUST_3D_TEMPERATURE_MAP_PAPER_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_PAPER_LOCATION"]
DUST_3D_TEMPERATURE_MAP_PLOTS_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_PLOTS_LOCATION"]
sys.path.insert(0, DUST_3D_TEMPERATURE_MAP_CODE_LOCATION)

import utils
import sky_analysis


class SkyPlots(sky_analysis.SkyAnalysis):
    def __init__(self, run_name, run_type):
        sky_analysis.SkyAnalysis.__init__(self, run_name, run_type)
   
    def plot_in_3D(self):
        ## there is a python package that plots 3D, I can use it if I want to
        g_xyz=utils.openFits(DUSTY_DATA_LOCATION+'/3D_dust_temperature/greg_fits/green_cartesian.fits', hdu=0)
    def healpy_test(self):
        for NSIDE in [1,2,4,8,16,32,64,128,256,512,1024]:
            NPIX = hp.nside2npix(NSIDE)
            print(NPIX)
            m = np.arange(NPIX)
            hp.mollview(m, title="Mollview image RING, NSIDE " +str(NSIDE),nest=False)
            hp.graticule()
            plt.savefig(DUST_3D_TEMPERATURE_MAP_CODE_LOCATION+"/../presentation/healpix_ring"+str(NSIDE)+".jpg")

            hp.mollview(m, title="Mollview image NESTED, NSIDE "+str(NSIDE), nest=True)
            hp.graticule()
            plt.savefig(DUST_3D_TEMPERATURE_MAP_CODE_LOCATION+"/../presentation/healpix_nested"+str(NSIDE)+".jpg")
        return
    def plot_healpix_mollview(self,data,pixel_index_array,total_sky_pixels,title,min=None,max=None,nest=True):
        plot_array = np.zeros(total_sky_pixels)
        plot_array[pixel_index_array]=data
        hp.mollview(plot_array,title=title,nest=nest,min=min,max=max)


    def plot_healpix_gnomview(self,data,pixel_index_array,total_sky_pixels,title,min=None,max=None,nest=True,rot=(0,20),pixels=500):
        plot_array = np.zeros(total_sky_pixels)
        plot_array[pixel_index_array]=data
        hp.gnomview(plot_array,title=title,nest=nest,min=min,max=max,rot=rot,xsize=pixels)
    ################################
    #### Planck functions
    ################################
    def planck_data_analytics(self):
        print("The maximum value in each Planck frequency map is ",np.max(self.planck,axis=1)," MJy/sr" )
        print("The minumum value in each Planck frequency map is ",np.min(self.planck,axis=1)," MJy/sr" )
        print("The median value in each Planck frequency map is ",np.median(self.planck,axis=1)," MJy/sr")
        print("The number of negative values in each freq map is ",np.sum(self.planck[:]<0,axis=1) )
        print("The total number of pixels in each map is ", len(self.planck[0]))
   

    def plot_planck(self,full_sky=True,min=0.,max= 15.):
        #nest = True for nested, False for Ring. Default is False
        for freq_index in range(len(self.freq_array)):
            freq_str = str(int(self.freq_array[freq_index]))
            if full_sky == True:
                ###Plotting the positions of the 0s
                #hp.mollview(self.planck[freq_index]<0, title="Planck dust emission "+str(self.freq_array[freq_index])+"GHz", nest=True, max= 1)
                
                hp.mollview(self.planck[freq_index], title="Dust emission at "+freq_str+"GHz", nest=True,min=min, max= max, unit='MJySr-1',xsize=1000)
                plt.savefig(DUST_3D_TEMPERATURE_MAP_PLOTS_LOCATION+"/planck_"+freq_str+".pdf",dpi=400)
                plt.savefig(self.optimizer_plots_folder+"/"+freq_str+"/planck_"+\
                        freq_str+".jpg",dpi=400)
            else:
                hp.gnomview(self.planck[freq_index], title="Dust emission at "+freq_str+"GHz", nest=True, min=min, max= max,rot=self.rot,xsize=self.xsize, unit='MJySr-1')
                plt.savefig(DUST_3D_TEMPERATURE_MAP_PLOTS_LOCATION+"/planck_zoom_"+str(self.xsize)+"_"+freq_str+".jpg")
                plt.savefig(self.optimizer_plots_folder+"/"+freq_str+"/planck_zoom_"+str(self.xsize)+"_"+freq_str+".jpg")
                plt.savefig(self.optimizer_plots_folder+"/planck_zoom_"+str(self.xsize)+"_"+freq_str+".jpg")


    def plot_smooth_planck(self):
        for freq_index in range(self.nfreq):
            freq_str = str(int(self.freq_array[freq_index]))

            hp.mollview(self.planck_smooth[freq_index], title="Smoothed Planck dust emission at "+freq_str+"GHz", nest=True,min=0, max=15., unit='MJySr-1')
            filename="planck_"+freq_str+"_smooth_"+str(int(self.full_maps_nside))+".jpg"
            plt.savefig(DUST_3D_TEMPERATURE_MAP_CODE_LOCATION+"/../presentation/"+filename)
            plt.savefig(self.optimizer_plots_folder+"/"+freq_str+"/"+filename)
            hp.gnomview(self.planck[freq_index], title="Planck dust emission at "+freq_str+"GHz", nest=True, min=0,max= 15.,rot=self.rot,xsize=self.xsize, unit='MJySr-1')
            zoom_original_filename ="planck_zoom_"+str(self.xsize)+"_"+freq_str+"_"+str(int(self.full_maps_nside))+".jpg"
            plt.savefig(DUST_3D_TEMPERATURE_MAP_CODE_LOCATION+"/../presentation/"+zoom_original_filename)
            plt.savefig(self.optimizer_plots_folder+"/"+freq_str+"/"+zoom_original_filename)
            plt.savefig(self.optimizer_plots_folder+"/"+zoom_original_filename)
            hp.gnomview(self.planck_smooth[freq_index], title="Smoothed Planck dust emission at "+freq_str+"GHz", nest=True, min=0,max=15.,rot=self.rot,xsize=self.xsize, unit='MJySr-1')
            zoom_filename="planck_zoom_"+str(self.xsize)+"_"+freq_str+"_smooth_"+str(int(self.full_maps_nside))+".jpg"
            plt.savefig(DUST_3D_TEMPERATURE_MAP_CODE_LOCATION+"/../presentation/"+zoom_filename)
            plt.savefig(self.optimizer_plots_folder+"/"+freq_str+"/"+zoom_filename)
            plt.savefig(self.optimizer_plots_folder+"/"+zoom_filename)
    def plot_planck_histograms(self):
        for i in range(self.nfreq):
            plt.hist(self.planck[i],bins=500, range=(-0.1,2),log=True)
            plt.xlabel("Inu [MJy/sr]")
    ################################
    ### Bayestar functions
    ################################

    def plot_ebv(self):
        #nest = True for nested, False for Ring. Default is False
        for ds_index in range(self.bayestar_ndistances):
            hp.mollview(self.ebv[ds_index], title="Cummulative E(B-V) till distance slice "+str(ds_index) +\
                        " at "+'{:.2f}'.format(self.bayestar_distances[ds_index])+" kpc", nest=True, max=1)
    def plot_dEBV(self):
        #nest = True for nested, False for Ring. Default is False
        for ds_index in range(self.model_nslices):
            hp.mollview(self.dEBV[ds_index], title="Differential E(B-V) at distance slice "+str(ds_index) +\
                         " at "+'{:.2f}'.format(self.model_dist_slices[ds_index])+" kpc", nest=True, max=1)
    def plot_ebv_hist(self):    
        plt.hist(self.ebv[-1,:],bins=500);
        plt.hist(self.ebv[0,:],bins=500);
        plt.xlim(0,1)

    def plot_fixed_smooth_dEBV(self):
        """
        to be used for the case that smoothing was applied to dEBV
        """
        for ds_index in range(self.model_nslices):
            nested_map = self.fixed_dEBV_smooth[ds_index]
            not_smooth_title = "Differential E(B-V) at distance slice "+str(ds_index) +" at "+'{:.2f}'.format(self.model_dist_slices[ds_index])+" kpc"
            title="Fixed Smoothed Differential E(B-V) at distance slice "+str(ds_index) +" at "+'{:.2f}'.format(self.model_dist_slices[ds_index])+" kpc"
            hp.mollview(nested_map, title=title, nest=True, max=1.)
            filename="dEBV_smooth_slice_"+str(ds_index)+".jpg"
    #         plt.savefig(DUST_3D_TEMPERATURE_MAP_CODE_LOCATION+"/../presentation/"+filename)
    #         plt.savefig(self.optimizer_plots_folder+"/"+freq_str+"/"+filename)
            hp.gnomview(self.dEBV[ds_index], title=not_smooth_title, nest=True, max= 1.,rot=self.rot,xsize=self.xsize)

            hp.gnomview(self.fixed_dEBV_smooth[ds_index], title=title, nest=True, max=1.,rot=self.rot,xsize=self.xsize)
            zoom_filename="dEBV_smooth_zoom_"+str(self.xsize)+"_slice_"+str(ds_index)+"_fixed.jpg"
            #plt.savefig(DUST_3D_TEMPERATURE_MAP_CODE_LOCATION+"/../presentation/"+zoom_filename)
            plt.savefig(self.optimizer_plots_folder+"/"+zoom_filename)

    def plot_smooth_dEBV(self):
        """
        to be used for the case that smoothing was applied to EBV
        """
        for ds_index in range(self.model_nslices):
            nested_map = self.dEBV_smooth[ds_index]
            not_smooth_title = "Differential E(B-V) at distance slice "+str(ds_index) +" at "+'{:.2f}'.format(self.model_dist_slices[ds_index])+" kpc"
            title="Smoothed Differential E(B-V) at distance slice "+str(ds_index) +" at "+'{:.2f}'.format(self.model_dist_slices[ds_index])+" kpc"
            hp.mollview(nested_map, title=title, nest=True, max=1.)
            filename="dEBV_smooth_slice_"+str(ds_index)+".jpg"
    #         plt.savefig(DUST_3D_TEMPERATURE_MAP_CODE_LOCATION+"/../presentation/"+filename)
    #         plt.savefig(self.optimizer_plots_folder+"/"+freq_str+"/"+filename)
            hp.gnomview(self.dEBV[ds_index], title=not_smooth_title, nest=True, max= 1.,rot=self.rot,xsize=self.xsize)

            hp.gnomview(self.dEBV_smooth[ds_index], title=title, nest=True, max=1.,rot=self.rot,xsize=self.xsize)
            zoom_filename="dEBV_smooth_zoom_"+str(self.xsize)+"_slice_"+str(ds_index)+".jpg"
            #plt.savefig(DUST_3D_TEMPERATURE_MAP_CODE_LOCATION+"/../presentation/"+zoom_filename)
            plt.savefig(self.optimizer_plots_folder+"/"+zoom_filename)
    #################################
    ### Optimizer functions
    #################################

    
    def plot_final_optimized_functions(self,data_dict):
        optimized_function_array = data_dict["final_optimized_functions_array"]
        final_chi_square_array = data_dict["final_chi_square_array"]
        super_pixels_index_array = data_dict["super_pixels_index_array"]
        self.plot_healpix_mollview(optimized_function_array,super_pixels_index_array,self.nr_of_super_pixels,\
                              title=r"Optimized function $f$",max=100)
        plt.savefig(self.optimizer_plots_folder+"/final_optimized_functions.jpg")

        self.plot_healpix_mollview(final_chi_square_array,super_pixels_index_array,self.nr_of_super_pixels,\
                              title=r"$\chi^2$",max=100)
        plt.savefig(self.optimizer_plots_folder+"/final_chi_square.jpg")

    def plot_optimizer_sky_parameters(self,data_dict):
        parameters = data_dict["final_parameters_array"]
        super_pixels_index_array = data_dict["super_pixels_index_array"]
        offsets,rhos,betas,Ts = self.separate_sky_optimizer_parameters(parameters)
        ##### Plotting the offsets
        if self.fit_for_the_offsets == True:
            for freq_index in range(self.nfreq):
                freq_str = str(int(self.freq_array[freq_index]))

                self.plot_healpix_mollview(offsets[:,freq_index],super_pixels_index_array,self.nr_of_super_pixels,\
                                      title=r"Offset "+freq_str+" GHz",max=50)
                plt.savefig(self.optimizer_plots_folder+"/offset_"+freq_str+".jpg")
                self.plot_healpix_gnomview(offsets[:,freq_index],super_pixels_index_array,self.nr_of_super_pixels,\
                                      title=r"Offset "+freq_str+" GHz",max=50,rot=self.rot,pixels=self.xsize)
                plt.savefig(self.optimizer_plots_folder+"/"+"offset_zoom_"+freq_str+".jpg")
        #####  Plotting the rhos
        ## case for rho varying in each voxel
        if self.fixed_rho_along_sightline == False:
            for ds_index in range(self.model_nslices):
                self.plot_healpix_mollview(rhos[:,ds_index],super_pixels_index_array,self.nr_of_super_pixels,\
                                      title=r"$\rho$ at distance slice "+str(ds_index) +\
                                      " at "+'{:.2f}'.format(self.model_dist_slices[ds_index])+" kpc",max=2E-4)
                plt.savefig(self.optimizer_plots_folder+"/rho_at_distance_slice_"+str(ds_index)+".jpg")
                self.plot_healpix_gnomview(rhos[:,ds_index],super_pixels_index_array,self.nr_of_super_pixels,\
                                      title=r"$\rho$ at distance slice "+str(ds_index) +\
                                      " at "+'{:.2f}'.format(self.model_dist_slices[ds_index])+" kpc",max=2E-4,rot=self.rot,pixels=self.xsize)
                plt.savefig(self.optimizer_plots_folder+"/rho_zoom_at_distance_slice_"+str(ds_index)+".jpg")
        ## case for rho only varying in each superpixel, or fixed across the sky
        else:                          
            self.plot_healpix_mollview(rhos,super_pixels_index_array,self.nr_of_super_pixels,\
                                      title=r"$\rho$",max=2E-4)
            plt.savefig(self.optimizer_plots_folder+"/rho.jpg")
            self.plot_healpix_gnomview(rhos,super_pixels_index_array,self.nr_of_super_pixels,\
                                      title=r"$\rho$",max=2E-4,rot=self.rot,pixels=self.xsize)
            plt.savefig(self.optimizer_plots_folder+"/rho_zoom.jpg")
        ##### Plotting the betas
        ## case for beta varying in each voxel
        if self.fixed_beta_along_sightline == False:
            for ds_index in range(self.model_nslices):
                self.plot_healpix_mollview(betas[:,ds_index],super_pixels_index_array,self.nr_of_super_pixels,\
                                      title=r"$\beta$ at distance slice "+str(ds_index) +\
                                      " at "+'{:.2f}'.format(self.model_dist_slices[ds_index])+" kpc",max=2.5)
                plt.savefig(self.optimizer_plots_folder+"/beta_at_distance_slice_"+str(ds_index)+".jpg")
                self.plot_healpix_gnomview(betas[:,ds_index],super_pixels_index_array,self.nr_of_super_pixels,\
                                      title=r"$\beta$ at distance slice "+str(ds_index) +\
                                      " at "+'{:.2f}'.format(self.model_dist_slices[ds_index])+" kpc",max=2.5,rot=self.rot,pixels=self.xsize)
                plt.savefig(self.optimizer_plots_folder+"/beta_zoom_at_distance_slice_"+str(ds_index)+".jpg")

        ## case for beta only varying in each superpixel, or fixed across the sky
        else:
            self.plot_healpix_mollview(betas,super_pixels_index_array,self.nr_of_super_pixels,\
                                       title=r"$\beta$",max=2.5)
            plt.savefig(self.optimizer_plots_folder+"/beta.jpg")
            self.plot_healpix_gnomview(betas,super_pixels_index_array,self.nr_of_super_pixels,\
                                       title=r"$\beta$",max=2.5,rot=self.rot,pixels=self.xsize)
            plt.savefig(self.optimizer_plots_folder+"/beta_zoom.jpg")
        ##### Plotting the Ts
        ## case for T varying in each voxel
        if self.fixed_T_along_sightline == False: 
            for ds_index in range(self.model_nslices):
                self.plot_healpix_mollview(Ts[:,ds_index],super_pixels_index_array,self.nr_of_super_pixels,\
                                      title=r"$T$ at distance slice "+str(ds_index) +\
                                      " at "+'{:.2f}'.format(self.model_dist_slices[ds_index])+" kpc",min=10,max=25, unit='K')
                plt.savefig(self.optimizer_plots_folder+"/T_at_distance_slice_"+str(ds_index)+".jpg")
                self.plot_healpix_gnomview(Ts[:,ds_index],super_pixels_index_array,self.nr_of_super_pixels,\
                                      title=r"$T$ at distance slice "+str(ds_index) +\
                                      " at "+'{:.2f}'.format(self.model_dist_slices[ds_index])+" kpc",min=10,max=25,rot=self.rot,pixels=self.xsize,unit='K')
                plt.savefig(self.optimizer_plots_folder+"/T_zoom_at_distance_slice_"+str(ds_index)+".jpg")

        ## case for T only varying in each superpixel
        else:
            self.plot_healpix_mollview(Ts,super_pixels_index_array,self.nr_of_super_pixels,\
                                      title=r"$T$",max=50,unit='K')
            plt.savefig(self.optimizer_plots_folder+"/T.jpg")
            self.plot_healpix_gnomview(Ts,super_pixels_index_array,self.nr_of_super_pixels,\
                                      title=r"$T$",max=50,rot=self.rot,pixels=self.xsize,unit='K')
            plt.savefig(self.optimizer_plots_folder+"/T_zoom.jpg")


    def plot_reconstructed_total_emission(self,data_dict):
        total_emission_array =data_dict["total_emission_array"]
        full_resolution_pixel_index_array = data_dict["full_resolution_pixel_index_array"]
        for freq_index in range(self.nfreq):
            freq_str = str(int(self.freq_array[freq_index]))
            title = "Total Reconstructed Emission at "+ freq_str+" GHz"
            print("Title is",title)
            self.plot_healpix_mollview(total_emission_array[:,freq_index],full_resolution_pixel_index_array,\
                                    self.nr_of_super_pixels*self.super_pixel_size,title=title,min=0,max=15,unit='MJySr-1')
            plt.savefig(self.optimizer_plots_folder+"/total_reconstructed_emisssion_"+freq_str+".jpg")
            plt.savefig(self.optimizer_plots_folder+"/"+freq_str+"/total_reconstructed_emisssion_"+\
                        freq_str+".jpg")
            plt.savefig(self.optimizer_plots_folder+"/total_reconstructed_emisssion_"+\
                        freq_str+".pdf")
            self.plot_healpix_gnomview(total_emission_array[:,freq_index],full_resolution_pixel_index_array,\
                                    self.nr_of_super_pixels*self.super_pixel_size,title=title,min=0,max=15,rot=self.rot,pixels=self.xsize,unit='MJySr-1')
            filename="planck_zoom_"+str(self.xsize)+"_"+freq_str+"_smooth_reconstructed_total_emission_"+str(int(self.full_maps_nside))+".jpg"
            plt.savefig(self.optimizer_plots_folder+"/"+filename)
            plt.savefig(self.optimizer_plots_folder+"/"+freq_str+"/"+filename)



    def plot_total_difference_emission(self,data_dict):
        total_difference_array =data_dict["total_difference_array"]
        full_resolution_pixel_index_array = data_dict["full_resolution_pixel_index_array"]
        for freq_index in range(self.nfreq):
            freq_str = str(int(self.freq_array[freq_index]))
            title = "Total Difference Emission at "+ freq_str+" GHz"
            print("Title is",title)
            self.plot_healpix_mollview(total_difference_array[:,freq_index],full_resolution_pixel_index_array,\
                                    self.nr_of_super_pixels*self.super_pixel_size,title=title,min=0,max=15,unit='MJySr-1')
            full_sky_name  = "total_difference_emisssion_"+freq_str+".jpg"
            plt.savefig(self.optimizer_plots_folder+"/"+full_sky_name)
            plt.savefig(self.optimizer_plots_folder+"/"+freq_str+"/"+full_sky_name)
            plt.savefig(self.optimizer_plots_folder+"/"+full_sky_name)
            self.plot_healpix_gnomview(total_difference_array[:,freq_index],full_resolution_pixel_index_array,\
                                    self.nr_of_super_pixels*self.super_pixel_size,title=title,min=0,max=15,rot=self.rot,pixels=self.xsize,unit='MJySr-1')
            filename="planck_zoom_"+str(self.xsize)+"_"+freq_str+"_smooth_total_difference_emission_"+str(int(self.full_maps_nside))+".jpg"
            plt.savefig(self.optimizer_plots_folder+"/"+filename)
            plt.savefig(self.optimizer_plots_folder+"/"+freq_str+"/"+filename)


if __name__ == "__main__":
    start_time = time.time()

    #p = SkyAnalysis("powell_nside_32_priors_T_unfixed")
    p = SkyPlots("test_bayestar2019",run_type='optimizer')
    p = SkyPlots("bayestar_2019_full_sky_beta_fixed",run_type='optimizer')
    #p = SkyPlots("tiny_cepheus_beta_varying",run_type='optimizer')
    p.set_up_analysis()
    p.load_data()
    p.run_optimizer()
    data_dict = p.load_optimizer_sky_data()
    p.plot_smooth_dEBV()
    p.plot_final_optimized_functions(data_dict)
    p.plot_optimizer_sky_parameters(data_dict)
    p.plot_reconstructed_total_emission(data_dict)
    p.plot_total_difference_emission(data_dict)


    time_string = utils.end_time(start_time)


