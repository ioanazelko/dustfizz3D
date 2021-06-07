import os
import sys

from configparser import ConfigParser                  ### For parsing the configuration file
import healpy as hp
import h5py
import matplotlib.pyplot as plt
from multiprocessing import Process, Manager  ### To run the processes in parallel on multiple corres
import numpy as np
import str2bool                      ### For converting the read value from the configuration parser from string to bool
import time


DUST_3D_TEMPERATURE_MAP_DATA_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_DATA_LOCATION"]
DUST_3D_TEMPERATURE_MAP_CODE_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_CODE_LOCATION"]
DUST_3D_TEMPERATURE_MAP_PAPER_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_PAPER_LOCATION"]
DUST_3D_TEMPERATURE_MAP_PLOTS_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_PLOTS_LOCATION"]
sys.path.insert(0, DUST_3D_TEMPERATURE_MAP_CODE_LOCATION)


import data_processing
import model
import optimizer
import utils
import sampler



class SkyAnalysis():
    def __init__(self, run_name,run_type):
        self.run_type= run_type
        if run_type == "sampler":
            self.sampler_run_name=run_name
            sampler_parser = ConfigParser()
            sampler_parser.read(DUST_3D_TEMPERATURE_MAP_CODE_LOCATION+"/configurations/sampler/"+self.sampler_run_name+".cfg")
            self.sampler_parser = sampler_parser
            self.optimizer_run_name=self.sampler_parser.get('Sampler_configuration','optimizer_run_name')
        elif run_type == "optimizer":
            self.optimizer_run_name=run_name

        optimizer_parser = ConfigParser()
        optimizer_parser.read(DUST_3D_TEMPERATURE_MAP_CODE_LOCATION+"/configurations/optimizer/"+self.optimizer_run_name+".cfg")
        
        self.optimizer_parser = optimizer_parser
        self.set_analysis_parameters()
        #self.load_data()
        
    ### Model and Analysis Set Up

    def set_analysis_parameters(self):
        """ Define the parameters of the analysis"""
        ######### General sky parameters
        self.bayestar_version = self.optimizer_parser.get('Analysis_configuration','bayestar_version')
        self.full_maps_nside=int(self.optimizer_parser.get('Analysis_configuration','full_maps_nside'))
        
        self.super_pixel_nside=int(self.optimizer_parser.get('Analysis_configuration','super_pixel_nside'))
        self.nr_of_super_pixels=self.super_pixel_nside**2*12
        self.super_pixel_size = int((self.full_maps_nside/self.super_pixel_nside)**2)
        self.sky_area = self.optimizer_parser.get('Analysis_configuration','sky_area')
        self.zoom_in_area = self.optimizer_parser.get('Analysis_configuration','zoom_in_area')
        if self.sky_area == "rho_ophiuchi" or "cepheus":
            sky_area_dict=utils.get_sky_area_parameters(self.sky_area, self.super_pixel_nside)
            self.start_super_pixel=sky_area_dict['start_super_pixel']
            self.end_super_pixel=sky_area_dict['end_super_pixel']
        elif self.sky_area == "full_sky":
            self.start_super_pixel=0
            self.end_super_pixel=self.nr_of_super_pixels

        else:
            raise ValueError("Please specify the right region of the sky")
        ########## Zooming in parameters
        # this is not necessarily the same as sky area, because if sky area is the full sky, we may still want to have a way to zoom in
        zoom_in_dict=utils.get_sky_area_zoom_in_parameters(self.zoom_in_area)
        self.rot = zoom_in_dict['rot']
        self.xsize= zoom_in_dict['xsize']
        ######### EBV map parameters
        self.use_smooth_EBV = str2bool.str2bool(self.optimizer_parser.get('Analysis_configuration','use_smooth_EBV'))
        self.first_distance_slice=int(self.optimizer_parser.get('Analysis_configuration','first_distance_slice'))### Not doing the first couple of distance slices because there might no be enough dust in them; but I should still try for curiousity to see how much is in them
        self.last_distance_slice=int(self.optimizer_parser.get('Analysis_configuration','last_distance_slice')) ## not sure why Albert was not using 31
        self.distance_slice_step=int(self.optimizer_parser.get('Analysis_configuration','distance_slice_step')) ### taking every other second distance slice to save time for now
        self.model_nslices=len(range(self.first_distance_slice, self.last_distance_slice,self.distance_slice_step))


        ######### Emission map parameters
        self.use_smooth_planck = str2bool.str2bool(self.optimizer_parser.get('Analysis_configuration','use_smooth_planck'))
        self.freq_array =np.array([217.,353.,545.,857.,2998.]) #GHz
        self.nfreq = len(self.freq_array)
        self.get_the_sigma_for_each_frequency_band()
        ######### Data structure
        ## Make the folders
        self.optimizer_plots_folder = DUST_3D_TEMPERATURE_MAP_PLOTS_LOCATION+"/optimizer/"+self.optimizer_run_name
        if os.path.isdir(self.optimizer_plots_folder):
            pass
        else:
            os.mkdir(self.optimizer_plots_folder)
            for i in range(self.nfreq):
                os.mkdir(self.optimizer_plots_folder+"/"+str(int(self.freq_array[i])))
        self.optimizer_data_folder = DUST_3D_TEMPERATURE_MAP_DATA_LOCATION+"/optimizer_fits/"+self.optimizer_run_name
        
        if os.path.isdir(self.optimizer_data_folder):
            pass
        else:
            os.mkdir(self.optimizer_data_folder)

        if self.run_type == "sampler":
            self.sampler_plots_folder = DUST_3D_TEMPERATURE_MAP_PLOTS_LOCATION+"/sampler/"+self.sampler_run_name
            if os.path.isdir(self.sampler_plots_folder):
                pass
            else:
                os.mkdir(self.sampler_plots_folder)

            self.sampler_data_folder =  DUST_3D_TEMPERATURE_MAP_DATA_LOCATION+"/sampler/"+self.sampler_run_name
            if os.path.isdir(self.sampler_data_folder):
                pass
            else:
                os.mkdir(self.sampler_data_folder)   

        ######### Make dictionary
        self.analysis_configuration_dictionary = {'full_maps_nside':self.full_maps_nside,'super_pixel_nside':self.super_pixel_nside,
                                               'nr_of_super_pixels':self.nr_of_super_pixels,'super_pixel_size':self.super_pixel_size,
                                               'first_distance_slice':self.first_distance_slice, 'last_distance_slice':self.last_distance_slice,
                                               'distance_slice_step':self.distance_slice_step, 'model_nslices':self.model_nslices,
                                               'freq_array': self.freq_array, 'nfreq':self.nfreq,
                                               'sigma_emission':self.sigma_emission, 'optimizer_run_name':self.optimizer_run_name
                                               }


    def sampler_options(self):
        self.nwalkers = int(self.sampler_parser.get('Sampler_configuration','nwalkers'))
        self.ntemps = int(self.sampler_parser.get('Sampler_configuration','ntemps'))
        self.nruns= int(self.sampler_parser.get('Sampler_configuration','nruns'))
        self.nthreads = int(self.sampler_parser.get('Sampler_configuration','nthreads'))
        self.thinning = int(self.sampler_parser.get('Sampler_configuration','thinning'))
        self.use_priors_in_sampler = str2bool.str2bool(self.sampler_parser.get('Sampler_configuration','use_priors_in_sampler'))
        self.sampler_configuration_parameters = {"nwalkers":self.nwalkers,"ntemps":self.ntemps,"nruns":self.nruns,
                                                "nthreads":self.nthreads, "thinning":self.thinning,
                                                'use_priors_in_sampler':self.use_priors_in_sampler}

    def model_options(self):
        self.fit_for_the_offsets = str2bool.str2bool(self.optimizer_parser.get('Model_configuration','fit_for_the_offsets'))
        self.fit_for_the_rhos = str2bool.str2bool(self.optimizer_parser.get('Model_configuration','fit_for_the_rhos'))
        self.fit_for_the_betas = str2bool.str2bool(self.optimizer_parser.get('Model_configuration','fit_for_the_betas'))
        self.fit_for_the_Ts = str2bool.str2bool(self.optimizer_parser.get('Model_configuration','fit_for_the_Ts'))
        self.fixed_rho_along_sightline = str2bool.str2bool(self.optimizer_parser.get('Model_configuration','fixed_rho_along_sightline'))
        self.fixed_beta_along_sightline = str2bool.str2bool(self.optimizer_parser.get('Model_configuration','fixed_beta_along_sightline'))
        self.fixed_T_along_sightline = str2bool.str2bool(self.optimizer_parser.get('Model_configuration','fixed_T_along_sightline'))
        self.fixed_offset_array = np.zeros(self.nfreq)
        self.model_configuration_dictionary = {"super_pixel_size":self.super_pixel_size,
                                               "model_nslices":self.model_nslices,
                                               "fit_for_the_offsets":self.fit_for_the_offsets,
                                               "fit_for_the_rhos":self.fit_for_the_rhos,
                                               "fit_for_the_betas":self.fit_for_the_betas,
                                               "fit_for_the_Ts":self.fit_for_the_Ts,
                                               "fixed_rho_along_sightline":self.fixed_rho_along_sightline,
                                               "fixed_beta_along_sightline":self.fixed_beta_along_sightline,
                                               "fixed_T_along_sightline":self.fixed_T_along_sightline,
                                               "fixed_offset_array":self.fixed_offset_array}

    def optimizer_options(self):
        #optimizer_method="SLSQP" #might have some fortran errors like in here: https://github.com/scipy/scipy/issues/10740
        #optimizer_method="Powell" #no gradient, at least gives an answer
        #optimizer_method="dogleg" #hessian is required
        #optimizer_method="CG" #weird error about being below precision; it cannot handle bounds
        #optimizer_method="BFGS" #can't handle bounds
        #optimizer_method="Newton-CG" #can't handle bounds
        #optimizer_method="TNC" #at least gave a result; much faster than Powell,but much less precise
        #optimizer_method="trust-ncg"
        #optimizer_method="trust-exact"
        #optimizer_method="trust-constr"
        optimizer_method=self.optimizer_parser.get('Optimizer_configuration','optimizer_method')
        if optimizer_method in ["CG","BFGS","Newton-CG"]:
            use_bounds = False
        else:
            use_bounds = True
        scale_variables = str2bool.str2bool(self.optimizer_parser.get('Optimizer_configuration','scale_variables'))
        use_gradient = str2bool.str2bool(self.optimizer_parser.get('Optimizer_configuration','use_gradient'))
        use_hessian = str2bool.str2bool(self.optimizer_parser.get('Optimizer_configuration','use_hessian'))
        use_priors = str2bool.str2bool(self.optimizer_parser.get('Optimizer_configuration','use_priors'))
        print_output = str2bool.str2bool(self.optimizer_parser.get('Optimizer_configuration','print_output'))
        self.optimizer_configuration_parameters = {"optimizer_method":optimizer_method,
                                                   "scale_variables":scale_variables,
                                                   "use_bounds":use_bounds,
                                                   "use_gradient":use_gradient,
                                                   "use_hessian":use_hessian,
                                                   "use_priors":use_priors,
                                                   "print_output":print_output}
    ### Loading data
    def load_data(self):
        data_options_dict ={'bayestar_nside':self.full_maps_nside,
                            'bayestar_version':self.bayestar_version,
                            'planck_nside' : self.full_maps_nside,
                            'planck_version':"albert"}
        data = data_processing.DataProcessing(data_options_dict)
        self.planck = data.load_planck()
        self.planck_smooth = data.load_smooth_planck()
        ### bayestar2019 or bayestar2017
        self.ebv = data.load_bayestar()
        self.ebv_smooth = data.load_smooth_ebv()
        ### distances
        self.bayestar_distances = data.load_distances()
        self.bayestar_ndistances = len(self.bayestar_distances)
        self.model_dist_slices=self.bayestar_distances[self.first_distance_slice:self.last_distance_slice:self.distance_slice_step] # get the selected distances as well
        self.analysis_configuration_dictionary['model_dist_slices']=self.model_dist_slices
        self.dEBV = self.ebv_for_each_bin(self.ebv)
        self.dEBV_smooth = self.ebv_for_each_bin(self.ebv_smooth)
        if self.use_smooth_EBV == True:
            self.analysis_ebv = self.ebv_smooth
            self.analysis_dEBV= self.dEBV_smooth
        else:
            self.analysis_ebv = self.ebv
            self.analysis_dEB = self.dEBV
        if self.use_smooth_planck ==  True:
            self.analysis_planck = self.planck_smooth
        else:
            self.analysis_planck = self.planck
        return

    def print_analysis_parameters(self):
        print("Nr of SuperPixels: ", self.nr_of_super_pixels)
        print("Nr of mini pixels in a Super Pixel: ",self.super_pixel_size)
        print("Nr of distance slices: ", self.model_nslices)
    #### Emission functions
    def get_the_sigma_for_each_frequency_band(self):
        sigma_type = self.optimizer_parser.get('Analysis_configuration','sigma_freq_type')
        if sigma_type == "albert_small":
            self.sigma_emission = [.005,.02, .05, .1, .1]
        elif sigma_type == "albert_large":
            self.sigma_emission = [.5,3.,9.,25.,30.]
        else:
            raise ValueError("Specify the emission error bars")
    def get_uncertainties_for_one_super_pixel(self):
        # Set uncertainties in Planck bands and 100um, for one SuperPixel
        super_pixel_sigma = np.zeros([self.super_pixel_size ,self.nfreq])
        super_pixel_sigma[:,:] = self.sigma_emission # For each sightline in the pixel, assign the 5 frequency errors
        self.super_pixel_sigma_array = super_pixel_sigma

    #### EBV functions   
    def ebv_for_each_bin(self,ebv): 
        ##### Calculating the EBV contribution for each distance bin
        ### get the EBV till the selected distance bins
        select_dist_EBV = ebv[self.first_distance_slice:self.last_distance_slice:self.distance_slice_step]
        ### differentiate to get get the individual contribution from each distance bin (bayestar is cumulative)
        diff_EBV = select_dist_EBV[1:,:] - select_dist_EBV[:-1,:]
        ### append the cumulative EBV till the first distance bin choice to the first value
        dEBV = np.append(ebv[self.first_distance_slice:self.first_distance_slice+1],diff_EBV,axis=0)
        return dEBV

    #ef custom_distance_slice(custom_distance_array_type):
    #    if custom_distance_array_type == "cepheus":


    def set_up_analysis(self):      
        ### Set Analysis Parameters
        #self.set_analysis_parameters()
        self.model_options()
        self.get_uncertainties_for_one_super_pixel()
        self.print_analysis_parameters()


    def save_optimizer_sky_data(self,data_dict,file_index=""):
        with h5py.File(self.optimizer_data_folder+"/"+self.optimizer_run_name+file_index+".hdf5", "w") as f:
            f.create_dataset("final_parameters_array",data=data_dict["final_parameters_array"] ,dtype='float64')
            f.create_dataset("final_optimized_functions_array",data=data_dict["final_optimized_functions_array"],dtype='float64')
            f.create_dataset("final_chi_square_array",data=data_dict["final_chi_square_array"],dtype='float64')
            f.create_dataset("super_pixels_index_array",data=data_dict["super_pixels_index_array"],dtype='int')
            f.create_dataset("voxel_emission_array",data=data_dict["voxel_emission_array"],dtype='float64')
            f.create_dataset("total_emission_array",data=data_dict["total_emission_array"],dtype='float64')
            f.create_dataset("total_difference_array",data=data_dict["total_difference_array"],dtype='float64')
            f.create_dataset("full_resolution_pixel_index_array",data=data_dict["full_resolution_pixel_index_array"],dtype='int')
            f.attrs["optimizer_configuration"] = str(self.optimizer_configuration_parameters)
            f.attrs["model_configuration"] = str(self.model_configuration_dictionary)
            f.attrs["analysis_configuration"]= str(self.analysis_configuration_dictionary)
            f.close()
    def load_optimizer_sky_data(self,file_index=""):
        data_dict = {}
        with h5py.File(self.optimizer_data_folder+"/"+self.optimizer_run_name+file_index+".hdf5", "r") as g:    
            data_dict["final_parameters_array"]=g["final_parameters_array"][()]
            data_dict["final_optimized_functions_array"] = g["final_optimized_functions_array"][()]
            data_dict["final_chi_square_array"] = g["final_chi_square_array"][()]
            data_dict["super_pixels_index_array"] = g["super_pixels_index_array"][()]
            data_dict["voxel_emission_array"] = g["voxel_emission_array"][()]
            data_dict["total_emission_array"] = g["total_emission_array"][()]
            data_dict["total_difference_array"] = g["total_difference_array"][()]
            data_dict["full_resolution_pixel_index_array"] = g["full_resolution_pixel_index_array"][()]
            stored_optimizer_configuration = g.attrs["optimizer_configuration"]
            stored_model_configuration = g.attrs["model_configuration"]
            stored_analysis_configuration = g.attrs["analysis_configuration"]
            g.close()
        return data_dict

    def separate_sky_optimizer_parameters(self,parameters):
        """
        This function parses the sky array of parameters resulting from the optimizer based on the options specified for the model
        """
        ### The order of parameters is:
        ### Offsets (0 or 5, one for each frequency), rho (could be 0,1, or nslices)
        ### beta (could be 0, 1, or nslices), and T (could be nslices)
        n_super_pixels = len(parameters)
        param_count=0 
        ### Parsing the offsets
        if self.fit_for_the_offsets == True:
            offsets = parameters[:,0:self.nfreq]
            param_count+=self.nfreq
        else:
            offsets = np.tile(self.fixed_offset_array,(n_super_pixels,1)) ### 5 nr corresponding to the 5 frequencies
            ### we don't update the index count because there is no value stored in the parameters for the offsets
        ###  Parsing the rhos
        if self.fit_for_the_rhos == True:
            ### case for rho varying in each voxel
            if self.fixed_rho_along_sightline == False:
                rhos = parameters[:,param_count:param_count+self.model_nslices]
                param_count+=self.model_nslices
            ### case for rho only varying in each superpixel
            else:
                rhos = parameters[:,param_count] # reading the value from the parameters since each superpixel now has one varying rho
                param_count+=1
        ### case for rho being fixed across the sky
        else:
            preset_rho = self.model_configuration_dictionary["preset_rho"]
            rhos = np.array([preset_rho]*n_super_pixels)
            ### we don't update the index count because there is no value stored in the parameters for rho

        ### Parsing the betas
        if self.fit_for_the_betas == True:
            ### case of beta varying in each voxel
            if self.fixed_beta_along_sightline == False:
                betas = parameters[:,param_count:param_count+self.model_nslices]
                param_count+=self.model_nslices
            ### case for beta only varying in each superpixel
            else:
                betas = parameters[:,param_count]
                param_count+=1
        ### case for beta being fixed across the sky
        else:
            preset_beta = self.model_configuration_dictionary["preset_beta"]
            betas = np.array([preset_beta]*n_super_pixels)
            ### we don't update the index count because there is no value stored in the parameters for beta
        ### Parsing the Ts
        if self.fit_for_the_Ts == True:
            ### case of T varying in each voxel
            if self.fixed_T_along_sightline == False:
                Ts = parameters[:,param_count:param_count+self.model_nslices]
                param_count+=self.model_nslices
            else:
                Ts = parameters[:,param_count]
        else:
            raise ValueError("What are you doing? You want fixed T? implement it")

        return offsets,rhos,betas,Ts 


    def calculate_reconstructed_emission_data(self, parameters_array,super_pixels_index_array):
        n_chosen_super_pix = len(super_pixels_index_array) #Number of superpixels
        total_emission_array = np.zeros((n_chosen_super_pix*self.super_pixel_size,self.nfreq)) #the emission map at each frequency (pixel_index,frequency)
        voxel_emission_array = np.zeros((n_chosen_super_pix*self.super_pixel_size,self.model_nslices,self.nfreq)) #the individual contribution of each voxel to the emission at each frequency (pixel_index, distance_slice,frequency)
        total_difference_array = np.zeros((n_chosen_super_pix*self.super_pixel_size,self.nfreq)) #the difference between reconstruction and original emission map at each frequency (pixel_index,frequency)
        full_resolution_pixel_index_array = np.zeros(n_chosen_super_pix*self.super_pixel_size,dtype=int)
        for i in range(n_chosen_super_pix):
            # selecting the right data for the SuperPixel from the entire sky array
            super_pixel_dEBV = self.analysis_dEBV[:,super_pixels_index_array[i]*self.super_pixel_size:(super_pixels_index_array[i]+1)*self.super_pixel_size].transpose()
            ### after the transpose, the dimmensions are now (pixel_index,distance_slice)
            super_pixel_emission = self.analysis_planck[:,super_pixels_index_array[i]*self.super_pixel_size:(super_pixels_index_array[i]+1)*self.super_pixel_size].transpose()
            ### after the transpose, the dimmensions are now (pixel_index,frequency)
            data_dictionary = {"super_pixel_dEBV":super_pixel_dEBV,
                               "super_pixel_emission": super_pixel_emission ,
                               "freq_array":self.freq_array,
                               "super_pixel_sigma_array":self.super_pixel_sigma_array}
            m = model.Model(model_configuration_dictionary=self.model_configuration_dictionary,
                           data_dictionary = data_dictionary)
            parameters=parameters_array[i,:]
            voxel,total = m.calculate_voxel_and_total_emission_for_reconstruction(parameters)
            voxel_emission_array[i*self.super_pixel_size:(i+1)*self.super_pixel_size] = voxel
            total_emission_array[i*self.super_pixel_size:(i+1)*self.super_pixel_size] = total
            full_resolution_pixel_index_array[i*self.super_pixel_size:(i+1)*self.super_pixel_size] = \
                np.array(range(super_pixels_index_array[i]*self.super_pixel_size,(super_pixels_index_array[i]+1)*self.super_pixel_size))
            total_difference = super_pixel_emission-total
            total_difference_array[i*self.super_pixel_size:(i+1)*self.super_pixel_size] = total_difference
        return voxel_emission_array, total_emission_array,total_difference_array, full_resolution_pixel_index_array

    


    ### Run the optimizer
    def run_optimizer(self):
        self.optimizer_options()
        # def test_run_sampler(self):
        #     #Selecting the super pixels for the fit
        print("Starting the optimizer") 

        start_pixel = self.start_super_pixel
        #end_pixel = int(self.nr_of_super_pixels/100)
        #end_pixel = self.nr_of_super_pixels 
        end_pixel= self.end_super_pixel

        super_pixels_index_array = np.array(range(start_pixel,end_pixel))
        n_chosen_super_pix = len(super_pixels_index_array) #Number of superpixels
        
        nr_of_parallel_processes = 30
        if n_chosen_super_pix%nr_of_parallel_processes !=0:
            print("number of super pixels chosen for fit ", n_chosen_super_pix)
            print("number of parallel process ", nr_of_parallel_processes)
            raise ValueError("Wrong nr of parralel processes or super pixels!!!")
        part_n_super_pixels = int(n_chosen_super_pix/nr_of_parallel_processes)


        start_time = time.time()

        processes = []

        for process_index in range(nr_of_parallel_processes):
            def do_it():
                sys.stdout.flush()
                ### Selecting the superpixels that each parallel process should take care of
                part_super_pixels_index_array = super_pixels_index_array[process_index*part_n_super_pixels:(process_index+1)*part_n_super_pixels]
                parameters_list = []
                optimized_function_list = []
                final_chi_square_list = []
                for i in range(part_n_super_pixels):
                    print("Running the optimizer for super_pixel ",i, 'process index', process_index)
                    ### To implement: disregard subpixels that are masked out
                    ### To implement: only fit superpixels with enough coverage by the extinction map

                    # transpose arrays here so that matrix operations work during model evaluation
                    # selecting the right data for the SuperPixel from the entire sky array
                    super_pixel_dEBV = self.analysis_dEBV[:,part_super_pixels_index_array[i]*self.super_pixel_size:(part_super_pixels_index_array[i]+1)*self.super_pixel_size].transpose()
                    ### after the transpose, the dimmensions are now (pixel_index,distance_slice)
                    super_pixel_emission = self.analysis_planck[:,part_super_pixels_index_array[i]*self.super_pixel_size:(part_super_pixels_index_array[i]+1)*self.super_pixel_size].transpose()
                    ### after the transpose, the dimmensions are now (pixel_index,frequency)

                    #print(super_pixel_dEBV)
                    data_dictionary = {"super_pixel_dEBV":super_pixel_dEBV,
                                       "super_pixel_emission": super_pixel_emission ,
                                       "freq_array":self.freq_array,
                                       "super_pixel_sigma_array":self.super_pixel_sigma_array}

                    o = optimizer.Optimizer(optimizer_configuration_parameters = self.optimizer_configuration_parameters,
                                            model_configuration_dictionary = self.model_configuration_dictionary,
                                            data_dictionary=data_dictionary)

                    parameters, optimized_function_result, final_chi_square= o.run_optimizer()
                    parameters_list.append(parameters)
                    optimized_function_list.append(optimized_function_result)
                    final_chi_square_list.append(final_chi_square)

                final_parameters_array = np.array(parameters_list)
                final_optimized_function_array = np.array(optimized_function_list)
                final_chi_square_array = np.array(final_chi_square_list)

                print("I am calculating the final emission accross the sky and voxels")
                voxel_emission_array, total_emission_array, total_difference_array, full_resolution_pixel_index_array=\
                    self.calculate_reconstructed_emission_data(final_parameters_array,part_super_pixels_index_array)
                partial_data_dict = {}
                partial_data_dict["final_parameters_array"]=final_parameters_array
                partial_data_dict["final_optimized_functions_array"] = final_optimized_function_array
                partial_data_dict["final_chi_square_array"] = final_chi_square_array
                partial_data_dict["super_pixels_index_array"] = part_super_pixels_index_array
                partial_data_dict["voxel_emission_array"] = voxel_emission_array
                partial_data_dict["total_emission_array"] = total_emission_array
                partial_data_dict["total_difference_array"] = total_difference_array
                partial_data_dict["full_resolution_pixel_index_array"] = full_resolution_pixel_index_array

                self.save_optimizer_sky_data(partial_data_dict,file_index="_"+str(process_index))
                #print(parallel_dict.keys())
                    

            p = Process(target = do_it)
            p.start()
            processes.append(p)
        for job in processes:# this makes the program wait for all the processes to be done
            job.join()
        print("I am done running the optimizers for superpixel ",start_pixel," to ", end_pixel)
        time_string = utils.end_time(start_time)
        # print("Saving the data")  
        parallel_dict={}

        for process_index in range(nr_of_parallel_processes):

            parallel_dict[str(process_index)]=self.load_optimizer_sky_data(file_index="_"+str(process_index))

        ### Now I need to concatenate the data from the processes into an individual coherent chain
        data_dict={}
        for key in parallel_dict['0'].keys():
            ## start the numpy arrays by reading in the first slice
            data_array = parallel_dict['0'][key]
            if nr_of_parallel_processes>1:
                for process_index in range(1, nr_of_parallel_processes):
                    data_array = np.concatenate((data_array,parallel_dict[str(process_index)][key]))
            data_dict[key]=data_array
        #### Saving the data
        self.save_optimizer_sky_data(data_dict)

    ############################################
    #### Sampler functions
    ############################################

    #### Run the sampler
    
    
    def run_sampler(self):
        #Selecting the super pixels for the fit
        self.sampler_options()
        print("running the sampler now")
        optimizer_data_dict = self.load_optimizer_sky_data()
        optimizer_parameters = optimizer_data_dict["final_parameters_array"]
        optimizer_chi_square = optimizer_data_dict["final_chi_square_array"]
        start_pixel  =self.start_super_pixel
        end_pixel = self.end_super_pixel 
        
        super_pixels_index_array = np.array(range(start_pixel,end_pixel))
        n_chosen_super_pix = len(super_pixels_index_array) #Number of superpixels

        nr_of_parallel_processes = 1
        if n_chosen_super_pix%nr_of_parallel_processes !=0:
            raise ValueError("Wrong nr of parralel processes or super pixels!!!")
        part_n_super_pixels = int(n_chosen_super_pix/nr_of_parallel_processes)
        
        

        start_time = time.time()

        processes = []
        for process_index in range(nr_of_parallel_processes):
            def do_it():
                sys.stdout.flush()
                part_super_pixels_index_array = super_pixels_index_array[process_index*part_n_super_pixels:(process_index+1)*part_n_super_pixels]
                part_optimizer_parameters = optimizer_parameters[process_index*part_n_super_pixels:(process_index+1)*part_n_super_pixels]
                part_optimizer_chi_square = optimizer_chi_square[process_index*part_n_super_pixels:(process_index+1)*part_n_super_pixels]

                parameters_list = []
                optimized_function_list = []
                final_chi_square_list = []
                #for i in range(part_n_super_pixels):
                for i in range(201,202):
                    super_pixel_index=part_super_pixels_index_array[i]
                    #print("Running the sampler for super_pixel ",super_pixel_index)
                    ### To implement: disregard subpixels that are masked out
                    ### To implement: only fit superpixels with enough coverage by the extinction map
                    print("original optimizer chi square at super pixel index",super_pixel_index, part_optimizer_chi_square[i])

                    # transpose arrays here so that matrix operations work during model evaluation
                    # selecting the right data for the SuperPixel from the entire sky array
                    super_pixel_dEBV = self.analysis_dEBV[:,super_pixel_index*self.super_pixel_size:(super_pixel_index+1)*self.super_pixel_size].transpose()
                    ### after the transpose, the dimmensions are now (pixel_index,distance_slice)
                    super_pixel_emission = self.analysis_planck[:,super_pixel_index*self.super_pixel_size:(super_pixel_index+1)*self.super_pixel_size].transpose()
                    ### after the transpose, the dimmensions are now (pixel_index,frequency)
                    initial_optimizer_parameters = part_optimizer_parameters[i]
                    ### you could also select based on whether the optimizer conditions were ok
                    
                    #print(super_pixel_dEBV)
                    data_dictionary = {"super_pixel_dEBV":super_pixel_dEBV,
                                       "super_pixel_emission": super_pixel_emission ,
                                       "freq_array":self.freq_array,
                                       "super_pixel_sigma_array":self.super_pixel_sigma_array}

                    sampler_super_pixel_specific_configuration ={"initial_optimizer_parameters":initial_optimizer_parameters,
                                                                "super_pixel_index":super_pixel_index,
                                                                "sampler_plots_folder":self.sampler_plots_folder+"/"+str(super_pixel_index),
                                                                "sampler_data_folder": self.sampler_data_folder+"/"+str(super_pixel_index)}
                    

                    s = sampler.Map3DSampler(sampler_configuration_parameters = self.sampler_configuration_parameters,
                        model_configuration_dictionary = self.model_configuration_dictionary,
                        data_dictionary=data_dictionary,sampler_super_pixel_specific_configuration=sampler_super_pixel_specific_configuration)
    #                 print(s.calculate_chi_square_with_priors(initial_optimizer_parameters))
    #                 print("original optimizer ln_prior: ",s.ln_prior(initial_optimizer_parameters))
    #               s.get_initial_positions(initial_optimizer_parameters)

                  
                    #check_initial_positions(s)
                    s.sample()
                    #process_sampler_data(s)

                    data_sampler=s.load_run_data()
                    s.plot_chain_corner_plots(data_sampler,throw_away_burnout=True)
                    s.plot_chain_corner_plots(data_sampler,throw_away_burnout=False)
                    

            p = Process(target = do_it)
            p.start()
            processes.append(p)
        for job in processes:# this makes the program wait for all the processes to be done
            job.join()
        print("I am done running the optimizers for superpixel ",start_pixel," to ", end_pixel)
        time_string = utils.end_time(start_time)
        # print("Saving the data")  
        parallel_dict={}



if __name__ == "__main__":
    start_time = time.time()
    #plot_the_paper_plots()
    #make_the_analysis()
    p = SkyAnalysis("cepheus_beta_varying",run_type="optimizer")        
    #p = SkyAnalysis("test_bayestar2019",run_type="optimizer")        
    #p = SkyAnalysis("powell_nside_128_priors_T_unfixed")        
    #p = SkyAnalysis("powell_nside_64_priors_T_unfixed")    
    #p = SkyAnalysis("powell_nside_32_priors_T_unfixed")
    #p = SkyAnalysis("powell_nside_32_no_priors")
    p.set_up_analysis()
    p.load_data()
    p.run_optimizer()

    #making_the_tables_of_the_paper()
    time_string = utils.end_time(start_time)



   
