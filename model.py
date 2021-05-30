### Created 2020.11.08 By Ioana Zelko
### Makes the model for the sampler to use

from __future__ import division,print_function

import os
import sys


import numpy as np
DUST_3D_TEMPERATURE_MAP_DATA_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_DATA_LOCATION"]
DUST_3D_TEMPERATURE_MAP_CODE_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_CODE_LOCATION"]
DUST_3D_TEMPERATURE_MAP_PAPER_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_PAPER_LOCATION"]
DUST_3D_TEMPERATURE_MAP_PLOTS_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_PLOTS_LOCATION"]
sys.path.insert(0, DUST_3D_TEMPERATURE_MAP_CODE_LOCATION)

from configparser import ConfigParser                  ### For parsing the configuration file
UNIVERSAL_CONSTANTS_LOCATION = os.environ["UNIVERSAL_CONSTANTS_LOCATION"]

parser = ConfigParser()
parser.read(UNIVERSAL_CONSTANTS_LOCATION+"/constant_configuration.cfg")
c = np.float64(parser.get('Universal Constants','c')) ## m*s-1
h = np.float64(parser.get('Universal Constants','h')) ## J*s
k = np.float64(parser.get('Universal Constants','k')) ## J*K-1
#T_cmb = np.float64(parser.get('Universal Constants','T_cmb')) ## K

import utils

class Model():
    def __init__(self, model_configuration_dictionary,data_dictionary):

        ### Loading the Data
        self.super_pixel_dEBV=data_dictionary["super_pixel_dEBV"] #(pixel_index,distance_slice)
        self.super_pixel_emission=data_dictionary["super_pixel_emission"] #(pixel_index,frequency)
        self.super_pixel_sigma_array=data_dictionary["super_pixel_sigma_array"] #(pixel_index,frequency)
        self.freq_array = data_dictionary["freq_array"] #[217,353,545,857,2998] #GHz
        self.nfreq = len(self.freq_array)
        ### Model Configuration
        self.model_configuration_dictionary = model_configuration_dictionary
        self.super_pixel_size = model_configuration_dictionary["super_pixel_size"]
        self.model_nslices =  model_configuration_dictionary["model_nslices"]
        self.fit_for_the_offsets = model_configuration_dictionary["fit_for_the_offsets"]
        self.fit_for_the_rhos = model_configuration_dictionary["fit_for_the_rhos"]
        self.fit_for_the_betas = model_configuration_dictionary["fit_for_the_betas"]
        self.fit_for_the_Ts = model_configuration_dictionary["fit_for_the_Ts"]
        ####
        self.T0 = 18.8    #K
        self.sigma_T = 9. #K

        ### Make the model attributes
        self.make_model_attributes()

    def array_ISM_dust_MBB_with_deriv(self,rhos, betas,Ts):
        """
        nu_array is expected to be in GHz. Make sure it is floats! Otherwise bad stuff will happen when you take the power of 3.
        rho : array of rho for each distance slice
        betas: array of beta for each distance slice
        Ts: array of T for each distance slice
        """
        nu0 = 353.*1E9 # reference nu, Hz
        scaled_nu_array = self.freq_array*10**9 #to transform from GHz to Hz
        x_D = np.outer(1/Ts,h*scaled_nu_array/k) #this creates an array (distance, freq)
        B = 1/(np.exp(x_D)-1)*2*h*scaled_nu_array**3/c**2*1E26 *1E-6# 1E26 is to transform from J/m^2/sr to Jy/sr,1E-6 to move to mega 
        ones = np.ones(self.model_nslices) #used for the outer product to propagate nu_array for all distance slices
        I_over_ebv = (rhos*((np.outer(scaled_nu_array/nu0,ones))**betas)).T*B #(distance,freq) #the transpose moves (freq,distance) to (distance,freq)
        
        # This would be correct, but it does some more operations than necessary, so for a faster code I am using the other less clear expression
        #derivBdT = x_D*np.outer(1/Ts,np.ones(self.nfreq))*np.exp(x_D)/(np.exp(x_D)-1)**2*2*h*scaled_nu_array**3/c**2*1E26 *1E-6
        #dI_ver_ebv_dT =  (rhos*((np.outer(scaled_nu_array/nu0,ones))**betas)).T*derivBdT#(distance,freq) #the transpose moves (freq,distance) to (distance,freq)
        dI_over_ebv_drho = I_over_ebv/np.outer(rhos,np.ones(self.nfreq))
        dI_over_ebv_dbeta = I_over_ebv *np.log(scaled_nu_array/nu0)
        dI_over_ebv_dT = I_over_ebv *x_D*np.outer(1/Ts,np.ones(self.nfreq))/(np.exp(x_D)-1)
        return I_over_ebv, dI_over_ebv_drho, dI_over_ebv_dbeta, dI_over_ebv_dT 

    def calculate_chi_square_and_its_gradient(self,parameters):
        offsets,rhos,betas,Ts = self.separate_parameters(parameters)
        I_over_ebv, dI_over_ebv_drho, dI_over_ebv_dbeta, dI_over_ebv_dT =self.array_ISM_dust_MBB_with_deriv(rhos, betas,Ts) #(distance, frequency)
        
        ### Calculating the total intensity
        total_intensity_no_offsets = np.dot(self.super_pixel_dEBV,I_over_ebv) #matrix multiplication of  (pixel,distance)X(distance,frequency) = (pixel, frequency)
        total_intensity = total_intensity_no_offsets + offsets #(pixel, frequency)
        chi_square=np.sum(((total_intensity-self.super_pixel_emission )/self.super_pixel_sigma_array)**2) # summing over the frequencies and the pixels in the superpixel

        ### Calculating the gradients
        g_k_nu = 2*(total_intensity-self.super_pixel_emission)/self.super_pixel_sigma_array**2 #(pixel,frequency)

        gradient = np.array([])
        if self.fit_for_the_offsets == True:
            offset_gradient = np.sum(g_k_nu, axis=0)  # nshape = (nfreq)
            gradient = np.concatenate((gradient,offset_gradient))
        #### rho
        if self.fit_for_the_rhos == True:
            fixed_rho_along_sightline = self.model_configuration_dictionary["fixed_rho_along_sightline"]
            ### case for rho varying in each voxel
            if fixed_rho_along_sightline == False:
                #rho_gradient[j] =      g_k_nu    *    self.super_pixel_dEBV[:,j] *   dI_over_ebv_drho[j,:]
                #sum_{pixel, frequency) (pixel, frequency) * (pixel,distance )    *  (distance, frequency))
                #np.sum((dI_over_ebv_drho[j,:] *g_k_nu).T *self.super_pixel_dEBV[:,j],axis=0 )
                rho_gradient =  np.sum(np.dot(self.super_pixel_dEBV.T, g_k_nu )*dI_over_ebv_drho,axis=1) #nshape = (model_nslices)
                gradient = np.concatenate((gradient,rho_gradient))
            else:
                rho_gradient =  np.sum(g_k_nu*np.dot(self.super_pixel_dEBV,dI_over_ebv_drho))# 1 number
                gradient = np.append(gradient,rho_gradient)
        #print("After rho the gradient len is", len(gradient))
       ### betas
        if self.fit_for_the_betas == True:
            fixed_beta_along_sightline=self.model_configuration_dictionary["fixed_beta_along_sightline"]
            ### case for beta varying in each voxel
            if fixed_beta_along_sightline == False:
                beta_gradient = np.sum(np.dot(self.super_pixel_dEBV.T,g_k_nu)*dI_over_ebv_dbeta,axis=1)#nshape = (model_nslices)
                gradient = np.concatenate((gradient,beta_gradient)) 
            else:
                beta_gradient = np.sum(g_k_nu*np.dot(self.super_pixel_dEBV,dI_over_ebv_dbeta)) # 1 number
                gradient = np.append(gradient,beta_gradient) 
        #print("After beta the gradient len is", len(gradient))

        ### Ts
        if self.fit_for_the_Ts == True:
            ### case for T varying in each voxel
            fixed_T_along_sightline=self.model_configuration_dictionary["fixed_T_along_sightline"]

            if fixed_T_along_sightline == False:
                T_gradient = np.sum(np.dot(self.super_pixel_dEBV.T,g_k_nu)*dI_over_ebv_dT,axis=1) #nshape = (model_nslices)
                gradient = np.concatenate((gradient,T_gradient))
            else:
                T_gradient = np.sum(g_k_nu*np.dot(self.super_pixel_dEBV,dI_over_ebv_dT)) # 1 number
                gradient = np.append(gradient,T_gradient) 


        return chi_square,gradient


    def calculate_chi_square_with_priors_and_its_gradient(self,parameters):
        offsets,rhos,betas,Ts = self.separate_parameters(parameters)
        I_over_ebv, dI_over_ebv_drho, dI_over_ebv_dbeta, dI_over_ebv_dT =self.array_ISM_dust_MBB_with_deriv(rhos, betas,Ts) #(distance, frequency)
        
        ### Calculating the total intensity
        total_intensity_no_offsets = np.dot(self.super_pixel_dEBV,I_over_ebv) #matrix multiplication of  (pixel,distance)X(distance,frequency) = (pixel, frequency)
        total_intensity = total_intensity_no_offsets + offsets #(pixel, frequency)
        
        fixed_T_along_sightline=self.model_configuration_dictionary["fixed_T_along_sightline"]
        if fixed_T_along_sightline == False:
            prior = np.sum(((Ts-self.T0)/self.sigma_T)**2)
        else:
            prior = ((Ts[0]-self.T0)/self.sigma_T)**2
        chi_square=np.sum(((total_intensity-self.super_pixel_emission )/self.super_pixel_sigma_array)**2) # summing over the frequencies and the pixels in the superpixel
        chi_square_prior =chi_square+prior
        ### Calculating the gradients
        g_k_nu = 2*(total_intensity-self.super_pixel_emission)/self.super_pixel_sigma_array**2 #(pixel,frequency)

        gradient = np.array([])
        if self.fit_for_the_offsets == True:
            offset_gradient = np.sum(g_k_nu, axis=0)  # nshape = (nfreq)
            gradient = np.concatenate((gradient,offset_gradient))
        #### rho
        if self.fit_for_the_rhos == True:
            fixed_rho_along_sightline = self.model_configuration_dictionary["fixed_rho_along_sightline"]
            ### case for rho varying in each voxel
            if fixed_rho_along_sightline == False:
                #rho_gradient[j] =      g_k_nu    *    self.super_pixel_dEBV[:,j] *   dI_over_ebv_drho[j,:]
                #sum_{pixel, frequency) (pixel, frequency) * (pixel,distance )    *  (distance, frequency))
                #np.sum((dI_over_ebv_drho[j,:] *g_k_nu).T *self.super_pixel_dEBV[:,j],axis=0 )
                rho_gradient =  np.sum(np.dot(self.super_pixel_dEBV.T, g_k_nu )*dI_over_ebv_drho,axis=1) #nshape = (model_nslices)
                gradient = np.concatenate((gradient,rho_gradient))
            else:
                rho_gradient =  np.sum(g_k_nu*np.dot(self.super_pixel_dEBV,dI_over_ebv_drho))# 1 number
                gradient = np.append(gradient,rho_gradient)
        #print("After rho the gradient len is", len(gradient))
       ### betas
        if self.fit_for_the_betas == True:
            fixed_beta_along_sightline=self.model_configuration_dictionary["fixed_beta_along_sightline"]
            ### case for beta varying in each voxel
            if fixed_beta_along_sightline == False:
                beta_gradient = np.sum(np.dot(self.super_pixel_dEBV.T,g_k_nu)*dI_over_ebv_dbeta,axis=1)#nshape = (model_nslices)
                gradient = np.concatenate((gradient,beta_gradient)) 
            else:
                beta_gradient = np.sum(g_k_nu*np.dot(self.super_pixel_dEBV,dI_over_ebv_dbeta)) # 1 number
                gradient = np.append(gradient,beta_gradient) 
        #print("After beta the gradient len is", len(gradient))

        ### Ts
        if self.fit_for_the_Ts == True:
            fixed_T_along_sightline=self.model_configuration_dictionary["fixed_T_along_sightline"]
            ### case for T varying in each voxel

            if fixed_T_along_sightline == False:
                T_gradient = np.sum(np.dot(self.super_pixel_dEBV.T,g_k_nu)*dI_over_ebv_dT,axis=1) #nshape = (model_nslices)
                T_prior = 2*(Ts-self.T0)/(self.sigma_T**2)
                T_gradient += T_prior
                gradient = np.concatenate((gradient,T_gradient))
            else:
                T_gradient = np.sum(g_k_nu*np.dot(self.super_pixel_dEBV,dI_over_ebv_dT)) # 1 number
                T_prior = 2*(Ts[0]-self.T0)/(self.sigma_T**2)
                T_gradient += T_prior
                gradient = np.append(gradient,T_gradient) 

        return chi_square_prior,gradient

    def calculate_chi_square_with_priors(self,parameters):
        offsets,rhos,betas,Ts = self.separate_parameters(parameters)   
        fixed_T_along_sightline=self.model_configuration_dictionary["fixed_T_along_sightline"]

        if fixed_T_along_sightline == False:
            prior = np.sum(((Ts-self.T0)/self.sigma_T)**2)
        else:
            prior = ((Ts[0]-self.T0)/self.sigma_T)**2
        prior = np.sum(((Ts-self.T0)/self.sigma_T)**2) # (1 nr)
        chi_square=self.calculate_chi_square(parameters) 
        chi_square_prior =chi_square+prior
        return chi_square_prior

    def calculate_chi_square(self,parameters):
        """This does not have priors"""
        
        total_intensity = self.calculating_total_radiation_from_slices(parameters)
        chi_square=np.sum(((self.super_pixel_emission - total_intensity)/self.super_pixel_sigma_array)**2) # summing over the frequencies and the pixels in the superpixel
        return chi_square

    def array_ISM_dust_MBB(self,rhos, betas,Ts):
        """
        nu_array is expected to be in GHz. Make sure it is floats! Otherwise bad stuff will happen when you take the power of 3.
        rho : array of rho for each distance slice
        betas: array of beta for each distance slice
        Ts: array of T for each distance slice
        """
        nu0 = 353.*1E9 # reference nu, Hz
        scaled_nu_array = self.freq_array*10**9 #to transform from GHz to Hz
        x_D = np.outer(1/Ts,h*scaled_nu_array/k) #this creates an array (distance, freq)
        B = 1/(np.exp(x_D)-1)*2*h*scaled_nu_array**3/c**2*1E26 *1E-6# 1E26 is to transform from J/m^2/sr to Jy/sr,1E-6 to move to mega 
        ones = np.ones(self.model_nslices) #used for the outer product to propagate nu_array for all distance slices
        I_over_ebv = (rhos*((np.outer(scaled_nu_array/nu0,ones))**betas)).T*B #(distance,freq) #the transpose moves (freq,distance) to (distance,freq)
        return I_over_ebv

    def calculating_total_radiation_from_slices(self,parameters):
        offsets,rhos,betas,Ts = self.separate_parameters(parameters) 
        d_emission_over_dEBV = self.array_ISM_dust_MBB(rhos,betas,Ts) #(distance,freq) 
        total_intensity = np.dot(self.super_pixel_dEBV,d_emission_over_dEBV) #matrix multiplication of  (pixel,distance)X(distance,frequency) = (pixel, frequency)
        total_intensity += offsets #(pixel, frequency)
        return total_intensity

    def calculate_voxel_and_total_emission_for_reconstruction(self, parameters):
        """
        Used to make the reconstruction plots after the fits, not during the fits
        """
        offsets,rhos,betas,Ts = self.separate_parameters(parameters) 
        d_emission_over_dEBV = self.array_ISM_dust_MBB(rhos,betas,Ts) #(distance,freq) 
        voxel_intensity = np.zeros((self.super_pixel_size,self.model_nslices ,self.nfreq)) #(pixel,distance,freq) 
        for freq_index in range(self.nfreq):
            voxel_intensity[:,:,freq_index] = self.super_pixel_dEBV*d_emission_over_dEBV[:,freq_index]
        total_intensity = np.dot(self.super_pixel_dEBV,d_emission_over_dEBV) #matrix multiplication of  (pixel,distance)X(distance,frequency) = (pixel, frequency)
        total_intensity += offsets #(pixel, frequency)
        return voxel_intensity, total_intensity


    def ln_likelihood(self, parameters):
        delta_chi_square=self.calculate_chi_square(parameters)
        ln_L = -0.5*delta_chi_square
        return  ln_L 
    def ln_prior_only_boundaries(self, parameters):
        boundaries_met =  np.all(parameters>=self.parameter_boundaries[:,0]) and np.all(parameters<=self.parameter_boundaries[:,1])

        if boundaries_met == True:
            return 0
        else:
            return -np.inf
        return 
    def ln_prior(self,parameters):
        """adds the Gaussian priors on the parameters as well"""
        boundaries_met =  np.all(parameters>=self.parameter_boundaries[:,0]) and np.all(parameters<=self.parameter_boundaries[:,1])

        if boundaries_met == True:
            offsets,rhos,betas,Ts = self.separate_parameters(parameters)
            fixed_T_along_sightline=self.model_configuration_dictionary["fixed_T_along_sightline"]

            if fixed_T_along_sightline == False:
                prior = -0.5* np.sum(((Ts-self.T0)/self.sigma_T)**2)
            else:
                prior = -0.5* ((Ts[0]-self.T0)/self.sigma_T)**2
            return prior
        else:
            return -np.inf
        return 



    def parameter_type_limits(self):
        self.dictionary_of_parameter_limits = {
        'offset': np.array([-100.,50.]), # this is random now, fix it
        'rho'   : np.array([1E-6,1E3]),   # this was by looking at my paper 1 and seeing the range of values
        'beta'  : np.array([1.0,3.]),    # can make this tighter, but let's start like this for now
        'T'     : np.array([5.,100.])   # can make this tighter, but let's start like this for now
        }
    def fiducial_values_for_the_parameter_type(self):
        self.dictionary_fiducial_values = {
        'offset': 0.5, # this can be changed to be different for each offset
        'rho'   : 1E-4,  # this is just by looking at my paper 1  
        'beta'  : 1.53,  # I can change this to actually be the average across the sky later
        'T'     : self.T0    # K
        }

    def make_model_attributes(self):
        """
        This function initializes the list of parameters based on the options specified for the model
        """
        self.parameter_type_limits()
        self.fiducial_values_for_the_parameter_type()
        variable_list = []
        latex_variable_list = []
        latex_units_list = []
        parameter_boundaries = []
        fiducial_values_list = []
        param_count=0
        ### Parsing the offsets
        if self.fit_for_the_offsets == True:
            variable_list += ["Of_217","Of_353","Of_545","Of_857","Of_2998"]
            latex_variable_list += [r"$Offset_{217}$",r"$Offset_{353}$",r"$Offset_{545}$",r"$Offset_{857}$",r"$Offset_{2998}$"]
            latex_units_list += ["MJy/sr"]*self.nfreq
            parameter_boundaries += [self.dictionary_of_parameter_limits["offset"]]*self.nfreq
            fiducial_values_list+=[self.dictionary_fiducial_values["offset"]]*self.nfreq
            param_count += self.nfreq
        ### Parsing the rhos
        if self.fit_for_the_rhos == True:
            fixed_rho_along_sightline = self.model_configuration_dictionary["fixed_rho_along_sightline"]
            ### case for rho varying in each voxel
            if fixed_rho_along_sightline == False:
                variable_list += ['rho_'+str(i) for i in range(self.model_nslices)]
                latex_variable_list += [r"$\rho_{"+str(i)+"}$" for i in range(self.model_nslices)]
                latex_units_list += ["-"]*self.model_nslices
                parameter_boundaries += [self.dictionary_of_parameter_limits["rho"]]*self.model_nslices
                fiducial_values_list += [self.dictionary_fiducial_values["rho"]]*self.model_nslices
                param_count += self.model_nslices
            ### case for rho only varying in each superpixel
            else:
                variable_list += ["rho"]
                latex_variable_list += [r'$\rho$']
                latex_units_list += ["-"]
                parameter_boundaries += [self.dictionary_of_parameter_limits["rho"]]
                fiducial_values_list += [self.dictionary_fiducial_values['rho']]
                param_count += 1
        ### Parsing the betas
        if self.fit_for_the_betas == True:
            fixed_beta_along_sightline=self.model_configuration_dictionary["fixed_beta_along_sightline"]
            ### case for beta varying in each voxel
            if fixed_beta_along_sightline == False:
                variable_list += ['beta_'+str(i) for i in range(self.model_nslices)]
                latex_variable_list += [r"$\beta_{"+str(i)+"}$" for i in range(self.model_nslices)]
                latex_units_list += ["-"]*self.model_nslices
                parameter_boundaries += [self.dictionary_of_parameter_limits["beta"]]*self.model_nslices
                fiducial_values_list += [self.dictionary_fiducial_values["beta"]]*self.model_nslices
                param_count += self.model_nslices
            ### case for beta only varying in each superpixel
            else:
                variable_list += ["beta"]
                latex_variable_list += [r'$\beta$']
                latex_units_list += ["-"]
                parameter_boundaries += [self.dictionary_of_parameter_limits["beta"]]
                fiducial_values_list += [self.dictionary_fiducial_values['beta']]
                param_count += 1
        ### Parsing the Ts
        if self.fit_for_the_Ts == True:
            fixed_T_along_sightline=self.model_configuration_dictionary["fixed_T_along_sightline"]

            ### case for T varying in each voxel
            if fixed_T_along_sightline == False:
                variable_list += ['T_'+str(i) for i in range(self.model_nslices)]
                latex_variable_list += [r"$T_{"+str(i)+"}$" for i in range(self.model_nslices)]
                latex_units_list += ["K"]*self.model_nslices
                parameter_boundaries += [self.dictionary_of_parameter_limits['T']]*self.model_nslices
                fiducial_values_list += [self.dictionary_fiducial_values['T']]*self.model_nslices
                param_count += self.model_nslices
            else:
                variable_list += ["T"]
                latex_variable_list += ["T"]
                latex_units_list += ["K"]
                parameter_boundaries += [self.dictionary_of_parameter_limits["T"]]
                fiducial_values_list += [self.dictionary_fiducial_values['T']]
                param_count += 1

        self.variable_list = variable_list
        self.latex_variable_list = latex_variable_list
        self.latex_units_list = latex_units_list
        self.parameter_boundaries = np.array(parameter_boundaries)
        self.fiducial_values = np.array(fiducial_values_list)
        self.ndim = param_count



    def separate_parameters(self,parameters):
        """
        This function parses the 1D array of parameters based on the options specified for the model
        """
        ### The order of parameters is:
        ### Offsets (0 or 5, one for each frequency), rho (could be 0,1, or nslices)
        ### beta (could be 0, 1, or nslices), and T (could be nslices)
        param_count=0 
        ### Parsing the offsets
        if self.fit_for_the_offsets == True:
            offsets = parameters[0:self.nfreq]
            param_count+=self.nfreq
        else:
            offsets = self.model_configuration_dictionary["fixed_offset_array"] ### 5 nr corresponding to the 5 frequencies
            ### we don't update the index count because there is no value stored in the parameters for the offsets
        ###  Parsing the rhos
        if self.fit_for_the_rhos == True:
            fixed_rho_along_sightline = self.model_configuration_dictionary["fixed_rho_along_sightline"]
            ### case for rho varying in each voxel
            if fixed_rho_along_sightline == False:
                rhos = parameters[param_count:param_count+self.model_nslices]
                param_count+=self.model_nslices
            ### case for rho only varying in each superpixel
            else:
                rho = parameters[param_count] # reading the value from the parameters since each superpixel now has one varying rho
                rhos = np.array([rho]*self.model_nslices)
                param_count+=1
        ### case for rho being fixed across the sky
        else:
            preset_rho = self.model_configuration_dictionary["preset_rho"]
            rhos = np.array([preset_rho]*self.model_nslices)
            ### we don't update the index count because there is no value stored in the parameters for rho

        ### Parsing the betas
        if self.fit_for_the_betas == True:
            fixed_beta_along_sightline = self.model_configuration_dictionary["fixed_beta_along_sightline"]
            ### case of beta varying in each voxel
            if fixed_beta_along_sightline == False:
                betas = parameters[param_count:param_count+self.model_nslices]
                param_count+=self.model_nslices
            ### case for beta only varying in each superpixel
            else:
                beta = parameters[param_count]
                betas = np.array([beta]*self.model_nslices)
                param_count+=1
        ### case for beta being fixed across the sky
        else:
            preset_beta = self.model_configuration_dictionary["preset_beta"]
            betas = np.array([preset_beta]*self.model_nslices)
            ### we don't update the index count because there is no value stored in the parameters for beta
        ### Parsing the Ts
        if self.fit_for_the_Ts == True:
            fixed_T_along_sightline = self.model_configuration_dictionary["fixed_T_along_sightline"]

            ### case of T varying in each voxel
            if fixed_T_along_sightline == False:
                Ts = parameters[param_count:param_count+self.model_nslices]
                param_count+=self.model_nslices
            else:
                T = parameters[param_count]
                Ts = np.array([T]*self.model_nslices)

        else:
            raise ValueError("What are you doing? You want fixed T? implement it")

        return offsets,rhos,betas,Ts 


    # def ln_posterior(self, scaled_parameters):
    """ This does not seem to be used for anything
    """
    #     ln_prior = self.ln_prior(scaled_parameters)
    #     if not np.isfinite(ln_prior):
    #         return -np.inf
    #     ln_like = self.ln_likelihood(scaled_parameters)
    #     return ln_prior+ln_like

