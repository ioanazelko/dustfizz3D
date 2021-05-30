#### Created 2020.11.08 by Ioana Zelko
#### Purpose: sampler to be used in the creation of the 3D dust temperature map
#### Matches the tau,beta,T integrated over the distance bins, and compares it to the emission
#### from the Planck data





from __future__ import division,print_function

import os
import sys

from configparser import ConfigParser                  ### For parsing the configuration file
import functools                     ### To be able to pass class methods to ptemcee when threading (mupliprocess doesn't work with class in 2.7; only in python 3)
import h5py
import matplotlib as mpl             ### For controlling the plotting options
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection  # For plotting
import numpy as np
import pandas as pd                  ### To use for data processing
import ptemcee                       ### To use parallel tempering mcmc
import str2bool                      ### For converting the read value from the configuration parser from string to bool
import subprocess                    ### To be able to read the git HASH
import time                          ### To keep track of the program run time

DUST_3D_TEMPERATURE_MAP_DATA_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_DATA_LOCATION"]
DUST_3D_TEMPERATURE_MAP_CODE_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_CODE_LOCATION"]
DUST_3D_TEMPERATURE_MAP_PAPER_LOCATION = os.environ["DUST_3D_TEMPERATURE_MAP_PAPER_LOCATION"]
sys.path.insert(0, DUST_3D_TEMPERATURE_MAP_CODE_LOCATION)

import model
import plot_utils
import utils



def ln_likelihood_helper(instance, scaled_parameters):   #### Used to pass the class function to the Sampler object
    #unscale parameters!
    unscaled_parameters = scaled_parameters*np.abs(instance.fiducial_values)
    return instance.ln_likelihood(unscaled_parameters)
def ln_prior_helper(instance, scaled_parameters):
    #unscale parameters!
    unscaled_parameters = scaled_parameters*np.abs(instance.fiducial_values)
    if instance.use_priors_in_sampler == True:
        return instance.ln_prior(unscaled_parameters)
    else:
        return instance.ln_prior_only_boundaries(unscaled_parameters)


class Map3DSampler(model.Model):
    ### Fits a superpixel
    def __init__(self, sampler_configuration_parameters,model_configuration_dictionary,data_dictionary,sampler_super_pixel_specific_configuration):
        model.Model.__init__(self,model_configuration_dictionary,data_dictionary)
        self.nwalkers = sampler_configuration_parameters["nwalkers"]
        self.ntemps = sampler_configuration_parameters["ntemps"]
        self.nthreads = sampler_configuration_parameters["nthreads"]
        self.thinning = sampler_configuration_parameters["thinning"]
        self.nruns = sampler_configuration_parameters["nruns"]
        self.use_priors_in_sampler = sampler_configuration_parameters["use_priors_in_sampler"]
        self.initial_optimizer_parameters =sampler_super_pixel_specific_configuration["initial_optimizer_parameters"]
        self.super_pixel_index = sampler_super_pixel_specific_configuration["super_pixel_index"]
        self.sampler_plots_folder = sampler_super_pixel_specific_configuration["sampler_plots_folder"]
        self.sampler_data_folder = sampler_super_pixel_specific_configuration["sampler_data_folder"]
        if os.path.isdir(self.sampler_plots_folder):
            pass
        else:
            os.mkdir(self.sampler_plots_folder)

        if os.path.isdir(self.sampler_data_folder):
            pass
        else:
            os.mkdir(self.sampler_data_folder)   

    def get_initial_positions(self):
        np.random.seed(123)
        initial_positions_array = np.zeros((self.ntemps,self.nwalkers,self.ndim))
        for var_index in range(self.ndim):
            var = self.variable_list[var_index]
            var_min=self.parameter_boundaries[var_index,0]
            var_max=self.parameter_boundaries[var_index,1]
            optimizer_var= self.initial_optimizer_parameters[var_index]
            if var[0] in ["O","r","b"]:#for offsets, rhos, betas
                initial_positions_array[:,:,var_index] = np.clip(np.random.uniform(low=optimizer_var/1.01,\
                                high=optimizer_var*1.01,size=(self.ntemps,self.nwalkers)),a_min=var_min,\
                                                                                  a_max=var_max)
            elif var[0] == "T":#for T
                initial_positions_array[:,:,var_index] = np.clip(np.random.normal(optimizer_var, 1E-2*optimizer_var,\
                                                        size=(self.ntemps,self.nwalkers)),a_min=var_min,
                                                                a_max=var_max)
            else:
                raise ValueError("wrong variable type")
        self.scaled_initial_positions = initial_positions_array/np.abs(self.fiducial_values)
    def check_initial_positions(self):
        descaled_initial_positions =self.scaled_initial_positions*np.abs(self.fiducial_values) 
        likelihood_array =np.zeros((self.ntemps,self.nwalkers))
        prior_array = np.zeros((self.ntemps,self.nwalkers))
        for temp_index in range(self.ntemps):
            for walker_index in range(self.nwalkers):
                #likelihood_array[temp_index,walker_index]=self.ln_likelihood(descaled_initial_positions[temp_index,walker_index])
                likelihood_array[temp_index,walker_index]=self.calculate_chi_square_with_priors(descaled_initial_positions[temp_index,walker_index])
                prior_array[temp_index,walker_index]=self.ln_prior(descaled_initial_positions[temp_index,walker_index])
        print(likelihood_array)
        print(prior_array) 


    def sample(self):
        self.start_time = time.time()
        self.get_initial_positions()
        np.random.seed(123)
        sampler = ptemcee.Sampler(nwalkers = self.nwalkers,dim = self.ndim,
                                  logl= functools.partial(ln_likelihood_helper,self),
                                  logp= functools.partial(ln_prior_helper,self),
                                  ntemps=self.ntemps,threads=self.nthreads)
        state = sampler.run_mcmc(self.scaled_initial_positions,iterations=self.nruns,thin=self.thinning)

        scaled_values=sampler.chain[0]
        data_dict={}
        data_dict["chain_unscaled_values"]=scaled_values*np.abs(self.fiducial_values)
        #self.eoc_unscaled_values = self.chain_unscaled_values[:,-1,:]
        sampler_acceptance_fraction = sampler.acceptance_fraction
        data_dict["mean_acceptance_fraction"]=np.mean(sampler_acceptance_fraction)
        print("The mean of the acceptance fraction was: "+ str(np.mean(sampler_acceptance_fraction)) +"\n")
        data_dict["chain_log_likelihood"] = sampler.loglikelihood[0]
        data_dict["chain_log_prob"] = sampler.logprobability[0]
        print("length of thinned chain is",len(data_dict["chain_unscaled_values"]))
        my_R,my_R_square = utils.gelman_rubin_convergence_test(data_dict["chain_unscaled_values"])
        self.save_run_data(data_dict)
        run_time = utils.end_time(self.start_time)
        

    def save_run_data(self,data_dict):
        with h5py.File(self.sampler_data_folder+"/"+str(self.super_pixel_index)+"_sampler_run.hdf5", "w") as f:
            for key in data_dict.keys():
                f.create_dataset(key,data=data_dict[key],dtype='float64')
            f.close()
            
    def load_run_data(self):
        data_dict = {}
        with h5py.File(self.sampler_data_folder+"/"+str(self.super_pixel_index)+"_sampler_run.hdf5", "r") as g:
            for key in list(g.keys()):
                data_dict[key]=g[key][()]
            g.close()
        return data_dict
    def process_sampler_data(self):
        throw_index=int(self.nruns/self.thinning/2)
        data = self.chain_unscaled_values[:,throw_index:,:].reshape((self.nwalkers*throw_index,self.ndim))
        df = pd.DataFrame(data)
        df.columns=self.variable_list
        quantile_50=df.quantile(.50).values
        quantile_84=df.quantile(.84135).values ## to get the upper 34.135% thresholds
        quantile_16=df.quantile(.15865).values ## to get the lower 34.135% thresholds
        upperlim=quantile_84-quantile_50
        lowerlim=quantile_50-quantile_16
        for fiducial,var,center,up,down in list(zip( self.initial_optimizer_parameters , self.variable_list,quantile_50, upperlim,lowerlim)):
            print(var,'{:.1e}'.format(fiducial),'{:.1e}'.format(center),'{:.1e}'.format(up),'{:.1e}'.format(down))
        return quantile_50, upperlim, lowerlim
    def plot_chain_corner_plots(self,data_dict,throw_away_burnout=False,bins=100,plot_countours=False):
        print("Plotting the eoc corner plot")
        chain_data= data_dict["chain_unscaled_values"]
        if throw_away_burnout == False:
            data = chain_data.reshape((self.nwalkers*int(self.nruns/self.thinning),self.ndim))
        else:
            throw_index=int(self.nruns/self.thinning/2)

            data0 =  chain_data[:,throw_index:,:]
            data = data0.reshape((self.nwalkers*int(self.nruns/self.thinning/2),self.ndim))


        fig = plot_utils.plot_parameters_corner(data, self.latex_variable_list, self.initial_optimizer_parameters,bins=bins,plot_countours=plot_countours)
        if throw_away_burnout == False:
            fig.savefig(self.sampler_plots_folder+"/corner_chains_burnout_included.jpg")
        else:
            fig.savefig(self.sampler_plots_folder+"/corner_chains_burnout_not_included.jpg")
             #if self.save_plots_to_paper == True:
             #    fig.savefig(CMB_PAPER_LOCATION+"/"+self.MCMC_run_configuration_name+".pdf",bbox_inches = 'tight',pad_inches = 0.01)
        plt.clf()
        plt.close()
        print("Finish plotting the chains corner plot")