

from __future__ import division,print_function

import os
import sys


import functools                     ### To be able to pass class methods to ptemcee when threading (mupliprocess doesn't work with class in 2.7; only in python 3)
import matplotlib as mpl             ### For controlling the plotting options
import numpy as np
from scipy.optimize import  minimize ### To be used for the optimizer
import time                          ### To keep track of the program run time


import model
import utils


class Optimizer(model.Model):
    def __init__(self, optimizer_configuration_parameters,model_configuration_dictionary,data_dictionary):
        model.Model.__init__(self, model_configuration_dictionary,data_dictionary)
        self.optimizer_method=optimizer_configuration_parameters["optimizer_method"]
        self.scale_variables =optimizer_configuration_parameters["scale_variables"]
        self.use_bounds = optimizer_configuration_parameters["use_bounds"]
        self.use_gradient = optimizer_configuration_parameters["use_gradient"]
        self.use_hessian = optimizer_configuration_parameters["use_hessian"]
        self.use_priors = optimizer_configuration_parameters["use_priors"]
        self.print_output = optimizer_configuration_parameters["print_output"]

         
    def run_optimizer(self):

        options_dict = {"maxiter" : 300000, "ftol" : 1E-9}
        ### gtol parameter  is the threshold that all gradients must be below for termination

        
        if self.scale_variables ==  True:
            def scaled_chi_square(scaled_parameters):
                unscaled_parameters = scaled_parameters*np.abs(self.fiducial_values)
                return self.calculate_chi_square(unscaled_parameters)

            def scaled_chi_square_gradient(scaled_parameters):
                unscaled_parameters = scaled_parameters*np.abs(self.fiducial_values)
                chi_s , grad = self.calculate_chi_square_and_its_gradient(unscaled_parameters)
                return chi_s , grad*np.abs(self.fiducial_values)
            
            def scaled_chi_square_with_priors(scaled_parameters):
                unscaled_parameters = scaled_parameters*np.abs(self.fiducial_values)
                return self.calculate_chi_square_with_priors(unscaled_parameters)
            
            def scaled_chi_square_with_priors_and_its_gradient(scaled_parameters):
                unscaled_parameters = scaled_parameters*np.abs(self.fiducial_values)
                chi_s , grad = self.calculate_chi_square_with_priors_and_its_gradient(unscaled_parameters)
                return chi_s , grad*np.abs(self.fiducial_values)

        if self.scale_variables == True:
            initial_values = self.fiducial_values/np.abs(self.fiducial_values) 
            boundaries = (self.parameter_boundaries.T/np.abs(self.fiducial_values)).T
            if self.use_gradient ==  True:
                if self.use_priors == True:
                    fun_for_min = scaled_chi_square_with_priors_and_its_gradient
                else:
                    fun_for_min = scaled_chi_square_gradient
            else:
                if self.use_priors == True:
                    fun_for_min = scaled_chi_square_with_priors
                else:
                    fun_for_min = scaled_chi_square
        else:
            initial_values = self.fiducial_values
            boundaries = self.parameter_boundaries
            ## Run optimizer without scaling variables
            #print("Running optimizer without scaling variables")
            if self.use_gradient ==  True:
                if self.use_priors == True:
                    fun_for_min = self.calculate_chi_square_with_priors_and_its_gradient
                else:
                    fun_for_min = self.calculate_chi_square_and_its_gradient
            else:
                if self.use_priors == True:
                    fun_for_min = self.calculate_chi_square_with_priors
                else:
                    fun_for_min= self.calculate_chi_square
        if self.use_bounds == True:
            ret= minimize(fun=fun_for_min,x0=initial_values,\
                      bounds=boundaries,method=self.optimizer_method, jac=self.use_gradient)
        else:
            ret= minimize(fun=fun_for_min,x0=initial_values,\
                     method=self.optimizer_method, jac=self.use_gradient)
          
        popt = ret.x
        if self.scale_variables ==True:
            final_opt=popt*np.abs(self.fiducial_values)
        else:
            final_opt = popt

        if self.use_gradient == True:
            optimized_function_result, gradient  = fun_for_min(popt)
        else:
            optimized_function_result = fun_for_min(popt)
        final_chi_square = self.calculate_chi_square(final_opt)
        if self.print_output == True:
            print("I ran the optimizer")
            print("---------------------------")
            print("Messages from the optimizer")
            print("Success: ", ret.success)
            print("Message: ", ret.message)
            print("Nr of iterations ", ret.nit)
            print("Nr of function evaluations: ",ret.nfev)
            print("---------------------------")
            print("Final Optimized Function Result:",optimized_function_result)
            print("Final Chi Square is: ", final_chi_square)
            print("Chi Square per pixel is: ",final_chi_square/self.super_pixel_size)
            print("Chi Square per pixel per frequency is: ",final_chi_square/self.super_pixel_size/self.nfreq)
            offsets,rhos,betas,Ts  = self.separate_parameters(final_opt)
            print("offsets are ", offsets)
            print("rhos are ", rhos)
            print("betas are ", betas)
            print("temperatures are", Ts)
            print("descaled optimized values ", final_opt)
        return final_opt, optimized_function_result, final_chi_square