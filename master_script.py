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




def make_smooth_maps():
	#### Smooth the EBV data
	start_time = time.time()
    
    #### Making the smooth ebv maps; this takes 1h 10 min , and needs a lot of RAM, so better be done with a machine that has 
    ### more than 50GB of RAM
    data_proc_dict = {'bayestar_version':"bayestar2019",'bayestar_nside':1024,'planck_version':"albert",'planck_nside':1024}
    data_proc= data_processing.DataProcessing(data_proc_dict)
    data_proc.smooth_ebv()
    time_string = utils.end_time(start_time)

