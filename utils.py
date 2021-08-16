import os

import astropy.io.fits as pyfits
from datetime import datetime, timedelta
import healpy as hp
import numpy as np
import time

from configparser import ConfigParser                  ### For parsing the configuration file
UNIVERSAL_CONSTANTS_LOCATION = os.environ["UNIVERSAL_CONSTANTS_LOCATION"]

parser = ConfigParser()
parser.read(UNIVERSAL_CONSTANTS_LOCATION+"/constant_configuration.cfg")
c = np.float64(parser.get('Universal Constants','c')) ## m*s-1
h = np.float64(parser.get('Universal Constants','h')) ## J*s
k = np.float64(parser.get('Universal Constants','k')) ## J*K-1
#T_cmb = np.float64(parser.get('Universal Constants','T_cmb')) ## K


def one_set_ISM_dust_MBB(nu_array, MBB_params):
    tau=MBB_params[0]
    beta_D = MBB_params[1]
    T_D=MBB_params[2]
    """
    nu_arary is expected in [GHz]
    """

    #h = 6.62607004 * 1E-34 ## J*s
    #c = 299792458 ## m*s-1
    #k = 1.38064852 * 1E-23 ## J*K-1
    #1Jy = E-26 W/m^2/Hz =   E-26 J/m^2
    # J/m^2 = E26 Jy 
    nu0 = 353.*1E9# reference nu, Hz
    scaled_nu_array = nu_array*10**9 #to transform from GHz to Hz
    x_D = h*scaled_nu_array/k/T_D

    B = 2*h*scaled_nu_array**3/c**2 / (np.exp(x_D)-1) *1E26 *1E-6# 1E26 is to transform from J/m^2/sr to Jy/sr,1E-6 to move to mega 

    #print "planck length:" , len(p)
    I_D = tau*((scaled_nu_array/nu0)**beta_D)*B
    return I_D

# Save numpy array as Primary HDU in specified FITS file
def saveArray(array, filename='temp.fits'):
    hdu = pyfits.PrimaryHDU(array)
    hdulist= pyfits.HDUList([hdu])
    hdulist.writeto(filename)
    return True


# Open fits as array
def openFits(filename, hdu=0):
    hdulist = pyfits.open(filename)
    data = hdulist[hdu].data
    hdulist.close()
    return data


# Open multiple fits files and concatenate
# Assumes files were created using run_all.sh
def openBatch(path, verbose=0):
    # Get list of files and pixels to index data into output
    files = os.listdir(path)
    if path[-1] != '/':
        path = path+'/'
    pix0 = []
    pix1 = []
    for f in files:
        string = f.split('pix')
        string = string[-1].split('-')
        pix0 += [int(string[0])]
        pix1 += [int(string[-1].split('.')[0])]
    # Initialize output array and copy over data of first file
    data = openFits(path+files[0])  #np.zeros(0)
    dims = data.shape
    print(dims)
    data = np.zeros([12288,dims[1],dims[2],dims[3]])
    data[pix0[0]:pix1[0],:,:,:] = openFits(path+files[0])
    # Loops through files and copy over data
    for i,f in enumerate(files[1:]):
        if verbose >= 1:
            print('Loading '+f)
        print(i,f,pix0[i+1],pix1[i+1],data[pix0[i+1]:pix1[i+1],:,:,:].shape,openFits(path+f).shape)
        data[pix0[i+1]:pix1[i+1],:,:,:] = openFits(path+f) #np.append(data, openFits(path+f))
    return data

 ### healpy functions
def print_NSIDE_res(NSIDE):
    print("Approximate resolution at NSIDE {} is {:.2} deg".format(
            NSIDE, hp.nside2resol(NSIDE, arcmin=True) / 60))


### smoothing functions
def get_bayestar_psf(NSIDE):
    """
    returns psf in arcmin
    """
    return hp.nside2resol(NSIDE, arcmin=True) 

def get_rad_from_arcmin(arcmin):
    #2pi...360*60
    #x-----arcmin
    rad= arcmin*np.pi/180/60
    return rad

def calculate_smoothing_psf(data_type,target_psf=10.,NSIDE=1024):
    """
    target_psf is expected in arcmin
    
    target_spf = sqrt(smoothing_psf**2+original_psf**2)
    smoothing_psf = sqrt(target_psf**2-original_psf**2)
    """
    
    if data_type=="Planck":
        #https://www.wikiwand.com/en/Planck_(spacecraft)217GH
        ## I have no idea what PSF the SFD/IRAS map has, ask Doug
        original_psf = np.array([5.5,5.0,5.0,5.0,6.1])

    elif data_type =="Bayestar":
        original_psf = get_bayestar_psf(NSIDE)
    else:
        raise ValueError("Wrong map data type!")
    smoothing_psf = (target_psf**2-original_psf**2)**0.5
    return smoothing_psf


def end_time(start_time):
    """
    Function that prints the total elapsed time since the moment the script started running
    """
    end_time = time.time()
    nr_of_seconds = end_time-start_time
    d = datetime(1,1,1) + timedelta(seconds = nr_of_seconds)

    print("Total running time was:")
    print("DAYS:HOURS:MIN:SEC")
    time_string = "%d:%d:%d:%d" % (d.day-1, d.hour, d.minute, d.second)
    print(time_string)
    return time_string


def colorblind_color_list_15():
    """
    https://somersault1824.com/tips-for-designing-scientific-figures-for-color-blind-readers/
    http://mkweb.bcgsc.ca/biovis2012/

    These values are RGB
    """
    cb_black = (0,0,0)
    cb_dark_green=(0,73/255,73/255)
    cb_blue_green=(0,146/255,146/255)
    cb_blue=(0,109/255,219/255)
    cb_medium_blue=(109/255,182/255,255/255)
    cb_light_blue=(182/255,219/255,255/255)
    cb_bright_pink=(255/255,109/255,182/255)
    cb_light_pink=(255/255,182/255,119/255)
    cb_magenta=(182/255,109/255,255/255)
    cb_purple=(73/255,0,146/255)
    cb_red=(146/255,0,0)
    cb_brown=(146/255,73/255,0)
    cb_orange=(219/255,209/255,0)
    cb_bright_green=(36/255,255/255,36/255)
    cb_yellow=(255/255,255/255,109/255)
    color_list = [cb_black,cb_dark_green,cb_blue_green,cb_blue,cb_medium_blue,
                  cb_light_blue,cb_bright_pink,cb_light_pink,cb_magenta,cb_purple,
                  cb_red,cb_brown,cb_orange,cb_bright_green,cb_yellow]
    return color_list
    
def colorblind_color_dict_15():
    color_list = colorblind_color_list_15()
    color_name_dict = dict(
        cb_black= color_list[0],
        cb_dark_green= color_list[1],
        cb_blue_green= color_list[2],
        cb_blue= color_list[3],
        cb_medium_blue= color_list[4],
        cb_light_blue= color_list[5],
        cb_bright_pink= color_list[6],
        cb_light_pink= color_list[7],
        cb_magenta= color_list[8],
        cb_purple= color_list[9],
        cb_red= color_list[10],
        cb_brown= color_list[11],
        cb_orange= color_list[12],
        cb_bright_green= color_list[13],
        cb_yellow= color_list[14])
    return color_name_dict


def gelman_rubin_convergence_test(data_from_chains):
    ### Assumes the data_from_chains array has shape (nwalkers, nruns, nparams), nparams being the nr of size dist params (11)

    ### If you look in "~/miniconda2/lib/python2.7/site-packages/pymc/diagnostics.py"
    ### you'll see that the pymc function implements a different formula for the R than it is
    ### listed at https://pymc-devs.github.io/pymc/modelchecking.html, which is the formula I
    ### implemented; in particular, it does not take the square root, and it adds to V another
    ### 1/m *B_over_n factor, which because it is so big cannot be seen in the difference between
    ### mine and theirs up until the 6th digit
    (nwalkers, nruns, nparams) = data_from_chains.shape
    print("The number of walkers is:", nwalkers)
    print("The number of runs is:", nruns)
    print("The numbers of parameters is:", nparams)
    #print("The shape of the data is:",data_from_chains.shape)
    mean_for_each_walker  =  np.mean(data_from_chains, axis = 1) # shape = (nwalker, nparams)
    #print("The mean for each walker is for the 0'th parameter is:", mean_for_each_walker[:,0])
    #print("The shape for mean_for_each_walker is:", mean_for_each_walker.shape)
    mean_of_mean = np.mean(mean_for_each_walker, axis = 0) # shape  = (nparams)
    #print("mean_of_mean is", mean_of_mean)
    #print("len of mean_of_mean is:", len(mean_of_mean))
    #print("The shape for mean_of_mean is:", mean_of_mean.shape)
    B_array = np.sum((mean_for_each_walker-mean_of_mean)**2., axis =0) * nruns/ (nwalkers- 1.) #shape = (nparams)

    #print("The shape of B_array is:",B_array.shape)
    #print("B_array is:", B_array)
    replicated_mean_for_each_walker =  np.zeros((nwalkers, nruns, nparams))
    for i in range(nwalkers):
        replicated_mean_for_each_walker[i] = np.repeat( [mean_for_each_walker[i]]  ,nruns, axis=0)
    #print("the shape of the replicated_mean_for_each_walker is:",replicated_mean_for_each_walker.shape)
    W_array = np.sum(np.sum((data_from_chains - replicated_mean_for_each_walker)**2, axis =1 )/(float(nruns)-1.) , axis = 0) / nwalkers # shape = (nparams)
    #print("The shape of the W_array is:", W_array.shape)
    Var = (float(nruns)-1.)/float(nruns) *W_array + B_array/float(nruns)
    #print("The shape of the Var array is:", Var.shape)
    R = (Var/W_array)**0.5
    #print("The shape of the R array is:", R.shape)
    # I removed the pymc, pymc3 in september 2020 because of packaging issues
    #R_pymc = pymc3.gelman_rubin(data_from_chains)

    print("my R gives:", R)
    print("my R^2 is:", R**2)
    #print("R from pymc gives:",R_pymc)
    my_R = R
    my_R_square = R**2

    return my_R,my_R_square


#####################################
#### Sky areas parameters
#####################################

def get_sky_area_parameters(sky_area, super_pixel_nside):
            
    dict_zoom={}
    if sky_area == "rho_ophiuchi":
        if super_pixel_nside == 32:
            dict_zoom['start_super_pixel']=4864
            dict_zoom['end_super_pixel']=5120
        else:
            raise ValueError("Calculate the zoom in super pixel index")
    elif sky_area == "cepheus":
        if super_pixel_nside == 128:
            dict_zoom['start_super_pixel']=88558 
            dict_zoom['end_super_pixel']=96208
        else:
            raise ValueError("Calculate the zoom in super pixel index")
    elif sky_area == "tiny_cepheus": #, l,b=(99.5,10.2)
        if super_pixel_nside == 128:
            ### the center is at 94558
            dict_zoom['start_super_pixel']=94557
            dict_zoom['end_super_pixel']=94559
        elif super_pixel_nside == 64:
            ### the center is at 23639
            dict_zoom['start_super_pixel']=23638
            dict_zoom['end_super_pixel']=23640
        elif super_pixel_nside == 32:
            ### the center is at 5909
            dict_zoom['start_super_pixel']=5908
            dict_zoom['end_super_pixel']=5910
        else:
            raise ValueError("Calculate the zoom in super pixel index")
    elif sky_area == "lower_right_tiny_cepheus":
        if super_pixel_nside == 128:
            ### the center is at 94558
            dict_zoom['start_super_pixel']=94534
            dict_zoom['end_super_pixel']=94535
        elif super_pixel_nside == 64:
            ### the center is at 23639
            dict_zoom['start_super_pixel']=23633
            dict_zoom['end_super_pixel']=23634
        elif super_pixel_nside == 32:
            ### the center is at 5909
            dict_zoom['start_super_pixel']=5908
            dict_zoom['end_super_pixel']=5909
        else:
            raise ValueError("Calculate the zoom in super pixel index")
    else:
        raise ValueError("unknown sky area")

    return dict_zoom
def get_sky_area_zoom_in_parameters(zoom_in_area):
    dict_zoom={}
    if zoom_in_area == "rho_ophiuchi":
        dict_zoom['rot'] = (0.,20.) # for the Orion region above the galactic plane
        dict_zoom['xsize']= 500
        ### assuming a nested NSIDE 1024
    elif zoom_in_area == "cepheus":
        dict_zoom['rot']=(106.5,14) 
        dict_zoom['xsize']=840
    else:
        raise ValueError("unknown zoom-in area")

    return dict_zoom

def get_cepheus_healpix_index(map_nside, nested=True):
    cepheus_bi_cloud_point = (99.5,10.2)
    cepheus_bi_cloud_point = (99.5,10.2)
    l,b= cepheus_bi_cloud_point
    return l_and_b_to_healpix_index(map_nside=map_nside,l=l,b=b,nested=nested)
def get_cepheus_plot_center():
    return ()

def l_and_b_to_healpix_index(map_nside, l,b,nested=True):
    """
    l = galactic longitude, positive towards east
    b = galaactic latitude
    """
    healpix_index=hp.ang2pix(map_nside, theta=l, phi=b, nest=nested, lonlat=True)
    return healpix_index
