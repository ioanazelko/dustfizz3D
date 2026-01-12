# This is the master script that runs the entire code.
import data_processing
import model
import optimizer
import utils
import sampler
import sky_analysis
#import sky_plots

##################################
##### How to generate data with mathing PSFs
##################################

### Making the smooth ebv maps
#data_processing.make_smooth_ebv()
data_processing.make_smooth_ebv(dtype="float16")

## Making the smooth planck and IRAS maps
#data_processing.make_smooth_planck()

##########################################
##### How to generate a 3D temperature map
##########################################

sky_analysis.test_sky_analysis()

