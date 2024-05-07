# This is the master script that runs the entire code.
import data_processing
import model
import optimizer
import utils
import sampler
#import sky_plots

##################################
##### How to generate data with mathing PSFs
##################################

### Making the smooth ebv maps
data_processing.make_smooth_ebv()
## Making the smooth planck and IRAS maps
data_processing.make_smooth_planck()

##################################
##### How to generate the paper plots
##################################

#### Figure 7 (rho, and beta)


p= sky_plots.SkyPlots('bayestar_2019_full_sky_beta_fixed_nside_128_3D_5_steps',run_type='optimizer', nr_of_parallel_processes=32)
p.set_up_analysis()
data_dict = p.load_optimizer_sky_data()
# #### and Figure 8 (Ts), 128, and 64
# ##### need to do the ther runs at 128 and 64 that have the same steps
# p= sky_plots.SkyPlots('bayestar_2019_full_sky_beta_fixed_nside_128_3D_5_steps',run_type='optimizer', nr_of_parallel_processes=32)
# p.set_up_analysis()
# p.load_data()
# data_dict = p.load_optimizer_sky_data()
