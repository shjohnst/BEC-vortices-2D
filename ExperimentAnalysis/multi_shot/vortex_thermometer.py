###=========================vortex_thermometer.py===========================###
###                                                                         ###
### Calculates the vortex temperature for a given cluster and dipole        ###
### fraction and total vortex number, by finding the closest match in the   ###
### two cluster and dipole populations to Monte Carlo simulation data.      ###
### For details, see supplementary information for "Evolution of            ###
### large-scale flow from turbulence in a two-dimensional superfluid".      ###
###                                                                         ###
### Copyright (C) 2018 Shaun Johnstone and Andrew Groszek                   ###
###                                                                         ###
### This program is free software: you can redistribute it and/or modify    ###
### it under the terms of the GNU General Public License as published by    ###
### the Free Software Foundation, either version 3 of the License, or       ###
### (at your option) any later version.                                     ###
###                                                                         ###
### This program is distributed in the hope that it will be useful,         ###
### but WITHOUT ANY WARRANTY; without even the implied warranty of          ###
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           ###
### GNU General Public License for more details.                            ###
###                                                                         ###
### You should have received a copy of the GNU General Public License       ###
### along with this program.  If not, see <https://www.gnu.org/licenses/>.  ###
###                                                                         ###
### Contact: shaun.johnstone@monash.edu                                     ###
###                                                                         ###
###=========================================================================###


import os
from pylab import *
import scipy.io
from scipy.signal import savgol_filter
from labscript_utils.labconfig import LabConfig

## Get the directory containing the Monte Carlo data:
exp_config = LabConfig()
mc_data_storage_directory = exp_config.get('paths', 'simulation_data_storage')


### Load the Monte Carlo data for each vortex number, and generate smoothed versions ###


################ 4 ##################
data4 = scipy.io.loadmat(os.path.join(mc_data_storage_directory,'data_4_for_python.mat'), squeeze_me=True)

N_v_T4  = 4 # Total number of vortices

N_c_T4  = data4['N_clus']
N_d_T4  = data4['N_dip']
N_f_T4  = data4['N_free']
beta4 = data4['beta_vec']
window = 21
poly = 3
smooth_f_c4 = savgol_filter(N_c_T4/N_v_T4, window,poly)
smooth_f_d4 = savgol_filter(N_d_T4/N_v_T4,window,poly)
smooth_f_f4 = savgol_filter(N_f_T4/N_v_T4,window,poly)


################ 6 ##################
data6 = scipy.io.loadmat(os.path.join(mc_data_storage_directory,'data_6_for_python.mat'), squeeze_me=True)

N_v_T6  = 6 # Total number of vortices

N_c_T6  = data6['N_clus']
N_d_T6  = data6['N_dip']
N_f_T6  = data6['N_free']
beta6 = data6['beta_vec']
window = 21
poly = 3
smooth_f_c6 = savgol_filter(N_c_T6/N_v_T6, window,poly)
smooth_f_d6 = savgol_filter(N_d_T6/N_v_T6,window,poly)
smooth_f_f6 = savgol_filter(N_f_T6/N_v_T6,window,poly)


################ 8 ##################
data8 = scipy.io.loadmat(os.path.join(mc_data_storage_directory,'data_8_for_python.mat'), squeeze_me=True)

N_v_T8  = 8 # Total number of vortices

N_c_T8  = data8['N_clus']
N_d_T8  = data8['N_dip']
N_f_T8  = data8['N_free']
beta8 = data8['beta_vec']
window = 21
poly = 3
smooth_f_c8 = savgol_filter(N_c_T8/N_v_T8, window,poly)
smooth_f_d8 = savgol_filter(N_d_T8/N_v_T8,window,poly)
smooth_f_f8 = savgol_filter(N_f_T8/N_v_T8,window,poly)

################ 10 ##################
data10 = scipy.io.loadmat(os.path.join(mc_data_storage_directory,'data_10_for_python.mat'), squeeze_me=True)

N_v_T10  = 10 # Total number of vortices

N_c_T10  = data10['N_clus']
N_d_T10  = data10['N_dip']
N_f_T10  = data10['N_free']
beta = data10['beta_vec']
window = 21
poly = 3
smooth_f_c10 = savgol_filter(N_c_T10/N_v_T10, window,poly)
smooth_f_d10 = savgol_filter(N_d_T10/N_v_T10,window,poly)
smooth_f_f10 = savgol_filter(N_f_T10/N_v_T10,window,poly)
# sca(class_vs_t_ax)

################ 14 ##################
data14 = scipy.io.loadmat(os.path.join(mc_data_storage_directory,'data_14_for_python.mat'), squeeze_me=True)

N_v_T14  = 14 # Total number of vortices

N_c_T14  = data14['N_clus']
N_d_T14  = data14['N_dip']
N_f_T14  = data14['N_free']
beta14 = data14['beta_vec']
window = 21
poly = 3
smooth_f_c14 = savgol_filter(N_c_T14/N_v_T14, window,poly)
smooth_f_d14 = savgol_filter(N_d_T14/N_v_T14,window,poly)
smooth_f_f14 = savgol_filter(N_f_T14/N_v_T14,window,poly)

################ 20 ##################
data20 = scipy.io.loadmat(os.path.join(mc_data_storage_directory,'data_20_for_python.mat'), squeeze_me=True)

N_v_T20  = 20 # Total number of vortices

N_c_T20  = data20['N_clus']
N_d_T20  = data20['N_dip']
N_f_T20  = data20['N_free']
beta20 = data20['beta_vec']
window = 21
poly = 3
smooth_f_c20 = savgol_filter(N_c_T20/N_v_T20, window,poly)
smooth_f_d20 = savgol_filter(N_d_T20/N_v_T20,window,poly)
smooth_f_f20 = savgol_filter(N_f_T20/N_v_T20,window,poly)



## Define which MC data will be used depending on the number of vortices in the experiment:
# N_MC = 4 for 0<= N_v < 5
# N_MC = 6 for 5<= N_v < 7
# N_MC = 8 for 7<= N_v < 9
# N_MC = 10 for 9<= N_v < 12
# N_MC = 14 for 12<= N_v < 17
# N_MC = 20 for 17<= N_v
min_max_temp_N = {4:[0,5],6:[5,7],8:[7,9],10:[9,12],14:[12,17],20:[17,100]}


## Our thermometry function
def thermometer_cdN(f_c,f_d,N):
    ## Calculate \Delta_c^2 and \Delta_d^2
    d_c = (N<min_max_temp_N[4][1])[:,None]*(f_c[:,None] - smooth_f_c4[None,:])**2 + ((N>=min_max_temp_N[6][0])*(N<min_max_temp_N[6][1]))[:,None]*(f_c[:,None] - smooth_f_c6[None,:])**2 + ((N>=min_max_temp_N[8][0])*(N<min_max_temp_N[8][1]))[:,None]*(f_c[:,None] - smooth_f_c8[None,:])**2 + ((N>=min_max_temp_N[10][0])*(N<min_max_temp_N[10][1]))[:,None]*(f_c[:,None] - smooth_f_c10[None,:])**2 + ((N>=min_max_temp_N[14][0])*(N<min_max_temp_N[14][1]))[:,None]*(f_c[:,None] - smooth_f_c14[None,:])**2 + (N>=min_max_temp_N[20][0])[:,None]*(f_c[:,None] - smooth_f_c20[None,:])**2
    d_d = (N<min_max_temp_N[4][1])[:,None]*(f_d[:,None] - smooth_f_d4[None,:])**2 + ((N>=min_max_temp_N[6][0])*(N<min_max_temp_N[6][1]))[:,None]*(f_d[:,None] - smooth_f_d6[None,:])**2 + ((N>=min_max_temp_N[8][0])*(N<min_max_temp_N[8][1]))[:,None]*(f_d[:,None] - smooth_f_d8[None,:])**2 + ((N>=min_max_temp_N[10][0])*(N<min_max_temp_N[10][1]))[:,None]*(f_d[:,None] - smooth_f_d10[None,:])**2 + ((N>=min_max_temp_N[14][0])*(N<min_max_temp_N[14][1]))[:,None]*(f_d[:,None] - smooth_f_d14[None,:])**2 + (N>=min_max_temp_N[20][0])[:,None]*(f_d[:,None] - smooth_f_d20[None,:])**2
    
    ## Take the RMS
    residual = sqrt(d_c + d_d)
    
    ## Find the index of the minimum of this:
    minidx = argmin(residual,axis = 1)
    ## Return the temperature corresponding to this index
    return beta[minidx]