###=======================plot_thermometry_curves.py========================###
###                                                                         ###
### Plots the experimental values of cluster, diploe and free vortex        ###
### populations on the Monte Carlo thermometry curves at the calculated     ###
### temperatures. Used to generate Figures 3 and S2 for "Evolution of       ###
### large-scale flow from turbulence in a two-dimensional superfluid".      ###
###                                                                         ###
### Copyright (C) 2018 Shaun Johnstone                                      ###
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

import lyse
from pylab import *
from labscript_utils.labconfig import LabConfig
import os

from scipy.stats import sem
from scipy.interpolate import griddata
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerLine2D
from mpl_toolkits.axes_grid.inset_locator import InsetPosition


from vortex_thermometer import (thermometer_cdN, beta,
    smooth_f_c4, smooth_f_d4, smooth_f_f4, smooth_f_c6, smooth_f_d6, smooth_f_f6,
    smooth_f_c8, smooth_f_d8, smooth_f_f8, smooth_f_c10, smooth_f_d10, smooth_f_f10,
    smooth_f_c14, smooth_f_d14, smooth_f_f14,smooth_f_c20, smooth_f_d20, smooth_f_f20,
    min_max_temp_N, data4, data6, data8, data10, data14, data20)

smooth_c = {4:smooth_f_c4, 6:smooth_f_c6, 8:smooth_f_c8, 10:smooth_f_c10, 14:smooth_f_c14, 20:smooth_f_c20}
smooth_d = {4:smooth_f_d4, 6:smooth_f_d6, 8:smooth_f_d8, 10:smooth_f_d10, 14:smooth_f_d14, 20:smooth_f_d20}
smooth_f = {4:smooth_f_f4, 6:smooth_f_f6, 8:smooth_f_f8, 10:smooth_f_f10, 14:smooth_f_f14, 20:smooth_f_f20}
raw_mc_data = {4:data4, 6:data6, 8:data8, 10:data10, 14:data14, 20:data20}
    
## "Standard colours for Fig. 3 and uncertainty panel come from parameters, as do font settings:
from analysislib.johnstone_vortices_2018.parameters import *

## Fix the legend text alignment when there is latex:
mpl.rcParams['text.latex.preview'] = True

## Check where to save output:
exp_config = LabConfig()
results_path = exp_config.get('paths', 'analysis_output_folder')

## Generate colour maps to take colours from for each plot in Fig. S2 depending on N_v:
cluster_cdict = {'red': ((0.0, 157./255, 153./255),
                        (1.0, 0, 0)),
                'green': ((0.0,199./255,199./255),
                        (1.0,69./255,69./255)),
                'blue': ((0.0,255./255,255./255),
                        (1.0,153./255,153./255))}

cluster_cmap = LinearSegmentedColormap('cluster_cmap', cluster_cdict)

dipole_cdict = {'red': ((0.0, 247./255, 247./255),
                        (1.0, 164./255, 164./255)),
                'green': ((0.0,185./255,185./255),
                        (1.0,57./255,57./255)),
                'blue': ((0.0,161./255,161./255),
                        (1.0,14./255,14./255))}

dipole_cmap = LinearSegmentedColormap('dipole_cmap', dipole_cdict)

free_cdict = {'red': ((0.0, 102./255, 102./255),
                        (1.0, 0./255, 0./255)),
                'green': ((0.0,255./255,255./255),
                        (1.0,128./255,128./255)),
                'blue': ((0.0,173./255,173./255),
                        (1.0,60./255,60./255))}

free_cmap = LinearSegmentedColormap('free_cmap', free_cdict)


## Create Fig. 3 and axes:
fig_t = figure("Thermometry", figsize = (4.75,1.93))

gs_fig = GridSpec(1,2)
thermometry_example_ax = fig_t.add_subplot(gs_fig[0,0])
thermometry_example_n14_ax = fig_t.add_subplot(gs_fig[0,1])



## Create Fig. S2 and axes:
fig = figure("Multiple thermometers", figsize = (6,5))
gs_temps = GridSpec(3,3)
gs_bottom = GridSpec(3,3)
t4_ax = fig.add_subplot(gs_temps[0,0])
t6_ax = fig.add_subplot(gs_temps[0,1])
t8_ax = fig.add_subplot(gs_temps[0,2])
t10_ax = fig.add_subplot(gs_temps[1,0])
t14_ax = fig.add_subplot(gs_temps[1,1])
t20_ax = fig.add_subplot(gs_temps[1,2])
temp_axes = {4:t4_ax,6:t6_ax,8:t8_ax,10:t10_ax,14:t14_ax, 20:t20_ax}
t_compare_ax = fig.add_subplot(gs_bottom[2,0:2])
temp_uncertainty_ax = fig.add_subplot(gs_bottom[2,2])

## Get the data from lyse and only keep the good shots:
df = lyse.data(timeout=30)
df = df[df['review_vortex_locations','Good']==True]



#### Begin thermometry ####
# Divide data into plots for the following vortex numbers:
temp_Ns = [4,6,8,10,14,20]

# For the uncertainty example panel, pick a grid size and hold time to display:
uncertainty_example_time = 1.5
uncertainty_example_grid = 16 ##Note grid size is in pixels == 9.7 micron

# For the temperature assignment panel, pick a grid size and hold time to display:
temperature_example_time = 0.5
temperature_example_grid = 19 ## Note grid size is in pixels == 11.5 micron

# Get the colourmaps finalised:

cluster_colours = iter(cluster_cmap(np.linspace(0,1,len(temp_Ns))))
dipole_colours = iter(dipole_cmap(np.linspace(0,1,len(temp_Ns))))
free_colours = iter(free_cmap(np.linspace(0,1,len(temp_Ns))))


## loop over each sub-panel of Fig. S2:
for tN in temp_Ns:
    # Create lists to hold the data that will be used in this panel:
    grid_time_vortex_numbers = []
    grid_time_cluster_fractions = []
    u_grid_time_cluster_fractions = []
    grid_time_dipole_fractions = []
    u_grid_time_dipole_fractions = []
    grid_time_free_fractions = []
    u_grid_time_free_fractions = []
    
    # Group our data by grid size first:
    for grid, gridgroup in df.groupby('vortex_beam_radius'):
        # And then by hold time:
        for time, timegroup in gridgroup.groupby('vortex_spoon_wait_time'):
            # For the current grid and hold time, get the data:
            current_grid_vortex_numbers = timegroup[('vortex_signs_classification_and_statistics','N_v')]
            # If the mean number of vortices for this data is in the range for the current plot, we will use it,
            # and take the mean and s.e.m. of the populations:
            if mean(current_grid_vortex_numbers) >= min_max_temp_N[tN][0] and mean(current_grid_vortex_numbers) <min_max_temp_N[tN][1]:
                grid_time_vortex_numbers.append(mean(current_grid_vortex_numbers))
                
                ## We want the fractional populations, so before averaging,
                ## each absolute number of clusters has to be divided by the
                ## total number of vortices in that shot.
                current_grid_clustered_numbers = timegroup[('vortex_signs_classification_and_statistics','N_c')]
                current_grid_clustered_fraction = current_grid_clustered_numbers/current_grid_vortex_numbers
                grid_time_cluster_fractions.append(mean(current_grid_clustered_fraction))
                u_grid_time_cluster_fractions.append(sem(current_grid_clustered_fraction))
                
                current_grid_dipole_numbers = timegroup[('vortex_signs_classification_and_statistics','N_d')]
                current_grid_dipole_fraction = current_grid_dipole_numbers/current_grid_vortex_numbers
                grid_time_dipole_fractions.append(mean(current_grid_dipole_fraction))
                u_grid_time_dipole_fractions.append(sem(current_grid_dipole_fraction))
                
                current_grid_free_numbers = timegroup[('vortex_signs_classification_and_statistics','N_f')]
                current_grid_free_fraction = current_grid_free_numbers/current_grid_vortex_numbers
                grid_time_free_fractions.append(mean(current_grid_free_fraction))
                u_grid_time_free_fractions.append(sem(current_grid_free_fraction))
                
                if time == temperature_example_time and grid == temperature_example_grid:
                    thermometer_example_C = mean(current_grid_clustered_fraction)
                    u_thermometer_example_C = sem(current_grid_clustered_fraction)
                    thermometer_example_D = mean(current_grid_dipole_fraction)
                    u_thermometer_example_D = sem(current_grid_dipole_fraction)
                    thermometer_example_F = mean(current_grid_free_fraction)
                    u_thermometer_example_F = sem(current_grid_free_fraction)
                    thermometer_example_N = mean(current_grid_vortex_numbers)
                    thermometer_example_N_MC = tN
                    beta_measured = thermometer_cdN(array([thermometer_example_C]), array([thermometer_example_D]),array([thermometer_example_N]))
                    # And uncertainty in temperature:
                    beta_upper_error = thermometer_cdN(array([thermometer_example_C])+array([u_thermometer_example_C]), array([thermometer_example_D])-array([u_thermometer_example_D]),array([thermometer_example_N]))
                    beta_lower_error = thermometer_cdN(array([thermometer_example_C])-array([u_thermometer_example_C]), array([thermometer_example_D])+array([u_thermometer_example_D]),array([thermometer_example_N]))
                    thermometer_example_beta = beta_measured[0]
                    sca(thermometry_example_ax)
                    # errorbar(beta_measured, thermometer_example_F,  yerr = u_thermometer_example_F, xerr = (beta_measured - beta_upper_error, beta_lower_error - beta_measured),  c = colour_free, ls = '', fmt="o", ms=5, mew = 0, capsize = 1, capthick =1)
                    errorbar(beta_measured, thermometer_example_C,  yerr = u_thermometer_example_C, xerr = (beta_measured - beta_upper_error, beta_lower_error - beta_measured), c = colour_cluster, ls = '', fmt="o",ms=5, mew = 0, capsize = 1, capthick =1)
                    errorbar(beta_measured, thermometer_example_D,  yerr = u_thermometer_example_D, xerr = (beta_measured - beta_upper_error, beta_lower_error - beta_measured), c = colour_dipole, ls = '', fmt="o", ms=5, mew = 0, capsize = 1, capthick =1)
                    
                
                
                ## Check to see if we should generate the uncertainty example panel this round:
                if time == uncertainty_example_time and grid == uncertainty_example_grid:
                    uncertainty_example_C = mean(current_grid_clustered_fraction)
                    u_uncertainty_example_C = sem(current_grid_clustered_fraction)
                    uncertainty_example_D = mean(current_grid_dipole_fraction)
                    u_uncertainty_example_D = sem(current_grid_dipole_fraction)
                    uncertainty_example_F = mean(current_grid_free_fraction)
                    u_uncertainty_example_F = sem(current_grid_free_fraction)
                    uncertainty_example_N = mean(current_grid_vortex_numbers)
                    uncertainty_example_N_MC = tN
                    
                    sca(temp_uncertainty_ax)
                    # Calculate the temperature:
                    beta_measured = thermometer_cdN(array([uncertainty_example_C]), array([uncertainty_example_D]),array([uncertainty_example_N]))
                    # And uncertainty in temperature:
                    beta_upper_error = thermometer_cdN(array([uncertainty_example_C])+array([u_uncertainty_example_C]), array([uncertainty_example_D])-array([u_uncertainty_example_D]),array([uncertainty_example_N]))
                    beta_lower_error = thermometer_cdN(array([uncertainty_example_C])-array([u_uncertainty_example_C]), array([uncertainty_example_D])+array([u_uncertainty_example_D]),array([uncertainty_example_N]))
                    
                    # Plot the data point with error bars for both fractional population and temperature:
                    errorbar(beta_measured, uncertainty_example_F,  yerr = u_uncertainty_example_F, xerr = (beta_measured - beta_upper_error, beta_lower_error - beta_measured),  c = colour_free, ls = '',fmt="o", ms=4, mew = 0, capsize = 0)
                    errorbar(beta_measured, uncertainty_example_D,  yerr = u_uncertainty_example_D, xerr = (beta_measured - beta_upper_error, beta_lower_error - beta_measured), c = colour_dipole, ls = '',fmt="o", ms=4, mew = 0, capsize = 0)
                    errorbar(beta_measured, uncertainty_example_C,  yerr = u_uncertainty_example_C, xerr = (beta_measured - beta_upper_error, beta_lower_error - beta_measured), c = colour_cluster, ls = '',fmt="o", ms=4, mew = 0, capsize = 0)
                    
                    # Draw extra points showing the upper and lower temperature bounds for both cluster and dipole populations:
                    errorbar(beta_lower_error, uncertainty_example_C-u_uncertainty_example_C, yerr = ([0],[u_uncertainty_example_C]), xerr = (beta_lower_error - beta_measured, [0]), c = light_colour_cluster, ls = '',fmt="s", ms=4, mew = 0, capsize = 0)
                    errorbar(beta_upper_error, uncertainty_example_C+u_uncertainty_example_C, yerr = ([u_uncertainty_example_C],[0]), xerr = ([0],beta_measured-beta_upper_error), c = light_colour_cluster, ls = '',fmt="D", ms=4, mew = 0, capsize = 0)
                    
                    errorbar(beta_lower_error, uncertainty_example_D+u_uncertainty_example_D, yerr = ([u_uncertainty_example_D],[0]), xerr = (beta_lower_error - beta_measured, [0]), c = light_colour_dipole, ls = '',fmt="s", ms=4, mew = 0, capsize = 0)
                    errorbar(beta_upper_error, uncertainty_example_D-u_uncertainty_example_D, yerr = ([0],[u_uncertainty_example_D]), xerr = ([0],beta_measured-beta_upper_error), c = light_colour_dipole, ls = '',fmt="D", ms=4, mew = 0, capsize = 0)
    
    ##### We now have all data for parameters where the average number of vortices falls within the range of the current sub-figure #####
    
    ### Perform temperature measurements on this data ###
    beta_measured = thermometer_cdN(array(grid_time_cluster_fractions), array(grid_time_dipole_fractions),array(grid_time_vortex_numbers))
    beta_upper_error = thermometer_cdN(array(grid_time_cluster_fractions)+array(u_grid_time_cluster_fractions), array(grid_time_dipole_fractions)-array(u_grid_time_dipole_fractions),array(grid_time_vortex_numbers))
    beta_lower_error = thermometer_cdN(array(grid_time_cluster_fractions)-array(u_grid_time_cluster_fractions), array(grid_time_dipole_fractions)+array(u_grid_time_dipole_fractions),array(grid_time_vortex_numbers))
    
    ### Select the correct sub-figure and generate the plot:
    sca(temp_axes[tN])
    errorbar(beta_measured, grid_time_free_fractions,  xerr = (beta_measured - beta_upper_error, beta_lower_error - beta_measured),yerr = u_grid_time_free_fractions,  c = colour_free, ls = '',fmt="o", ms=4, mew = 0, capsize = 1, capthick =1)
    errorbar(beta_measured, grid_time_dipole_fractions,  xerr = (beta_measured - beta_upper_error, beta_lower_error - beta_measured),yerr = u_grid_time_dipole_fractions,  c = colour_dipole, ls = '',fmt="o", ms=4, mew = 0, capsize = 1, capthick =1)
    errorbar(beta_measured, grid_time_cluster_fractions, xerr = (beta_measured - beta_upper_error, beta_lower_error - beta_measured), yerr = u_grid_time_cluster_fractions,  c = colour_cluster, ls = '',fmt="o", ms=4, mew = 0, capsize = 1, capthick =1)
    
    ## Now plot the smoothed Monte Carlo curves under the data points:
    plot(beta,smooth_f[tN], color = free_colours.next(),zorder = -1)
    plot(beta,smooth_d[tN], color = dipole_colours.next(),zorder = -1)
    plot(beta,smooth_c[tN], color = cluster_colours.next(),zorder = -1)
    
    ### If this is the data for N_v = 14, we also want to generate Fig. 3:
    if tN == 14:
        sca(thermometry_example_n14_ax)
        errorbar(beta_measured, grid_time_free_fractions,  xerr = (beta_measured - beta_upper_error, beta_lower_error - beta_measured),yerr = u_grid_time_free_fractions,  c = colour_free, ls = '',fmt="o", ms=5, mew = 0, capsize = 1, capthick =1)
        errorbar(beta_measured, grid_time_dipole_fractions,  xerr = (beta_measured - beta_upper_error, beta_lower_error - beta_measured),yerr = u_grid_time_dipole_fractions,  c = colour_dipole, ls = '',fmt="o", ms=5, mew = 0, capsize = 1, capthick =1)
        errorbar(beta_measured, grid_time_cluster_fractions,  xerr = (beta_measured - beta_upper_error, beta_lower_error - beta_measured),yerr = u_grid_time_cluster_fractions,  c = colour_cluster, ls = '',fmt="o", ms=5, mew = 0, capsize = 1, capthick =1)
        plot(beta,smooth_f[tN], color = colour_free,zorder = -1)
        plot(beta,smooth_d[tN], color = colour_dipole,zorder = -1)
        plot(beta,smooth_c[tN], color = colour_cluster,zorder = -1)
        sca(temp_axes[tN])
    

    ## Fix the axes and add a label:
    temp_axes[tN].invert_xaxis()
    ylim(0.01,0.75)
    xlim(0.65,-0.78)
    text(0.3,0.6635,'$N_v^{\mathrm{MC}} = %s$'%tN)

## Now tidy up the Fig. 3 example axis:
sca(thermometry_example_n14_ax)
ylim(0.0,0.67)
xlim(0.3,-0.75)
ylabel('$p_i$')
xlabel("$\\beta / |\\beta_{\mathrm{BKT}}|$    |    $\\beta / |\\beta_{\mathrm{EBC}}|$")
fig_t_xlabel_x = 0.286
fig_t_xlabel_y = -0.1
fig_t_ylabel_x = -0.11
fig_t_ylabel_y = 0.5
thermometry_example_n14_ax.yaxis.set_label_coords(fig_t_ylabel_x, fig_t_ylabel_y)
thermometry_example_n14_ax.xaxis.set_label_coords(fig_t_xlabel_x, fig_t_xlabel_y)

## Add sub-figure labels
title_x = -0.08
title_y = 0.9
thermometry_example_ax.set_title('A', fontdict = blackfont,x=title_x,y=title_y)
thermometry_example_n14_ax.set_title('B', fontdict = blackfont,x=title_x,y=title_y)


## sub-plot formatting ##
# turn off tick labels for shared axes
t4_ax.xaxis.set_ticklabels([])
t6_ax.xaxis.set_ticklabels([])
t8_ax.xaxis.set_ticklabels([])

t6_ax.yaxis.set_ticklabels([])
t8_ax.yaxis.set_ticklabels([])
t14_ax.yaxis.set_ticklabels([])
t20_ax.yaxis.set_ticklabels([])

# Label positioning:
ylabel_x = -0.11
xlabel_x = 0.4545
xlabel_y = -0.12


sca(t4_ax)
ylabel('$p_i$')
t4_ax.yaxis.set_label_coords(ylabel_x, 0.5)
## Add a legend to this first sub-figure:

leg_cluster = mlines.Line2D([], [], color=colour_cluster, marker='o',
                          ms=5, mew = 0, ls = '', label='$p_c$')
leg_dipole = mlines.Line2D([], [], color=colour_dipole, marker='o',
                          ms=5, mew = 0, ls = '', label='$p_d$')
leg_free = mlines.Line2D([], [], color=colour_free, marker='o',
                          ms=5, mew = 0, ls = '', label='$p_f$')
legend( handles=[leg_cluster,leg_dipole,leg_free],
        #loc = 2,
        bbox_to_anchor = (0.95,0.953),
        ncol=3,
        borderpad = 0,
        borderaxespad= 0,
        handletextpad= 0.2,
        columnspacing = 0.5,
        handlelength = 1,
        frameon=False,
        handler_map = { leg_cluster:HandlerLine2D(numpoints=1,marker_pad=0),
                        leg_dipole:HandlerLine2D(numpoints=1,marker_pad=0),
                        leg_free:HandlerLine2D(numpoints=1,marker_pad=0)}
        )

sca(t10_ax)
ylabel('$p_i$')
t10_ax.yaxis.set_label_coords(ylabel_x, 0.5)
xlabel("$\\beta / |\\beta_{\mathrm{BKT}}|$    |    $\\beta / |\\beta_{\mathrm{EBC}}|$")
t10_ax.xaxis.set_label_coords(xlabel_x, xlabel_y)

sca(t14_ax)
xlabel("$\\beta / |\\beta_{\mathrm{BKT}}|$    |    $\\beta / |\\beta_{\mathrm{EBC}}|$")
t14_ax.xaxis.set_label_coords(xlabel_x, xlabel_y)

sca(t20_ax)
xlabel("$\\beta / |\\beta_{\mathrm{BKT}}|$    |    $\\beta / |\\beta_{\mathrm{EBC}}|$")
t20_ax.xaxis.set_label_coords(xlabel_x, xlabel_y)
#### End of placing data on thermometer curves ####


#### Show all MC curves on top of each other for comparison ####

## Re-instansiate the iterators for the colours, corresponding to the same colours used in the previous section
cluster_colours = iter(cluster_cmap(np.linspace(0,1,len(temp_Ns))))
dipole_colours = iter(dipole_cmap(np.linspace(0,1,len(temp_Ns))))
free_colours = iter(free_cmap(np.linspace(0,1,len(temp_Ns))))

sca(t_compare_ax)
# For each N_v used above:
for N_v_T in temp_Ns:
    # Get the raw MC data:
    # data = scipy.io.loadmat(r'Z:\\Data\vortex_thermometry\data_%i_for_python.mat'%N_v_T, squeeze_me=True)
    data = raw_mc_data[N_v_T]
    
    # Normalise:
    N_c_T  = data['N_clus']/N_v_T
    N_d_T  = data['N_dip']/N_v_T
    N_f_T  = data['N_free']/N_v_T
    beta = data['beta_vec']

    # Plot:
    plot(beta,N_f_T, color = next(free_colours))
    plot(beta,N_d_T, color = next(dipole_colours))
    plot(beta,N_c_T, color = next(cluster_colours))
# Format and label the axes:
ax = gca()
ax.invert_xaxis()
xlim(1.25,-1.5)
ylabel('$p_i$')
xlabel("$\\beta / |\\beta_{\mathrm{BKT}}|$    |    $\\beta / |\\beta_{\mathrm{EBC}}|$")
t_compare_ax.xaxis.set_label_coords(xlabel_x, xlabel_y)
t_compare_ax.yaxis.set_label_coords(ylabel_x/2, 0.5)
#### End thermometer comparison part ####

#### Thermometer example ####
## Earlier we plotted the data and error bars, but didn't plot the MC curves or tidy up the figure ##
sca(thermometry_example_ax)

example_curve_cluster = smooth_c[thermometer_example_N_MC]
example_curve_dipole = smooth_d[thermometer_example_N_MC]
# example_curve_free = smooth_f[thermometer_example_N_MC]

# plot(beta,example_curve_free, color = colour_free)
plot(beta,example_curve_dipole, color = colour_dipole)
plot(beta,example_curve_cluster, color = colour_cluster)

## now calculate the RMS difference between the dipole and cluster measurements and their MC curves:
## This is the same code as in the function thermometer_cdN, but it usually only returns the index corresponding to the minimum RMS value
f_c = array([thermometer_example_C])
f_d = array([thermometer_example_D])
N = array([thermometer_example_N])


d_c = (N<min_max_temp_N[4][1])[:,None]*(f_c[:,None] - smooth_f_c4[None,:])**2 + ((N>=min_max_temp_N[6][0])*(N<min_max_temp_N[6][1]))[:,None]*(f_c[:,None] - smooth_f_c6[None,:])**2 + ((N>=min_max_temp_N[8][0])*(N<min_max_temp_N[8][1]))[:,None]*(f_c[:,None] - smooth_f_c8[None,:])**2 + ((N>=min_max_temp_N[10][0])*(N<min_max_temp_N[10][1]))[:,None]*(f_c[:,None] - smooth_f_c10[None,:])**2 + ((N>=min_max_temp_N[14][0])*(N<min_max_temp_N[14][1]))[:,None]*(f_c[:,None] - smooth_f_c14[None,:])**2 + (N>=min_max_temp_N[20][0])[:,None]*(f_c[:,None] - smooth_f_c20[None,:])**2

d_d = (N<min_max_temp_N[4][1])[:,None]*(f_d[:,None] - smooth_f_d4[None,:])**2 + ((N>=min_max_temp_N[6][0])*(N<min_max_temp_N[6][1]))[:,None]*(f_d[:,None] - smooth_f_d6[None,:])**2 + ((N>=min_max_temp_N[8][0])*(N<min_max_temp_N[8][1]))[:,None]*(f_d[:,None] - smooth_f_d8[None,:])**2 + ((N>=min_max_temp_N[10][0])*(N<min_max_temp_N[10][1]))[:,None]*(f_d[:,None] - smooth_f_d10[None,:])**2 + ((N>=min_max_temp_N[14][0])*(N<min_max_temp_N[14][1]))[:,None]*(f_d[:,None] - smooth_f_d14[None,:])**2 + (N>=min_max_temp_N[20][0])[:,None]*(f_d[:,None] - smooth_f_d20[None,:])**2

residual = sqrt(d_c + d_d)
residual = residual[0,:]


plot(beta,residual, color = 'black')

## Add lines showing the measured values
axhline(y=thermometer_example_C, c= colour_cluster, ls='--', dashes=(0.75, 1.25), lw = 0.5, zorder = -1)
axhline(y=thermometer_example_D, c= colour_dipole, ls='--', dashes=(0.75, 1.25), lw = 0.5, zorder = -1)
# axhline(y=thermometer_example_F, c= colour_free, ls='--', dashes=(0.75, 1.25), lw = 0.5, zorder = -1)

axvline(x=thermometer_example_beta, c= '0.4', ls='-', lw = 0.5, zorder = -1)
text(thermometer_example_beta + 0.02, -0.052,r'$\beta_m$', color = '0.4')

# generate delta_c and delta_f
delta_c = example_curve_cluster - thermometer_example_C
delta_d = example_curve_dipole - thermometer_example_D

delta_max = 0.2
normalize = mpl.colors.Normalize(vmin=0, vmax=delta_max)
# find the indices where delta_c is closest to -delta_max and delta_max, we only want to fill in this region
delta_c_gradient_min = argmin(abs(delta_c - delta_max))
delta_c_gradient_max = argmin(abs(delta_c + delta_max))

beta_delta_c_gradient = beta[delta_c_gradient_min:delta_c_gradient_max]
delta_c_for_gradient = abs(delta_c[delta_c_gradient_min:delta_c_gradient_max])
example_curve_cluster_for_gradient = example_curve_cluster[delta_c_gradient_min:delta_c_gradient_max]


delta_cluster_cdict = {'red': ((0.0, 77./255, 77./255),
                        (1.0, 1, 1)),
                'green': ((0.0,167./255,167./255),
                        (1.0,1,1)),
                'blue': ((0.0,255./255,255./255),
                        (1.0,1,1))}

delta_cluster_cmap = LinearSegmentedColormap('delta_cluster_cmap', delta_cluster_cdict)


# Now colour in the gradient to indicate the magnitude of delta. Unfortunately, there is no easy built-in function to do this, so we essentially have to colour each slice of data by hand:
for i in range(len(delta_c_for_gradient)-1):
    plt.fill_between([beta_delta_c_gradient[i], beta_delta_c_gradient[i+1]], [example_curve_cluster_for_gradient[i], example_curve_cluster_for_gradient[i+1]],[thermometer_example_C,thermometer_example_C], color=delta_cluster_cmap(normalize(delta_c_for_gradient[i])),zorder = -2)

    
## Set up inset axes to place the colorbars on the thermometry example plot in nice positions
# colorbar_c_location = [0.35,0.1,0.22,0.03]
colorbar_c_location = [0.05,0.1,0.22,0.03]
col_c_ax = axes([0,0,1,1],frameon = False,label = 'c')
col_c_ip = InsetPosition(thermometry_example_ax,colorbar_c_location)
col_c_ax.set_axes_locator(col_c_ip)

## Now add a colorbar to show the magnitude of delta:
fakefig_c = figure('Colorbar C')
cbar_c_show = imshow(array([[0,delta_max]]),cmap=delta_cluster_cmap)
sca(thermometry_example_ax)
cb_c = colorbar(cbar_c_show, cax=col_c_ax,ticks=[0,delta_max],orientation="horizontal")
cb_c.outline.set_linewidth(0.5)
cb_c.ax.set_xticklabels(cb_c.ax.get_xticklabels(), fontsize=5, y=1)

text(colorbar_c_location[0]+0.5*colorbar_c_location[2],colorbar_c_location[1]+2*colorbar_c_location[3],r"$|\Delta_c(\beta)|$", fontsize = 5,ha ='center',transform=thermometry_example_ax.transAxes)
    
# find the indices where delta_d is closest to -delta_max and delta_max, we only want to fill in this region
delta_d_gradient_min = argmin(abs(delta_d + delta_max))
delta_d_gradient_max = argmin(abs(delta_d - delta_max))

beta_delta_d_gradient = beta[delta_d_gradient_min:delta_d_gradient_max]
delta_d_for_gradient = abs(delta_d[delta_d_gradient_min:delta_d_gradient_max])
example_curve_dipole_for_gradient = example_curve_dipole[delta_d_gradient_min:delta_d_gradient_max]


delta_dipole_cdict = {'red': ((0.0, 1, 1),
                        (1.0, 1, 1)),
                'green': ((0.0,157./255,157./255),
                        (1.0,1,1)),
                'blue': ((0.0,97./255,97./255),
                        (1.0,1,1))}

delta_dipole_cmap = LinearSegmentedColormap('delta_dipole_cmap', delta_dipole_cdict)

# Now colour in the gradient to indicate the magnitude of delta. Unfortunately, there is no easy built-in function to do this, so we essentially have to colour each slice of data by hand:
for i in range(len(delta_d_for_gradient)-1):
    plt.fill_between([beta_delta_d_gradient[i], beta_delta_d_gradient[i+1]], [example_curve_dipole_for_gradient[i], example_curve_dipole_for_gradient[i+1]],[thermometer_example_D,thermometer_example_D], color=delta_dipole_cmap(normalize(delta_d_for_gradient[i])),zorder = -2)
colorbar_d_location = [0.35,0.1,0.22,0.03]
col_d_ax = axes([0,0,1,1],frameon = False, label = 'd')
col_d_ip = InsetPosition(thermometry_example_ax,colorbar_d_location)
col_d_ax.set_axes_locator(col_d_ip)


## Now add a colorbar to show the magnitude of delta:
fakefig_d = figure('Colorbar D')
cbar_d_show = imshow(array([[0,delta_max]]),cmap=delta_dipole_cmap)
sca(thermometry_example_ax)
cb_d = colorbar(cbar_d_show, cax=col_d_ax,ticks=[0,delta_max],orientation="horizontal")
cb_d.outline.set_linewidth(0.5)
cb_d.ax.set_xticklabels(cb_d.ax.get_xticklabels(), fontsize=5, y=1)

text(colorbar_d_location[0]+0.5*colorbar_d_location[2],colorbar_d_location[1]+2*colorbar_d_location[3],r"$|\Delta_d(\beta)|$", fontsize = 5,ha ='center',transform=thermometry_example_ax.transAxes)


### Now draw some arrows to indicate the deltas:

## Get delta_c at point where arrow will be drawn
delta_c_arrow_beta = -0.15
delta_c_arrow_beta_index = min(range(len(beta)), key=lambda i: abs(beta[i]-delta_c_arrow_beta))
# And draw the arrow
annotate(s='', xy=(delta_c_arrow_beta,thermometer_example_C), xytext=(delta_c_arrow_beta,example_curve_cluster[delta_c_arrow_beta_index]), arrowprops=dict(arrowstyle='<->',shrinkA=0, shrinkB=0, ec = colour_cluster, lw = 0.5))
text(delta_c_arrow_beta+0.08, example_curve_cluster[delta_c_arrow_beta_index]+0.01, r'$\Delta_c$', color = colour_cluster)

## Get delta_d at point where arrow will be drawn
delta_d_arrow_beta = delta_c_arrow_beta
delta_d_arrow_beta_index = min(range(len(beta)), key=lambda i: abs(beta[i]-delta_d_arrow_beta))
# And draw the arrow
annotate(s='', xy=(delta_d_arrow_beta,thermometer_example_D), xytext=(delta_d_arrow_beta,example_curve_dipole[delta_d_arrow_beta_index]), arrowprops=dict(arrowstyle='<->',shrinkA=0, shrinkB=0, ec = colour_dipole, lw = 0.5))
text(delta_d_arrow_beta+.08, thermometer_example_D + 0.01, r'$\Delta_d$', color = colour_dipole)


## Make a legend:
leg_cluster = mlines.Line2D([], [], color=colour_cluster, marker='o',
                          ms=5, mew = 0, ls = '', label='$p_c$')
leg_dipole = mlines.Line2D([], [], color=colour_dipole, marker='o',
                          ms=5, mew = 0, ls = '', label='$p_d$')
leg_free = mlines.Line2D([], [], color=colour_free, marker='o',
                          ms=5, mew = 0, ls = '', label='$p_f$')
                          
leg_cluster_MC = mlines.Line2D([], [], color=colour_cluster,
                          ms=5, mew = 0, ls = '-', label=r'$p_c^{\mathrm{MC}}$')
leg_dipole_MC = mlines.Line2D([], [], color=colour_dipole, marker='',
                          ms=5, mew = 0, ls = '-', label=r'$p_d^{\mathrm{MC}}$')
leg_free_MC = mlines.Line2D([], [], color=colour_free, marker='',
                          ms=5, mew = 0, ls = '-', label=r'$p_f^{\mathrm{MC}}$')

leg_rms = mlines.Line2D([], [], color='black', marker='',
                          ms=5, mew = 0, ls = '-', label=r'$R_{c,d}$')
                          
leg = legend( handles=[leg_cluster_MC,leg_cluster,leg_dipole_MC,leg_dipole,leg_rms],
        bbox_to_anchor = (0.65,1),
        # loc = 'upper left',
        ncol=3,
        borderpad = 0.15,
        borderaxespad= 0.15,
        handletextpad= 0.2,
        columnspacing = 0.5,
        labelspacing = 0.15,
        handlelength = 0.8,
        frameon=False,
        handler_map = { leg_cluster:HandlerLine2D(numpoints=1,marker_pad=0),
                        leg_dipole:HandlerLine2D(numpoints=1,marker_pad=0),
                        leg_cluster_MC:HandlerLine2D(numpoints=1,marker_pad=0),
                        leg_dipole_MC:HandlerLine2D(numpoints=1,marker_pad=0),
                        leg_rms:HandlerLine2D(numpoints=1,marker_pad=0)
                        }
        )

thermometry_example_ax.invert_xaxis()
ylabel(r'$p_i,\ R$')
ylim(0,0.67)
xlim(0.3,-0.75)
xlabel("$\\beta / |\\beta_{\mathrm{BKT}}|$    |    $\\beta / |\\beta_{\mathrm{EBC}}|$")
thermometry_example_ax.yaxis.set_label_coords(fig_t_ylabel_x, fig_t_ylabel_y)
thermometry_example_ax.xaxis.set_label_coords(fig_t_xlabel_x, fig_t_xlabel_y)


#### End thermemeter example figure formatting ####

## Also add a legend to the N=14 window:
sca(thermometry_example_n14_ax)
legend( handles=[leg_cluster_MC,leg_cluster,leg_dipole_MC,leg_dipole,leg_free_MC,leg_free],
        bbox_to_anchor = (0.65,1),
        # loc = 'upper left',
        ncol=3,
        borderpad = 0.15,
        borderaxespad= 0.15,
        handletextpad= 0.2,
        columnspacing = 0.5,
        labelspacing = 0.15,
        handlelength = 0.8,
        frameon=False,
        handler_map = { leg_cluster:HandlerLine2D(numpoints=1,marker_pad=0),
                        leg_dipole:HandlerLine2D(numpoints=1,marker_pad=0),
                        leg_free:HandlerLine2D(numpoints=1,marker_pad=0),
                        leg_cluster_MC:HandlerLine2D(numpoints=1,marker_pad=0),
                        leg_dipole_MC:HandlerLine2D(numpoints=1,marker_pad=0),
                        leg_free_MC:HandlerLine2D(numpoints=1,marker_pad=0)
                        }
        )

#### Thermometer uncertainty example ####
## Earlier we plotted the data and error bars, but didn't plot the MC curves or tidy up the figure ##
sca(temp_uncertainty_ax)
plot(beta,smooth_f[uncertainty_example_N_MC], color = colour_free)
plot(beta,smooth_d[uncertainty_example_N_MC], color = colour_dipole)
plot(beta,smooth_c[uncertainty_example_N_MC], color = colour_cluster)

temp_uncertainty_ax.invert_xaxis()
ylabel('$p_i$')
ylim(0.03,0.75)
xlim(-0.43,-0.67)
xlabel("$\\beta / |\\beta_{\mathrm{EBC}}|$")
temp_uncertainty_ax.xaxis.set_label_coords(0.5, xlabel_y)

#### End thermemeter uncertainty figure formatting ####


### Figure labels for Fig. S2: ###
sub_label_x =0.15
sub_label_y = 0.84
t4_ax.set_title('A', fontdict = blackfont,x=sub_label_x,y=sub_label_y)
t6_ax.set_title('B', fontdict = blackfont,x=sub_label_x,y=sub_label_y)
t8_ax.set_title('C', fontdict = blackfont,x=sub_label_x,y=sub_label_y)
t10_ax.set_title('D', fontdict = blackfont,x=sub_label_x,y=sub_label_y)
t14_ax.set_title('E', fontdict = blackfont,x=sub_label_x,y=sub_label_y)
t20_ax.set_title('F', fontdict = blackfont,x=sub_label_x,y=sub_label_y)
sub_label_x =0.06
t_compare_ax.set_title('G', fontdict = blackfont,x=sub_label_x/2,y=sub_label_y)
temp_uncertainty_ax.set_title('H', fontdict = blackfont,x=sub_label_x,y=sub_label_y)

## Fix up figure spacing and save:
gs_temps.update(left=0.05, right=0.999, top = 0.999, bottom = 0.2, wspace=0., hspace = 0.)
gs_bottom.update(left=0.05, right=0.999, top = 0.999, bottom = 0.06, wspace=0.25, hspace = 0)
figure("Multiple thermometers")
savefig(os.path.join(results_path,'thermometer_comparison.pdf'),transparent = True)

figure("Thermometry")
gs_fig.update(left=0.07, right=0.999, top = 0.999, bottom = 0.14, wspace=0.2, hspace = 0.)
savefig(os.path.join(results_path,'thermometry.pdf'),transparent = True)

## If the script is being run from outside lyse, show the graphs:
if not lyse.spinning_top: show()