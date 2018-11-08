###====================generate_grid_comparison_plots.py====================###
###                                                                         ###
### Generates Figures 2 and S1 for "Evolution of large-scale flow from      ###
### turbulence in a two-dimensional superfluid".                            ###
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
from scipy.stats import sem
from labscript_utils.labconfig import LabConfig
import os
import StringIO
import PIL
import xlsxwriter
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerLine2D
from analysislib.johnstone_vortices_2018.parameters import *

from vortex_thermometer import thermometer_cdN, beta

## Check where to save output:
exp_config = LabConfig()
results_path = exp_config.get('paths', 'analysis_output_folder')

## Set up output to xlsx spreadsheet:
workbook = xlsxwriter.Workbook(os.path.join(results_path,'grid_data.xlsx'))
worksheet_early = workbook.add_worksheet('Early time grid data')
worksheet_late = workbook.add_worksheet('Late time grid data')
worksheet_spectra = workbook.add_worksheet('Early time energy spectra')
bold = workbook.add_format({'bold': True})
bold_center = workbook.add_format({'bold': True,'align':'center'})
boldshaded = workbook.add_format({'bold':True, 'bg_color': '#CCCCCC'})
shaded = workbook.add_format({'bg_color': '#CCCCCC','align':'center','left':1,'right':1})
leftborder = workbook.add_format({'left':1})
rightborder = workbook.add_format({'right':1})
boldleftborder = workbook.add_format({'bold': True,'left':1})
boldrightborder = workbook.add_format({'bold': True,'right':1})

## Add the headers
grid_data_table_headings = ['Grid radius (micron)',
                            'p_c',
                            'sem(p_c)',
                            'p_d',
                            'sem(p_d)',
                            'p_f',
                            'sem(p_f)',
                            'C_1',
                            'sem(C_1)',
                            'd',
                            'sem(d)',
                            'beta',
                            'beta^+',
                            'beta^-',
                            'E_k',
                            'sem(E_k)']
worksheet_early.write(0,0,'Early times grid data',bold)
worksheet_late.write(0,0,'Late times grid data',bold)
for head_i, heading in enumerate(grid_data_table_headings):
    worksheet_early.write(1,head_i,heading,bold)
    worksheet_late.write(1,head_i,heading,bold)
early_late_worksheets = [worksheet_early,worksheet_late]

worksheet_spectra.write(0,0,'Early times average energy spectra for each grid',bold)
worksheet_spectra.write(1, 0, 'Grid radius (micron):', boldshaded)
worksheet_spectra.write(2,0,'k vector',bold)

## We'll plot everything in units of microns, so convert units here
pixel_conversion = dmd_pixel_size * 1e6
xi = xi*1e6

## Get the dataframe from lyse
df = lyse.data(timeout=30)

## Throw out the bad shots
df = df[df['review_vortex_locations','Good']==True]

## Segregate the dataframe into "early times" and "late times" (df_l and df_h for low and high hold times)
df_l = df[df['vortex_spoon_wait_time']<2]
df_h = df[df['vortex_spoon_wait_time']>3.9]


### Set up the figure windows and axes for this grid ###
### Fig. 2:
fig = figure("Grid size",figsize=(2.25,3.45),  facecolor = (1,1,1,0))

## We use gridspec to arrange the subfigures
gs = GridSpec(6,1, height_ratios=[0.6, 0.8, 0.8, 0.8, 0.8, 0.8])

class_ax = fig.add_subplot(gs[1, 0])
corr_ax = fig.add_subplot(gs[2, 0])
dm_ax = fig.add_subplot(gs[3,0])
temp_ax = fig.add_subplot(gs[4,0])
e_ax = fig.add_subplot(gs[5, 0])
grid_ax = fig.add_subplot(gs[0, 0], axisbg=(1,1,1,0))
imshow_grid = gridspec.GridSpecFromSubplotSpec(1,5, gs[0, 0], wspace=0, hspace=0.2)
imshow_ax = []
for i in range(5):
    imshow_ax.append(fig.add_subplot(imshow_grid[0,i]))

### Fig. S1:
fig_specs = figure("Energy spectra", figsize = (4.75,3.5))
spectrum_ax = fig_specs.add_subplot(111)

## set up colours for the DMD images and the energy spectra
dmd_cmap = matplotlib.colors.ListedColormap([(0,0,0),(1,1,1)], name='dmd_map', N=2)
e_spec_marker_colors = iter(['red','gold','turquoise','mediumblue','darkviolet'])

### We want to loop over the two dataframes corresponding to short times and long times
### We use enumerate so that we can easily check which we're doing this round:
### i == 0 for short times and i ==1 for long times
for i, df_i in enumerate([df_l, df_h]):
    ### Create empty lists in which we can later store the average value and the uncertainty of each measurement for each grid
    grid_radii = []
    
    clustered_fraction = []
    u_clustered_fraction = []
    dipole_fraction = []
    u_dipole_fraction = []
    free_fraction = []
    u_free_fraction = []
    
    first_correlations = []
    u_first_correlations = []
    energy_per_n = []
    u_energy_per_n = []
    e_specs = []
    dipole_moment = []
    u_dipole_moment = []

    mean_vortex_number = []

    spectra_spreadsheet_column = 0
    ## Now, group the data by grid size, so we can calculate the average value of each measurement for each grid size
    for j, (grid_radius, group) in enumerate(df_i.groupby('vortex_beam_radius')):
        ## First, add the grid size to the list of grid sizes
        grid_radii.append(grid_radius * pixel_conversion)
        
        current_grid_vortex_numbers = group[('vortex_signs_classification_and_statistics','N_v')]
        mean_vortex_number.append(mean(current_grid_vortex_numbers))
        
        ## We want the fractional populations, so before averaging,
        ## each absolute number of clusters has to be divided by the
        ## total number of vortices in that shot.
        current_grid_clustered_numbers = group[('vortex_signs_classification_and_statistics','N_c')]
        current_grid_clustered_fraction = current_grid_clustered_numbers/current_grid_vortex_numbers
        clustered_fraction.append(mean(current_grid_clustered_fraction))
        u_clustered_fraction.append(sem(current_grid_clustered_fraction))
        
        current_grid_dipole_numbers = group[('vortex_signs_classification_and_statistics','N_d')]
        current_grid_dipole_fraction = current_grid_dipole_numbers/current_grid_vortex_numbers
        dipole_fraction.append(mean(current_grid_dipole_fraction))
        u_dipole_fraction.append(sem(current_grid_dipole_fraction))
        
        current_grid_free_numbers = group[('vortex_signs_classification_and_statistics','N_f')]
        current_grid_free_fraction = current_grid_free_numbers/current_grid_vortex_numbers
        free_fraction.append(mean(current_grid_free_fraction))
        u_free_fraction.append(sem(current_grid_free_fraction))
        
        current_grid_first_correlations = group[('vortex_signs_classification_and_statistics','C_1')]
        first_correlations.append(mean(current_grid_first_correlations))
        u_first_correlations.append(sem(current_grid_first_correlations))
        
        current_grid_energy_per_n = group[('energy_spectra','KE_per_vortex')]
        energy_per_n.append(mean(current_grid_energy_per_n))
        u_energy_per_n.append(sem(current_grid_energy_per_n))
        
        current_grid_dipole_moment = group[('vortex_signs_classification_and_statistics','D')]
        dipole_moment.append(mean(current_grid_dipole_moment))
        u_dipole_moment.append(sem(current_grid_dipole_moment))
        
        ## If we are looking at early times, then we want to generate Fig. S1, the comparison of early-time energy spectra
        if i ==0:
            ## We make a list containing arrays with the energy spectra of each shot, for this current grid size
            grid_e_specs = []
            for path in group['filepath']:
                run = lyse.Run(path,no_write=True)
                e_spec = run.get_result_array('energy_spectra','KE_per_vortex_spec')
                k_vec = run.get_result_array('energy_spectra','k_vec')
                grid_e_specs.append(e_spec)
            
            ## Now take that list of arrays and make it an array itself, so that we can take its average and sem across each shot at each k-vector
            grid_e_specs = array(grid_e_specs)
            grid_e_spec = mean(grid_e_specs,axis = 0)
            u_grid_e_spec = sem(grid_e_specs,axis = 0)
            e_specs.append(grid_e_spec)
            
            ## We're going to plot this in the loop, so that we don't have to worry about storing these values for each grid size independantly
            sca(spectrum_ax)
            col_now = next(e_spec_marker_colors)
            loglog(k_vec,grid_e_spec, c = col_now,label = "%.1f $\mu$m"%(grid_radius * pixel_conversion),zorder = 2)
            fill_between(k_vec[1:-1], grid_e_spec[1:-1]-u_grid_e_spec[1:-1], grid_e_spec[1:-1]+u_grid_e_spec[1:-1],color = col_now, alpha = 0.3, zorder=1)
            
            ## Also plot the current grid scale in the same colour
            plt.axvline(x = xi*pi/(grid_radius * pixel_conversion), c = col_now,ls = ':',zorder = 0)
            
            ## Add the spectra to the workbook
            for k_vec_i in range(len(k_vec)):
                worksheet_spectra.write(k_vec_i+3,0,k_vec[k_vec_i])
                worksheet_spectra.write(k_vec_i+3, 2*spectra_spreadsheet_column+1,grid_e_spec[k_vec_i],leftborder)
                worksheet_spectra.write(k_vec_i+3, 2*spectra_spreadsheet_column+2,u_grid_e_spec[k_vec_i],rightborder)
            
            worksheet_spectra.merge_range(1, 2*spectra_spreadsheet_column+1, 1, 2*spectra_spreadsheet_column+2,grid_radius * pixel_conversion,shaded)
            worksheet_spectra.write(2,2*spectra_spreadsheet_column+1,'E(k)',boldleftborder)
            worksheet_spectra.write(2,2*spectra_spreadsheet_column+2,'sem(E(k))',boldrightborder)
            spectra_spreadsheet_column += 1
            ### While we've got a single loop of the base for loop isolated, we also want to pull out the DMD image to display for each grid
            table_data = None
            with lyse.h5py.File(group['filepath'][0],'r') as hdf5_file:
                group = hdf5_file['/devices/dmd_0']
                if 'IMAGE_TABLE' in group:
                    table_data = group['IMAGE_TABLE'][:]
            ## We want to plot the 30th frame displayed on the DMD
            frame = table_data[30].tostring()
            ## Each frame is stored as a string containing a BMP image
            ## To get this into an actual image format we dump it into a memory buffer
            tempBuff = StringIO.StringIO()
            tempBuff.write(frame)
            
            ## We can then open the image using PIL from this buffer
            tempBuff.seek(0) # need to jump back to the beginning before handing it off to PIL
            img = PIL.Image.open(tempBuff)
            
            ## Due to the geometry of the DMD pixels, the raw image looks stretched. This is fixed by resizing the image.
            img = img.resize((608,342))

            ## rather than the orientation on the DMD display.
            ## First crop so that the trap is in the centre, then rotate.
            img = img.crop((222,102,419,299))
            img = img.rotate(195)
            
            ## Now that it has been rotated, we can crop further to isolate just the trapping region
            img = img.crop((31,30,166,166))
            
            ## Add the image to the figure with the x-axis location given by the grid size
            imshow_ax[j].imshow(img, interpolation='none', cmap = dmd_cmap, zorder=-1)
            imshow_ax[j].axis('off')
            
            ## If it is the central grid image, we also want to draw some arrows to indicate the direction of movement
            if j == 2:
                dx = 8*sin(15*pi/180)
                dy = 8*cos(15*pi/180)
                dx2 = (32*cos(-15*pi/180))
                dy2 = (32*sin(-15*pi/180))
                x0=66
                y0=63
                arrow_props = dict(headwidth = 3, headlength =2, width = 0.5, shrinkA=0, shrinkB=0, color = 'yellow')
                imshow_ax[j].annotate(s='', xy= (x0+2*dx,y0+2*dy), xycoords = 'data', xytext = (x0-dx,y0-dy), textcoords = 'data', arrowprops=arrow_props,zorder = 5)
                imshow_ax[j].annotate(s='', xy= (x0-dx2 -dx, y0- dy2-dy), xycoords = 'data', xytext = (x0-dx2+2*dx, y0-dy2+2*dy), textcoords = 'data', arrowprops=arrow_props,zorder = 5)
                imshow_ax[j].annotate(s='', xy= (x0+dx2 -dx, y0+dy2 -dy), xycoords = 'data', xytext = (x0+ dx2 +2*dx, y0+ dy2 +2*dy), textcoords = 'data', arrowprops=arrow_props,zorder = 5)
        
    #### DONE LOOPING OVER GRID SIZES ####
    
    ## Set up some colours for plotting, depending on whether we're on early or late times, so that you can tell them apart
    mfc_c = colour_cluster if i == 0 else light_colour_cluster
    mfc_d = colour_dipole if i == 0 else light_colour_dipole
    mfc_f = colour_free if i == 0 else light_colour_free

    mfc_i = [0,0,0] if i == 0 else [0.5,0.5,0.5]
    
    ## Calculate the temperatures associated with the average classifications:
    beta_measured = thermometer_cdN(array(clustered_fraction), array(dipole_fraction), array(mean_vortex_number))
    beta_upper_error = thermometer_cdN(array(clustered_fraction)+array(u_clustered_fraction), array(dipole_fraction)-array(u_dipole_fraction),array(mean_vortex_number))
    beta_lower_error = thermometer_cdN(array(clustered_fraction)-array(u_clustered_fraction), array(dipole_fraction)+array(u_dipole_fraction),array(mean_vortex_number))
    
    ## Plot the temperature
    sca(temp_ax)
    errorbar(grid_radii, beta_measured, yerr = array([beta_measured - beta_lower_error, beta_upper_error - beta_measured]), c=mfc_i,mfc = mfc_i,fmt="o", ms=5, mew = 0, capsize = 1, capthick =1, ls = '-', zorder = 2-i)
    
    ## Plot the classifications
    sca(class_ax)
    errorbar(grid_radii, free_fraction, yerr=u_free_fraction, c=mfc_f, mfc=mfc_f, mec = mfc_f,  fmt="o", ms=5, mew = 0, capsize = 1, capthick =1, ls = '-', label = "Free vortices", zorder = 2-i)
    errorbar(grid_radii, dipole_fraction, yerr=u_dipole_fraction, c=mfc_d, mfc=mfc_d, mec = mfc_d, fmt="o", ms=5, mew = 0, capsize = 1, capthick =1, ls = '-', label = "Dipole vortices", zorder = 2-i)
    errorbar(grid_radii, clustered_fraction, yerr=u_clustered_fraction, c = mfc_c, mfc=mfc_c, mec = mfc_c, fmt="o", ms=5, mew = 0, capsize = 1, capthick =1, ls = '-', label = "Clustered vortices", zorder = 2-i)
    
    ## Plot the correlation function
    sca(corr_ax)
    errorbar(grid_radii, first_correlations, yerr=u_first_correlations, c=mfc_i, mfc = mfc_i,fmt="o", label = "$C_1$", ms=5, mew = 0, capsize = 1, capthick =1, ls = '-', zorder = 2-i)

    ## Plot the energy
    sca(e_ax)
    errorbar(grid_radii, energy_per_n, yerr=u_energy_per_n, c=mfc_i, mfc = mfc_i, fmt="o", label = "$E$",ms=5,  mew = 0, capsize = 1, capthick =1, ls = '-', zorder = 2-i)
    
    ## Plot the dipole moment
    sca(dm_ax)
    errorbar(grid_radii,dipole_moment,yerr = u_dipole_moment, c=mfc_i, mfc = mfc_i, fmt="o", label = "$D$",ms=5,  mew = 0, capsize = 1, capthick =1, ls = '-', zorder = 2-i)
    
    ## Add data to the spreadsheet
    
    for grid_i, radius in enumerate(grid_radii):
        spreadsheet_row = grid_i +2
        early_late_worksheets[i].write(spreadsheet_row,0,radius)
        
        early_late_worksheets[i].write(spreadsheet_row,1,clustered_fraction[grid_i])
        early_late_worksheets[i].write(spreadsheet_row,2,u_clustered_fraction[grid_i])
        
        early_late_worksheets[i].write(spreadsheet_row,3,dipole_fraction[grid_i])
        early_late_worksheets[i].write(spreadsheet_row,4,u_dipole_fraction[grid_i])
        
        early_late_worksheets[i].write(spreadsheet_row,5,free_fraction[grid_i])
        early_late_worksheets[i].write(spreadsheet_row,6,u_free_fraction[grid_i])
        
        early_late_worksheets[i].write(spreadsheet_row,7,first_correlations[grid_i])
        early_late_worksheets[i].write(spreadsheet_row,8,u_first_correlations[grid_i])
        
        early_late_worksheets[i].write(spreadsheet_row,9,dipole_moment[grid_i])
        early_late_worksheets[i].write(spreadsheet_row,10,u_dipole_moment[grid_i])
        
        early_late_worksheets[i].write(spreadsheet_row,11,beta_measured[grid_i])
        early_late_worksheets[i].write(spreadsheet_row,12,beta_measured[grid_i] - beta_upper_error[grid_i])
        early_late_worksheets[i].write(spreadsheet_row,13,beta_lower_error[grid_i] - beta_measured[grid_i])
        
        early_late_worksheets[i].write(spreadsheet_row,14,energy_per_n[grid_i])
        early_late_worksheets[i].write(spreadsheet_row,15,u_energy_per_n[grid_i])
##### DONE LOOPING OVER EARLY/LATE TIMES ####

## Done writing to the spreadsheet, so close it
workbook.close()

### Time to tidy up the graphs ###



#### Create and align labels for each panel ####
## Some common label positioning values:
figlabelx = -0.19
figlabely = 0.7
labelx = -0.15
labely = -0.16
x_min = min(grid_radii)-(max(grid_radii)-min(grid_radii))/(2*(len(grid_radii)-1))
x_max = max(grid_radii)+(max(grid_radii)-min(grid_radii))/(2*(len(grid_radii)-1))

grid_ax.yaxis.set_label_coords(labelx, 0.45)
grid_ax.xaxis.set_label_coords(0.5, labely)
grid_ax.set_title('A', fontdict = blackfont,x=figlabelx,y=0.58)

class_ax.yaxis.set_label_coords(labelx, 0.5)
class_ax.xaxis.set_label_coords(0.5, labely)
class_ax.set_title('B', fontdict = blackfont,x=figlabelx,y=figlabely)

corr_ax.yaxis.set_label_coords(labelx, 0.5)
corr_ax.xaxis.set_label_coords(0.5, labely)
corr_ax.set_title('C', fontdict = blackfont,x=figlabelx,y=figlabely)
corr_ax.spines["top"].set_visible(False)

dm_ax.yaxis.set_label_coords(labelx, 0.5)
dm_ax.xaxis.set_label_coords(0.5, labely)
dm_ax.set_title('D', fontdict = blackfont,x=figlabelx,y=figlabely)
dm_ax.spines["top"].set_visible(False)

temp_ax.yaxis.set_label_coords(labelx, 0.5)
temp_ax.xaxis.set_label_coords(0.5, labely)
temp_ax.set_title('E', fontdict = blackfont,x=figlabelx,y=figlabely)
temp_ax.spines["top"].set_visible(False)

e_ax.yaxis.set_label_coords(labelx, 0.5)
e_ax.xaxis.set_label_coords(0.5, -0.3)
e_ax.set_title('F', fontdict = blackfont,x=figlabelx,y=figlabely)
e_ax.spines["top"].set_visible(False)

## DMD images - need to turn off the axis spines etc
sca(grid_ax)
tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='off')
tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    top='off',      # ticks along the bottom edge are off
    bottom='on',         # ticks along the top edge are off
    # length = 10,
    labelbottom='off')
xticks( grid_radii )
grid_ax.spines["top"].set_visible(False)
grid_ax.spines["right"].set_visible(False)
grid_ax.spines["left"].set_visible(False)
grid_ax.spines["bottom"].set_visible(False)

sca(grid_ax)
ylabel('Grid')
xlim(x_min,x_max)

## Classification graph formatting:
sca(class_ax)
xlim(x_min,x_max)
ylabel("$p_i$")
ylim(0.1,0.7)
xticks( grid_radii )
class_ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
yticks([0.2,0.3,0.4,0.5,0.6])
class_ax.tick_params(axis='x',which='minor',bottom='off', top="off")
class_ax.xaxis.set_ticklabels([])


sca(temp_ax)
xlim(x_min,x_max)
temp_ax.invert_yaxis()
ylabel("$\\beta$")
ylim(0.8,-0.8)
locator_params(axis='y',nbins=4)
temp_ax.tick_params(axis='x',which='minor',bottom='off', top="off")
axhline(y=0, c= 'k', ls='--', zorder = -1)
xticks( grid_radii )
temp_ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
temp_ax.xaxis.set_ticklabels([])

## Correlation function graph formatting:
sca(corr_ax)
axhline(y=0, c= 'k', ls='--')
ylim(-0.5,0.5)
ylabel(r"$C_1$")
xlim(x_min,x_max)
xticks( grid_radii )
corr_ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
corr_ax.xaxis.set_ticklabels([])

# Dipole moment graph formatting:
sca(dm_ax)
xlim(x_min,x_max)
xticks( grid_radii )
yticks([0.1,0.2,0.3])
ylabel(r'$d/R_\perp$')
ylim(0.05,0.35)
dm_ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
dm_ax.xaxis.set_ticklabels([])

# Energy graph formatting:
sca(e_ax)
xlim(x_min,x_max)
xticks( grid_radii )
ylabel(r"$E_k/N_v$")
ylim(2.25,3.75)
locator_params(axis='y',nbins=7)
e_ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# since this is the bottom graph, we also need an x-axis label!
xlabel(r'Grid semi-major axis $R_G$ ($\mu m$)')

spine_bound = 2.2
## Extend spines to help segregate y-labels on the stacked plots:
class_ax.spines['top'].set_bounds(spine_bound, class_ax.viewLim.intervalx[1])
class_ax.spines['bottom'].set_bounds(spine_bound, class_ax.viewLim.intervalx[1])
corr_ax.spines['bottom'].set_bounds(spine_bound, corr_ax.viewLim.intervalx[1])
dm_ax.spines['bottom'].set_bounds(spine_bound, dm_ax.viewLim.intervalx[1])
temp_ax.spines['bottom'].set_bounds(spine_bound, temp_ax.viewLim.intervalx[1])
e_ax.spines['bottom'].set_bounds(spine_bound, e_ax.viewLim.intervalx[1])

## Add a legend to the classification graph:
leg_cluster = mlines.Line2D([], [], color=colour_cluster, marker='o',
                          ms=5, mew = 0, ls = '-', label='$p_c$')
leg_dipole = mlines.Line2D([], [], color=colour_dipole, marker='o',
                          ms=5, mew = 0, ls = '-', label='$p_d$')
leg_free = mlines.Line2D([], [], color=colour_free, marker='o',
                          ms=5, mew = 0, ls = '-', label='$p_f$')
sca(class_ax)
legend( handles=[leg_cluster,leg_dipole,leg_free],
        #loc = 2,
        bbox_to_anchor = (0.55,0.92),
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

## Fix the margins and padding between subfigures
gs.update(left=0.18, right=0.999, top = 1, bottom = 0.08, wspace=0.2, hspace = 0.015)

## Make sure we've got the correct figure active and save it!
figure("Grid size")
savefig(os.path.join(results_path,'grid_plots.pdf'),transparent = True)

#### Now for the energy spectrum figure:
sca(spectrum_ax)
ylabel(r"$E(k)/N_v$")
xlabel(r"$k\xi$")

leg = legend(loc=3,title="Transverse grid width:",fontsize = 7)
leg.get_frame().set_linewidth(0.0)

plt.axvline(x = 1, c = '0.5', ls = '-',zorder = 0)
plt.axvline(x = xi*2*pi/(31*1.25*2), c = '0.5',ls = '-',zorder = 0)

plt.text(0.8 * xi * 2*pi/(31*1.25*2), 2e-8, r'$k = \pi/R_\perp$',zorder = 1,rotation = 90)
plt.text(0.8, 2e-8, r'$k = 1/\xi$',zorder = 1,rotation = 90)


## Generate guides to the eye at power-law scaling:
minus1 = 2e-7 * k_vec**(-1)
minus3 =2e-7 *k_vec**(-3)
minus53= 1.55e-7 *k_vec**(-5./3)

plt.loglog(k_vec,minus1, c='k', ls = '-', lw=0.5, zorder = 0)
plt.loglog(k_vec,minus3, c='k', ls = ':', lw=0.5, zorder = 0)
plt.loglog(k_vec,minus53, c='k', ls = '--', lw=0.5, zorder = 0) 

## Make sure we've got the correct figure active and save it!
figure("Energy spectra")
xlim(2e-2,5)
ylim(1e-10,4e-6)
tight_layout()
savefig(os.path.join(results_path,'kinetic_energy_spectra_per_grid.pdf'),transparent = True)

## If the script is being run from outside lyse, show the graphs:
if not lyse.spinning_top: show()