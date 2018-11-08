###==================plot_multiple_bragg_vortex_images.py===================###
###                                                                         ###
### Generates base for Figure 1 for "Evolution of large-scale flow from     ###
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

from __future__ import division

import os
import lyse
from pylab import *
from matplotlib.colors import Normalize
from copy import deepcopy
from analysislib.johnstone_vortices_2018.parameters import *
from labscript_utils.labconfig import LabConfig

## some new colours for dipole, vortex and antivortex:
plotred = (0.8, 0.3, 0.25)
plotblue = (0, 0.5, 0.85)
plotgreen = (0, 0.8, 0.4)

### redefine ROIs to be tighter:
plotROIs = deepcopy(ROIs)
for roiname, roidict in plotROIs.iteritems():
    roidict['y1']+=11
    roidict['x1']+=12
    roidict['y2'] -= 9
    roidict['x2'] -= 9

 

## Therefore have to redefine wx, wy, plot_condensate_centre too:
plotwx = plotROIs['roi0']['x2'] - plotROIs['roi0']['x1']
plotwy = plotROIs['roi0']['y2'] - plotROIs['roi0']['y1']
plot_condensate_centre = array([0.5*plotwx,0.5*plotwy])

exp_config = LabConfig()
shot_storage_directory = exp_config.get('paths', 'experiment_shot_storage')
results_path = exp_config.get('paths', 'analysis_output_folder')
## A list of files to pull the images from:
files = [r'20171004T093812_flat_trap_rb_quad_lev_240.h5',r'20171004T121555_flat_trap_rb_quad_lev_052.h5',r'20171004T144649_flat_trap_rb_quad_lev_109.h5']
full_paths = [os.path.join(shot_storage_directory,file) for file in files]



## Set up the figure environment
fig = figure('Example images',figsize=(4.75,3.575))
gs = GridSpec(3,4)

## Define the labels to use for each column of the figure:
labels0 = ['A','B','C']
labels1 = ['D','E','F']
labels2 = ['G','H','I']
labels3 = ['J','K','L']

## Loop over the three figures:
for i,file in enumerate(full_paths):
    with lyse.h5py.File(file, 'r') as h5_file:
        if "results/review_vortex_locations" in h5_file:
            for key, val in h5_file["results"]["review_vortex_locations"].iteritems():
                if key =="vortex_points":
                    points_of_interest = val[:]
        if "results/vortex_signs_classification_and_statistics" in h5_file:
            for key, val in h5_file["results"]["vortex_signs_classification_and_statistics"].iteritems():
                if key =="vortex_signs":
                    vortex_signs = val[:]
                elif key == "calculated_flow":
                    calculated_flow_field = val[:]
                elif key == "calculated_flow_vector":
                    calculated_vector_flow_field = val[:]
                elif key =="dipoles":
                    dipoles = val[:]
                elif key == "cluster_edges":
                    cluster_plot_lines = val[:]
        if "/results/eigenfringe_OD" in h5_file:
                for key, val in h5_file["results"]["eigenfringe_OD"].iteritems():
                    if key =="OD":
                        raw_OD = val[:]

    ### Now onto plotting everything ###
    
    ## First column is the OD image:
    ax = fig.add_subplot(gs[i, 0])
    
    plot_OD = array(raw_OD)
    ax.imshow(plot_OD[plotROIs['roi0']['y1']:plotROIs['roi0']['y2'],plotROIs['roi0']['x1']:plotROIs['roi0']['x2']], vmin = 0, vmax=2, cmap="gray",  interpolation="none")
    ax.autoscale(False)
    
    ax.axis('off')
    ax.text(4,8,labels0[i], fontdict = whitefont)
   
    ## Second column shows the differential Bragg signal:
    ax = fig.add_subplot(gs[i, 1])
    
    differential_OD = raw_OD[plotROIs['roi2']['y1']:plotROIs['roi2']['y2'],plotROIs['roi2']['x1']:plotROIs['roi2']['x2']] - raw_OD[plotROIs['roi1']['y1']:plotROIs['roi1']['y2'],plotROIs['roi1']['x1']:plotROIs['roi1']['x2']]
    X,Y = meshgrid(range(plotwx),range(plotwy))
    # We want to crop it outside the BEC:
    differential_OD[(X-0.5*plotwx)**2 + (Y-0.5*plotwy)**2 > (circ_radius+1)**2] = 0
    
    ax.imshow(differential_OD, vmin = -1, vmax = 1,cmap="seismic",  interpolation="none")
    ax.autoscale(False)
    
    # add a circle around the condensate:
    circ = Circle((plotwx/2.0,plotwy/2.0), circ_radius, fill=False, ls = "--")
    ax.add_patch(circ)
    
    # Put an 'x' at each vortex location:
    if len(points_of_interest):
        for point in points_of_interest:
            px = point[0][0] - plotROIs['roi0']['x1']
            py = point[0][1] - plotROIs['roi0']['y1']
            ax.scatter(px,py, color='k', marker='x', s=12)
    ax.axis('off')
    ax.text(4,8,labels1[i], fontdict = blackfont)
    
    
    ## Third column shows the projection of the calculated flow field:
    ax = fig.add_subplot(gs[i, 2])
    
    calculated_flow_field = calculated_flow_field[10:-10,11:-10]
    ax.imshow(calculated_flow_field,cmap="seismic",vmin = -0.0005, vmax = 0.0005)
    ax.autoscale(False)
    
    # add a circle around the condensate: 
    circ = Circle((plotwx/2.0,plotwy/2.0), circ_radius, fill=False, ls = "--")
    ax.add_patch(circ)
    
    ax.axis('off')
    ax.text(4,8,labels2[i], fontdict = blackfont)
    xlim(0,plotwx)
    ylim(plotwy,0)
    
    
    ## Fourth column shows the classified vortices and the flow contours:
    ax = fig.add_subplot(gs[i, 3])
    streamplot(array(range(plotwx)),array(range(plotwy)),-calculated_vector_flow_field[10:-10,11:-10,0],-calculated_vector_flow_field[10:-10,11:-10,1],density = 1.25,color = '0.6', linewidth = 0.4,arrowsize = 0.4,arrowstyle='-|>')
    
    ## Get the vortex locations into the correct coordinates for plotting:
    vortex_locations = []
    for point in points_of_interest:
            px = point[0][0] - plotROIs['roi0']['x1']
            py = point[0][1] - plotROIs['roi0']['y1']
            vortex_locations.append((px,py))
    vortex_locations_array = array(vortex_locations)
    
    ## An easy way of labelling positive vortices blue and negative vortices green:
    pos_neg_cmap = matplotlib.colors.ListedColormap([plotgreen,plotblue], name='from_list', N=2)

    ## First plot the outline of vortices that are in dipole pairs, and draw a line between them:
    for pair in dipoles:
        ax.plot([vortex_locations_array[pair[0],0],vortex_locations_array[pair[1],0]],[vortex_locations_array[pair[0],1],vortex_locations_array[pair[1],1]],c = plotred,lw = 1.5, solid_capstyle='round')
        ax.scatter(vortex_locations_array[pair[0],0],vortex_locations_array[pair[0],1],zorder=200,marker = 'o', lw = 1, s=11, c=[0,0,0,0],edgecolors=plotred)
        ax.scatter(vortex_locations_array[pair[1],0],vortex_locations_array[pair[1],1],zorder=200,marker = 'o', lw = 1, s=11, c=[0,0,0,0],edgecolors=plotred)
    
    ## Then draw all of the vortices with colours corresponding to their sign:
    ax.scatter(vortex_locations_array[:,0],vortex_locations_array[:,1],c=vortex_signs, marker = 'o', lw = 0, s=12, cmap = pos_neg_cmap, zorder = 100)

    # and finally, join the dots for the clusters
    for edge in cluster_plot_lines:
        if vortex_signs[edge[0]] == 1:
            col = plotblue
        else:
            col = plotgreen
        ax.plot([vortex_locations_array[edge[0],0],vortex_locations_array[edge[1],0]],[vortex_locations_array[edge[0],1],vortex_locations_array[edge[1],1]],c = col, lw = 1.5, solid_capstyle='round')
    
    # add a circle around the condensate: 
    circ = Circle((plotwx/2.0,plotwy/2.0), circ_radius, fill=False, ls = "--")        
    ax.add_patch(circ)

    ax.axis('off')
    ax.text(4,8,labels3[i], fontdict = blackfont)
    xlim(0,plotwx)
    ylim(0,plotwy)
    
    # This is the only panel that does not have an imshow, which means that its y-axis is upside down relative to all of the others!
    ax.invert_yaxis()
    
    
### Tidy up the figure spacing and save
gs.update(left=0.0, right=1, top = 1, bottom = 0, wspace=0., hspace = 0.)
savefig(os.path.join(results_path,'example_images.pdf'))


## If the script is being run from outside lyse, show the graphs:
if not lyse.spinning_top: show()