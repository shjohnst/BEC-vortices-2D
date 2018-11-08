###=======================review_vortex_locations.py========================###
###                                                                         ###
### Takes the results of "get_vortex_locations_automatically.py"            ###
### and displays the vortex locations on top of the optical density image.  ###
### Human intervention can optionally add or subtract vortices (and hence   ###
### adjust positions) before saving these locations, or reject a shot       ###
### entirely if the image is too distorted. It is the final saved locations ###
### of this script that are used in further analysis of the vortices.       ###
###                                                                         ###
### Copyright (C) 2018  Shaun Johnstone                                     ###
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
from pylab import *
from scipy import ndimage
from matplotlib.patches import Circle
import lyse
from lyse import *
from review_vortex_locations_figure_handler import VortexPoints
from analysislib.johnstone_vortices_2018.parameters import *

lyse.routine_storage.run = Run(lyse.path)

delay_results_return()
register_plot_class("Selected vortex locations", VortexPoints)

points_of_interest = []
auto_points_of_interest = []

with h5py.File(lyse.routine_storage.run.h5_path,'r') as h5_file:
    if "/results/eigenfringe_OD" in h5_file:
        for key, val in h5_file["results"]["eigenfringe_OD"].iteritems():
            if key =="OD":
                raw_OD = val[:]
    if "/results/get_vortex_locations_automatically" in h5_file:
        for key, val in h5_file["results"]["get_vortex_locations_automatically"].iteritems():
            if key =="vortex_points":
                auto_points_of_interest = val[:]
                
                
differential_OD = raw_OD[ROIs['roi2']['y1']:ROIs['roi2']['y2'],ROIs['roi2']['x1']:ROIs['roi2']['x2']] - raw_OD[ROIs['roi1']['y1']:ROIs['roi1']['y2'],ROIs['roi1']['x1']:ROIs['roi1']['x2']]


figure("Optical density reference")
imshow(raw_OD, vmin = 0, vmax=1.5, cmap="gray",  interpolation="none")
colorbar()
title("Reference image")
xlim(ROIs['roi0']['x1'],ROIs['roi0']['x2'])
ylim(ROIs['roi0']['y2'],ROIs['roi0']['y1'])
circ = Circle((ROIs['roi0']['x1']+wx/2.0,ROIs['roi0']['y1']+wy/2.0), circ_radius, fill=False, ls = "--", ec="Red")
ax = gca()
ax.add_patch(circ)

figure("Selected vortex locations")
imshow(raw_OD, vmin = 0, vmax=1.5, cmap="gray",  interpolation="none")
colorbar()
title("Select/review vortex locations")
xlim(ROIs['roi0']['x1'],ROIs['roi0']['x2'])
ylim(ROIs['roi0']['y2'],ROIs['roi0']['y1'])
circ = Circle((ROIs['roi0']['x1']+wx/2.0,ROIs['roi0']['y1']+wy/2.0), circ_radius, fill=False, ls = "--", ec="Red")
ax = gca()
ax.add_patch(circ)

## Plot the automatically found vortex locations
if len(auto_points_of_interest):
    points_of_interest = auto_points_of_interest
    for point in points_of_interest:
        plot(*point[0], color='g', marker='o', label = "POI", picker = 2.0)


figure("Differential signal")
imshow(differential_OD, vmin = -1, vmax = 1,cmap="RdBu",  interpolation="none")
colorbar()
title("Differential Bragg signal")
circ = Circle((wx/2.0,wy/2.0), circ_radius, fill=False, ls = "--")
ax = gca()
ax.add_patch(circ)
xlim(0,wx)
ylim(wy,0)