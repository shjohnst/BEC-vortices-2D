###==================get_vortex_locations_automatically.py==================###
###                                                                         ###
### Finds vortices in an optical density image of a Bose-Einstein           ###
### condensate, using a "blob detection" algorithm. Based on the algorithm  ###
### described in A. Rakonjac, et al., Physical Review A 93, 013607 (2016)   ###
### (DOI: https://doi.org/10.1103/PhysRevA.93.013607)                       ###
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
from matplotlib.patches import Ellipse
import lyse
from lyse import *
import cv2
from analysislib.johnstone_vortices_2018.parameters import *

run = Run(path)

with h5py.File(run.h5_path,'r') as h5_file:
    if "/results/eigenfringe_OD" in h5_file:
        for key, val in h5_file["results"]["eigenfringe_OD"].iteritems():
            if key =="OD":
                raw_OD = val[:]


X,Y = meshgrid(range(wx),range(wy))
cropped_OD = raw_OD[ROIs['roi0']['y1']:ROIs['roi0']['y2'],ROIs['roi0']['x1']:ROIs['roi0']['x2']]

sum_number = sum(cropped_OD)

figure('Cropped raw OD')
imshow(cropped_OD,vmin = 0, vmax=2, cmap="gray",  interpolation="none")
colorbar()

gaussian_blur_size = 1.4
blur = cv2.GaussianBlur(cropped_OD,(0,0),gaussian_blur_size)

figure('Gaussian Blurred OD')
imshow(blur,vmin = 0, vmax=2, cmap="gray",  interpolation="none")
colorbar()

laplacian = cv2.Laplacian(blur,cv2.CV_64F)

figure('Laplacian')
imshow(laplacian, cmap="gray",  interpolation="none")
colorbar()

crop_level = 0.09
ret,th = cv2.threshold(laplacian,crop_level,255,cv2.THRESH_BINARY)

th[(X-condensate_centre[0])**2 + (Y-condensate_centre[1])**2 > (circ_radius-2)**2] = 0
figure('Binary')
imshow(th,cmap="gray",  interpolation="none")
colorbar()


th2=cv2.convertScaleAbs(th)
img, contours, h = cv2.findContours(th2, cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
coordinates = []

figure('Detected vortices')
imshow(cropped_OD,vmin = 0, vmax=2, cmap="gray",  interpolation="none")
for cnt in contours:
    if len(cnt) > 2:
        
        M = cv2.moments(cnt)
        if M['m00'] != 0:
            cx = int(M['m10']/M['m00'])
            cy = int(M['m01']/M['m00'])
            
            dmu11 = M['mu11']/M['m00']
            dmu20 = M['mu20']/M['m00']
            dmu02 = M['mu02']/M['m00']
            
            if dmu20 != dmu02:
                # get largest eigenvector of covariance matrix
                l = sqrt(6*(dmu20+dmu02 + sqrt(4*dmu11**2 +(dmu20-dmu02)**2)))/2
                w = sqrt(6*(dmu20+dmu02 - sqrt(4*dmu11**2 +(dmu20-dmu02)**2)))/2
                theta = 0.5*arctan(2 * dmu11/(dmu20-dmu02)) + (dmu20<dmu02)* pi/2
            else:
                l = 1
                w = 1
                theta = 0
                
            if l > gaussian_blur_size*1.5:
                ## High aspect ratio ellipses are most likely due to close pairs of vortices
                print l/w
                print "l = %s"%l
                print dmu11
                print dmu20-dmu02
                theta = 0.5*arctan(2 * dmu11/(dmu20-dmu02)) + (dmu20<dmu02)* pi/2
                
                
                dx = 0.5 * l * cos(theta)
                dy = 0.5 * l * sin(theta)
                xx  = cx - dx
                yy = cy - dy
                coordinates.append([[xx+ROIs['roi0']['x1'],yy+ROIs['roi0']['y1']]])
                scatter(xx,yy, c='r')
                xx  = cx + dx
                yy = cy + dy
                coordinates.append([[xx+ROIs['roi0']['x1'],yy+ROIs['roi0']['y1']]])
                scatter(xx,yy, c='r', edgecolors  = None)
                
                circ = Ellipse((cx,cy), 2*l,2*w, angle = theta*180/pi, fill=False, ec='darkred',lw=2)
                ax = gca()
                ax.add_patch(circ)
            else:
                coordinates.append([[cx+ROIs['roi0']['x1'],cy+ROIs['roi0']['y1']]])
                
                circ = Ellipse((cx,cy), 2*l,2*w, angle = theta*180/pi, fill=False, ec='b',lw=2)
                ax = gca()
                ax.add_patch(circ)

xlim(0,wx)
ylim(wy,0)

run.save_result_array('vortex_points',coordinates)
run.save_result('N_v',len(coordinates))
run.save_result('pixel_sum',sum_number)


