###==============================parameters.py==============================###
###                                                                         ###
### Provides experimental parameters for use in vortex analysis scripts.    ###    
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


from numpy import array
import matplotlib
## Regions of interest ##

ROIs = {'roi2': {'y1': 228, 'x2': 358, 'x1': 270, 'y2': 315},
'roi0': {'y1': 180, 'x2': 313, 'x1': 225, 'y2': 267},
'roi1': {'y1': 133, 'x2': 268, 'x1': 180, 'y2': 220}}
wx = ROIs['roi0']['x2'] - ROIs['roi0']['x1']
wy = ROIs['roi0']['y2'] - ROIs['roi0']['y1']
circ_radius = 31.0
condensate_centre = array([0.5*wx,0.5*wy])


## Physical constants ##
hbar = 1.0545718e-34
mRb87 = 1.443160648e-25
pixel_size = 1.25e-6
dmd_pixel_size = 0.605e-6
xi = 0.8e-6

### Plotting settings ###

## Colours
colour_cluster = [0, 90./255, 200./255]
colour_dipole = [230./255, 80./255, 20./255]
colour_free = [0, 190./255, 90./255]
colour_532 = [101./255,255./255, 0]
light_colour_cluster = [77./255, 167./255, 255./255]
light_colour_dipole = [255./255, 157./255, 97./255]
light_colour_free = [77./255, 255./255, 167./255]

# Global font setting
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 7}

matplotlib.rc('font', **font)
blackfont = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'bold',
        'size': 9,
        }
        
whitefont = {'family': 'sans-serif',
        'color':  'white',
        'weight': 'bold',
        'size': 9,
        }
scalebarfont = {'family': 'sans-serif',
        'color':  'white',
        'weight': 'normal',
        'size': 7,
        }
tinyfont = font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 5.5}