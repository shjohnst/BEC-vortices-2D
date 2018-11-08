###============================energy_spectra.py============================###
###                                                                         ###
### Calculates the incompressible kinetic energy spectrum and total         ###
### incompressible kinetic energy due to vortices with known signs and      ###
### locations in a planar Bose-Einstein condensate.                         ###
###                                                                         ###
### Copyright (C) 2018  Shaun Johnstone and Andrew Groszek                  ###
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
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
from analysislib.johnstone_vortices_2018.parameters import *



Kappa = 2*np.pi*hbar/mRb87
units = 1.0e-6
points_of_interest = []
vortex_signs = []


### GET EXPERIMENTAL DATA ###
run = lyse.Run(lyse.path)
with h5py.File(run.h5_path,'r') as h5_file:
        if "results/review_vortex_locations" in h5_file:
            for key, val in h5_file["results"]["review_vortex_locations"].iteritems():
                if key =="vortex_points":
                    points_of_interest = val[:]
        if "results/vortex_signs_classification_and_statistics" in h5_file:
            for key, val in h5_file["results"]["vortex_signs_classification_and_statistics"].iteritems():
                if key =="vortex_signs":
                    vortex_signs = val[:]
                    
### Get vortex coordinates in microns from centre of trap.
if len(points_of_interest):
    vortex_locations = []
    for point in points_of_interest:
            px = pixel_size * (point[0][0] - ROIs['roi0']['x1'] - wx/2.)
            py = pixel_size * (point[0][1] - ROIs['roi0']['y1'] - wy/2.)
            vortex_locations.append((px,py))           
### END EXPERIMENTAL DATA COLLECTION ###


# Define function to calculate radial average of arbitrary 2D function F(x,y)
def radialaverage(F,X,Y,bin_vec,bin_width):

	R = np.sqrt(X**2 + Y**2)

	F_r = []

	for r in bin_vec:

		# Find all indices in this bin
		ind = (R < r + bin_width) & (R >= r)

		F_r.append( np.mean(F[ind]) );

	F_r = np.array(F_r)

	return F_r

# =======================================================
# Create grids
# =======================================================

# INPUT:

L     = 400 * units # System size (um) (2L = domain length)
N     = 2048#  # Grid resolution (NxN pixels)

# Spatial grid
dx    = 2.*L/N # Spacing
dy    = dx

# 1D grids
x     = np.linspace(-L,L-dx,N)
y     = x

# 2D grids
X, Y  = np.meshgrid(x,y) # Might need y, x not x, y?

# K-space grid
dk    = np.pi/L
k_max = np.pi/dx
# 1D grids
k_x   = np.linspace(-k_max, k_max-dk, N)
k_x   = np.fft.fftshift(k_x)
k_x = 2*np.pi*np.fft.fftfreq(N,dx)

k_y   = k_x

# 2D grids
K_y, K_x = np.meshgrid(k_x, k_y) # y, x not x, y?

# Need |k| matrix for energy spectrum calculation
K_mag = np.sqrt(K_x**2 + K_y**2)

# =======================================================
# INPUT:
# Density profile and velocity field
# =======================================================
R      = circ_radius * pixel_size
trap_edge = 2 * units #(um)
n0 = 1.0
n_atoms = 1e5
modpsisq = np.zeros([N,N])
modpsisq[X**2 + Y**2 < R**2] = 1
modpsisq = n0 * gaussian_filter(modpsisq,trap_edge/dx)
modpsi = np.sqrt(modpsisq)
# Imprint vortices... Multiply modpsi through by r/sqrt(r^2 + 2*xi^2) for each vortex


def vortex_core(x0,y0):
    sigma = np.sqrt((X-x0)**2 + (Y-y0)**2)/xi
    Lam = 0.82475449
    
    chi = sigma/np.sqrt(sigma**2 + Lam**(-2))
    return chi

def single_vortex_flow(vortex_loaction):
        px = vortex_loaction[0]
        py = vortex_loaction[1]
        Rx = (X -px )
        Ry = (Y - py)
        r2 = Rx**2 + Ry**2
        
        # vortex_radius = sqrt(r2) # want this for the dipole moment calc later
        fx = (hbar/ (mRb87 * r2)) * -Ry
        fy = (hbar/ (mRb87 * r2)) * Rx

        projection = -(fx+fy)
        vector_flow = np.dstack((fx,fy))
        
        # now for the image vortex. First, get its location
        pox = px
        poy = py
        por = np.sqrt(pox**2 + poy**2)
        poa = np.arctan2(poy,pox)
        # the radius of the image vortex is then:
        ir = R**2 / por
        iox = ir * np.cos(poa)
        ioy = ir * np.sin(poa)
        # and finally, add the origin of the BEC again to get back into our base coordinate system
        ix = iox
        iy = ioy
        
        Rx = (X - ix)
        Ry = (Y - iy)
        r2 = Rx**2 + Ry**2
        
        fx = (hbar/ (mRb87 * r2)) * -Ry
        fy = (hbar/ (mRb87 * r2)) * Rx
        
        projection += (fx+fy)
        vector_flow -= np.dstack((fx,fy))
        
        # and finally, mask the projection outside of the condensate:
        vector_flow[(X)**2 + (Y)**2 > R**2] = 0
        vector_flow = np.nan_to_num(vector_flow)
        return vector_flow    
    
v_x = np.zeros([N,N],'double')
v_y = np.zeros([N,N],'double')


if len(points_of_interest):
    vortex_coordinates = vortex_locations
    for i in range(len(vortex_coordinates)):
        modpsi *= vortex_core(*vortex_coordinates[i])
        vector_flow = single_vortex_flow(vortex_coordinates[i])
        v_x += vortex_signs[i] * vector_flow[:,:,0]
        v_y += vortex_signs[i] * vector_flow[:,:,1]

    ## renormalise psi
    psi_norm = np.sum(modpsi**2)*dx*dy/(n_atoms)
    modpsi =modpsi/np.sqrt(psi_norm)
    
    
    # =======================================================
    # Energy spectrum calculation
    # =======================================================
    
    # Density weighted velocity field
    u_x = modpsi * v_x
    u_y = modpsi * v_y
   
    # Fourier transform density-weighted velocity field
    u_k_x = np.fft.fft2(u_x)*dx*dy/(2*np.pi)
    u_k_y = np.fft.fft2(u_y)*dx*dy/(2*np.pi)

    # Incompressible kinetic energy density
    KE_inc_2d = 0.5 *mRb87* (np.abs(u_k_x)**2 + np.abs(u_k_y)**2)

    # Get rid of imaginary numerical noise
    KE_inc_2d = np.real(KE_inc_2d)

    # |k| bins
    # Number of bins in histogram
    n_bins = N/2
    # All radial bins beyond k_max fall partially outside of the grid
    k_vec  = np.linspace(0,np.sqrt(2)*k_max, n_bins+1)
    dk_vec = k_vec[1] - k_vec[0]
    
    # Multiply energy density by |k|
    KE_inc_times_k_2d = KE_inc_2d * K_mag

    # Radial average
    KE_inc_spec = radialaverage(KE_inc_times_k_2d, K_x, K_y, k_vec, dk_vec)
    KE_inc_spec =  np.nan_to_num(KE_inc_spec)

    rho = mRb87 / (psi_norm)
    rkp = (rho * Kappa**2 / (4*np.pi))

    KE_inc_spec /= rkp
    
    ## Now integrate the spectrum to get the total kinetic energy.
    
    KE_total = sum(KE_inc_spec) * dk_vec * 2 * np.pi
    run.save_result('KE_total',KE_total)
    run.save_result('KE_per_vortex', KE_total/len(vortex_locations))
    
    ### The following lines can be uncommented to check that the energy result is the same
    ### as that achieved by integrating over the velocity field.
    # KE_inc_2d_xy = 0.5 * mRb87 * (np.abs(u_x)**2 + np.abs(u_y)**2)
    # KE_total_v = np.sum(KE_inc_2d_xy) * dx * dy / rkp
    # print "Integrated spectrum energy: ",KE_total
    # print "Integrated velocity energy:", KE_total_v
    
    k_vec *= xi
    run.save_result_array('KE_inc_spec',KE_inc_spec)
    run.save_result_array('KE_per_vortex_spec',KE_inc_spec/len(vortex_locations))
    run.save_result_array('k_vec',k_vec)


    
