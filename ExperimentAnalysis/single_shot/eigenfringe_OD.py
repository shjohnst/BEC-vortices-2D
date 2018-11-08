###============================eigenfringe_OD.py============================###
###                                                                         ###
### Computes an optical density image based on a reference frame generated  ###
### from an eiginbasis of reference images.                                 ###
###                                                                         ###
### Copyright (C) 2018  Shaun Johnstone and Christopher Billington          ###
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
import h5py, lyse

def reconstruct(target, refs, mask):
    """Reconstruct a target image as a linear sum of reference images using
    weighted linear least squares with weights equal to 1 everywhere where
    mask is 1 and 0 where mask is 0. As such, the masked region is
    reconstructed based on the least squares solution in the unmasked
    region"""

    # Flatten images into the right shaped arrays
    n_refs, ny, nx = np.shape(refs)
    W = mask.flatten().astype(bool)
    a = target.flatten().astype(float)
    B = np.vstack(r.flatten().astype(float) for r in refs)

    # Solve the weighted linear least squares problem find the best reconstruction of a as a
    # linear sum of reference images, see https://en.wikipedia.org/wiki/Least_squares#Weighted_least_squares:
    #
    # We solve the following for x:
    #
    # B W B.T x = B W a
    #
    # B is a matrix with each row being one reference image, W is a diagonal matrix of weights
    # (only ones and zeros in our case and so equal to the mask), and a is the target image we are reconstructing.
    # The reconstructed image can then be computed as as:
    #
    # a_recon = x.T B
    #

    BW = B * W

    # The matrix of inner products between masked reference images. This is
    # the most expensive part of the algorithm. If you are wanting to use,
    # say, the most recent 100 frames or something, you may wish to store this
    # and update only one row and column of it when you replace the oldest
    # reference image with a newer one at position i (in both B and BW) like
    # this: P[:, i] = P[i, :] = np.dot(BW[i], B.T)
    
    P = np.dot(BW, B.T)

    # The RHS of the least squares equation:
    orig_projs = np.dot(BW, a.T)

    # Solve for x:
    x = np.linalg.solve(P, orig_projs)

    # Reconstruct the target image
    recon = np.dot(x.T, B)

    return recon.reshape(ny, nx)

run = lyse.Run(lyse.path)
dataset = lyse.data(lyse.path)    
    
try:
    lyse.routine_storage.refs
except:
    lyse.routine_storage.refs = None
    lyse.routine_storage.dark_frames = None


with h5py.File(lyse.path) as h5file:
    for orientation, images in h5file["images"].iteritems():
        if orientation == "bottom":
            for image_type, frameset in images.iteritems():
                if image_type.startswith("absorption"):
                    element = image_type.split("_")[1]
                    isat_name = "_".join([orientation, element, "saturation_intensity"])
                    if "Cameras" in h5file["globals"]:
                        if dataset.get(isat_name, False):
                            Isat = dataset.get(isat_name,5000)
                    atoms = float32(frameset["atoms"])
                    new_flat = float32(frameset["flat"])
                    new_dark = float32(frameset["dark"])
                    
                    
                    ## Add this shot's reference frame to our basis if needed
                    if lyse.routine_storage.refs is None:
                        lyse.routine_storage.refs = array([new_flat])
                        lyse.routine_storage.dark_frames = array([new_dark])
                    elif any([array_equal(lyse.routine_storage.refs[i,:],new_flat) for i in range(len(lyse.routine_storage.refs))]):
                        #nothing to do, as we're re-analysing a shot that is already in the basis
                        print "Existing shot"
                        pass
                    elif len(lyse.routine_storage.refs) == 200:
                        lyse.routine_storage.refs = roll(lyse.routine_storage.refs,-1, 0)
                        lyse.routine_storage.refs[-1] = new_flat
                        lyse.routine_storage.dark_frames = roll(lyse.routine_storage.dark_frames,-1,0)
                        lyse.routine_storage.dark_frames[-1] = new_dark                        
                    else:
                        lyse.routine_storage.refs = append(lyse.routine_storage.refs,[new_flat],0)
                        lyse.routine_storage.dark_frames = append(lyse.routine_storage.dark_frames, [new_dark],0)
                    print "Basis contains %s frames" %len(lyse.routine_storage.refs)
                    
                    mask = ones(atoms.shape,dtype=bool)
                    # mask regions roi0 (BEC location), roi1 & roi2 (Bragg-scattered atom locations)
                    for roi in ['roi0','roi1','roi2']:
                        x1,y1,x2,y2 = dataset[orientation, roi]
                        mask[y1:y2,x1:x2] = 0
                    lx, ly = atoms.shape
                    X, Y = np.ogrid[0:lx, 0:ly]
                    
                    ## also mask an elliptical region around the edges of the camera where no light is hitting
                    mask2 = ((Y - 245) / 240) ** 2 + ((X - 256)/270) ** 2 < 1
                    
                    mask = mask & mask2
                    imshow(mask*atoms)
                    title("masked atom frame")
                    reconstructed_flat = reconstruct(atoms,lyse.routine_storage.refs,mask)
                    average_dark = mean(lyse.routine_storage.dark_frames,0)
                    atoms -= average_dark
                    atoms *= mask2
                    reconstructed_flat -= average_dark
                    reconstructed_flat *= mask2
                    OD = (-log(atoms/reconstructed_flat) + (reconstructed_flat-atoms)/Isat)
                    OD[isnan(OD)] = 0.0
                    OD[isinf(OD)] = 0.0
                    figure()
                    imshow(OD,vmin=0, vmax = 2, cmap="gray", interpolation="none")
                    colorbar()
                    title("Optical density based on reconstructed flat frame")
                    figure()
                    imshow(atoms)
                    title("Raw atoms frame")
                    figure()
                    imshow(new_flat)
                    title("Latest flat field")
                    figure()
                    imshow(reconstructed_flat)
                    title("Reconstructed flat field")
                    figure()
                    imshow(new_dark)
                    colorbar()
                    title('Latest dark frame')
                    
                    
                    run.save_result_array("OD", OD)

                


