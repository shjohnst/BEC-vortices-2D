###==============vortex_signs_classification_and_statistics.py==============###
###                                                                         ###
### Calculates the best matching vortex signs for a given Bragg             ###
### differential signal and list of vortex locations.                       ###
### Uses the calculated signs to compute statistics about the vortex        ###
### distribution such as correlation functions and the dipole moment,       ###
### as well as classifying the vortices as belonging to clusters, dipole    ###
### pairs, or free vortices.                                                ###    
###                                                                         ###
### Copyright (C) 2018 Shaun Johnstone and Philip Starkey                   ###
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
import time
print 'script started at: ', time.time()
import os
from pylab import *
from matplotlib.patches import Circle
import lyse
from lyse import *
import itertools

import networkx as NX
from networkx.algorithms.components import number_connected_components

from analysislib.johnstone_vortices_2018.parameters import *


## This script can be run directly from lyse, however, it can also be run from another script.
## We use this to save time by splitting the vortex configuration space across multiple CPU cores to run in parallel
## if lyse.spinning_top == True we are running direct, otherwise if lyse.spinning_top == False the code is running via the multiprocessing code.
## On our computer, with 8 processors, the benefits of multiprocessing appear to outweigh the overheads of setting it up if there are 15 or more vortices.
n_vortices_single_process = 14

if not lyse.spinning_top:
    import sys
    import ast
    path = sys.argv[1]
    initial_vortex_config = ast.literal_eval(sys.argv[2])
    print("First two vortex states are set as: ")
    print(initial_vortex_config)
else:
    initial_vortex_config = None
    

run = Run(path)

def run_code(config_list = None, run = run):
    st = time.time()
    points_of_interest = []
    with h5py.File(run.h5_path,'r') as h5_file:
        if "/results/eigenfringe_OD" in h5_file:
            for key, val in h5_file["results"]["eigenfringe_OD"].iteritems():
                if key =="OD":
                    raw_OD = val[:]
        
        
        if "results/review_vortex_locations" in h5_file:
            for key, val in h5_file["results"]["review_vortex_locations"].iteritems():
                if key =="vortex_points":
                    points_of_interest = val[:]
        
    
    
    
    differential_OD = raw_OD[ROIs['roi2']['y1']:ROIs['roi2']['y2'],ROIs['roi2']['x1']:ROIs['roi2']['x2']] - raw_OD[ROIs['roi1']['y1']:ROIs['roi1']['y2'],ROIs['roi1']['x1']:ROIs['roi1']['x2']]
    

    print 'init time: ', time.time()-st
    st = time.time()
    bragg_signs = sign(differential_OD)



    def single_vortex_flow(vortex_loaction):
        '''Calculates the flow due to a single vortex and its image vortex using the point vortex model.
        The sign of the field returned assumes that the vortex charge is 1, (and hence the image has charge -1).
        This flow field can then be multiplied by the charge to obtain the field of a vortex with known charge.'''
        
        X,Y = meshgrid(range(wx),range(wy))
        px = vortex_loaction[0]
        py = vortex_loaction[1]
        Rx = pixel_size * (X -px )
        Ry = pixel_size * (Y - py)
        r2 = Rx**2 + Ry**2
        
        
        fx = (hbar/ (mRb87 * r2)) * -Ry
        fy = (hbar/ (mRb87 * r2)) * Rx
        
        projection = -(fx+fy)
        vector_flow = dstack((fx,fy))
        
        # now for the image vortex. First, get its location
        pox = px - condensate_centre[0]
        poy = py - condensate_centre[1]
        por = sqrt(pox**2 + poy**2)
        poa = arctan2(poy,pox)
        # the radius of the image vortex is then:
        ir = circ_radius**2 / por
        iox = ir * cos(poa)
        ioy = ir * sin(poa)
        # and finally, add the origin of the BEC again to get back into our base coordinate system
        ix = iox + condensate_centre[0]
        iy = ioy + condensate_centre[1]
        
        Rx = pixel_size * (X - ix)
        Ry = pixel_size * (Y - iy)
        r2 = Rx**2 + Ry**2
        
        fx = (hbar/ (mRb87 * r2)) * -Ry
        fy = (hbar/ (mRb87 * r2)) * Rx
        
        projection += (fx+fy)
        vector_flow -= dstack((fx,fy))
        
        # and finally, mask the projection outside of the condensate:
        projection[(X-condensate_centre[0])**2 + (Y-condensate_centre[1])**2 > circ_radius**2] = 0
        vector_flow[(X-condensate_centre[0])**2 + (Y-condensate_centre[1])**2 > circ_radius**2] = 0
        
        return vector_flow
        
    # Points corresponding to the vortex locations were stored in a format which is not particularly user friendly, and relative to the whole camera chip.
    # Here we will convert these into a more useful format.
    # This new list of points will be in the same coordinate system as the differential signal.
    if len(points_of_interest):
        vortex_locations = []
        image_vortex_locations = []
        
        vortex_mask = zeros((wy,wx))
        X,Y = meshgrid(range(wx),range(wy))
        
        # We are only going to compare the Bragg signal within 4 pixels of a vortex.
        vortex_roi = 4.0
        
        for point in points_of_interest:
            px = point[0][0] - ROIs['roi0']['x1']
            py = point[0][1] - ROIs['roi0']['y1']
            vortex_locations.append((px,py))
            vortex_mask[(X-px)**2 + (Y-py)**2 < vortex_roi**2] = 1
        
        n_vortices = len(vortex_locations)
        
        print 'There are %s vortices' %n_vortices
        
        ## This estimated run time will be dependant on your CPU!
        ## Note that for each extra vortex, the calculation time is doubled, so it is useful to know in advance how long to expect the code to take.
        estimated_duration = 2**n_vortices * 0.0013 / 120
        estimated_duration = 2**n_vortices * 3.82e-6
        
        print 'prep time: ', time.time()-st
        ###### Check if we should run the calculations on this experiment ######
        
        # Since this script scales like 2^n_vortices, you may want to limit the analysis to only low vortex numbers during testing, by changing the limits on the following line.
        # Here, we're limiting analysis to 2 to 13 vortices for the single-process version of the code (see comment at top of script)
        if n_vortices <= n_vortices_single_process and not lyse.spinning_top:
                print "Skipping shot as it has %i vortices: should be analysied with single process"%n_vortices
                sys.exit(0)
        elif (n_vortices < 2 or n_vortices > n_vortices_single_process) and lyse.spinning_top and config_list is None:
                print "Skipping shot as it has %i vortices"%n_vortices
        ###### This part only runs if there are a suitable number of vortices (i.e. gets past the check above) ######
        else:
            print "Iterating over %s possible vortex sign configurations, estimated duration is %s minutes." %(2**n_vortices,estimated_duration)
        
            # firstly, lets make a big array of "base" flow projections, for each vortex.
            vortex_vector_fields = map(single_vortex_flow, vortex_locations)
            
            # only want the projections along [-1,-1] for most of this:
            projected_vortex_fields = -sum(vortex_vector_fields, axis = 3)
            
                        
            ###### Loop over all possible combinations of the vortex fields, and keep the one that matches the Bragg signal the best ######
            goodness = 0
            start_time = time.time()
            
            if config_list is None:
                if not lyse.spinning_top:
                    n_vortices -= len(initial_vortex_config)
                    config_list = itertools.product([-1, 1], repeat=n_vortices)
                else:
                    config_list = itertools.product([-1, 1], repeat=n_vortices)
                optimise_for_duplicates = 2**n_vortices
                
                
            else:
                optimise_for_duplicates = None
                
                
            section_a = 0.0
            section_b = 0.0
            section_c = 0.0
            section_d = 0.0
            section_e = 0.0
            
            current_iteration = 1
            goodness_list_for_plot = []
            for current_sign_configuration in config_list:
                st = time.time()
                
                current_sign_configuration = list(current_sign_configuration)
                if not lyse.spinning_top:
                    current_sign_configuration =  current_sign_configuration + initial_vortex_config
                
                section_a+=time.time()-st
                st = time.time()
                
                new_full_projected_field = einsum('i,ijk->jk',current_sign_configuration,projected_vortex_fields)
                
                section_b+=time.time()-st
                st = time.time()
                
                sign_field = vortex_mask * sign(new_full_projected_field)
                
                section_c+=time.time()-st
                st = time.time()
                
                new_goodness1 = (sign_field == bragg_signs).sum()
                goodness_list_for_plot.append(new_goodness1)
                
                # Since the vortex configuration space is symmetric, we can check the match of the configuration that has all signs opposite
                # compared with the current configuration by just multiplying it by -1, rather than recomputing it from scratch.
                # This significantly reduces the time taken for the full calculation.
                if optimise_for_duplicates is not None:
                    new_goodness2 = (-sign_field == bragg_signs).sum()
                    goodness_list_for_plot.append(new_goodness2)
                section_d+=time.time()-st
                st = time.time()
                
                
                ## Check to see if the current configuration (or the opposite of it) is a better match than the current best configuration
                ## If so, save the current values as our best configuration
                if new_goodness1 > goodness:
                    goodness = new_goodness1
                    full_projected_field = new_full_projected_field
                    vortex_signs = current_sign_configuration
                    full_vector_field = einsum('i,ijkl->jkl',vortex_signs,vortex_vector_fields)
                if optimise_for_duplicates is not None: 
                    if new_goodness2 > goodness:
                        goodness = new_goodness2
                        full_projected_field = -new_full_projected_field
                        vortex_signs = [-el for el in current_sign_configuration]
                        full_vector_field = einsum('i,ijkl->jkl',vortex_signs,vortex_vector_fields)
                        
                section_e+=time.time()-st
                        
                # Since the vortex configuration space is symmetric, we can break out after half the iterations
                if optimise_for_duplicates is not None and current_iteration >= optimise_for_duplicates/2:
                    print "Optimised for duplicates"
                    break
                current_iteration += 1
            
            ## Optional: Print out timing for code optimisation ##
            # print '---'
            # print "section times"
            # print section_a
            # print section_b
            # print section_c
            # print section_d
            # print section_e
            # print '---'
            end_time = time.time()
            print "Took %s seconds to loop over all possibilities" %(end_time-start_time)
            
            print "Best vortex configuration:"
            ## Note: This must be the last print statement before the sys.exit() when run from the multiprocessing script
            print vortex_signs
            
            ## If we were running this code from the multiprocessing script, it's now time to exit.
            if not lyse.spinning_top:
                sys.exit(0)
            
            

            ## Optional: plot and save the data about how good each possible vortex configuration matches the Bragg signal (used to generate Fig. S8)
            ## Turn this on by setting plot_goodness = True
            plot_goodness = False
            if plot_goodness:
                font = {'family' : 'sans-serif',
                        'weight' : 'normal',
                        'size'   : 7}

                matplotlib.rc('font', **font)
                # fig = figure(figsize=(4.75,3.5))
                fig = figure("Vortex configuration match")
                
                goodness_list_for_plot.sort()
                
                ## Save the data for the "plot of goodness" to a csv file for later use
                
                savetxt("%s_vortex_configuration_match.csv"%os.path.splitext(os.path.basename(run.h5_path))[0], goodness_list_for_plot, delimiter=",")
                
                scatter(range(len(goodness_list_for_plot)),goodness_list_for_plot,linewidths =0)
                ylabel("Number of matching pixels")
                xlabel("Configuration")
                axhline(y=sum(vortex_mask), c= 'k', ls='--',lw = 1)
                axhline(y=sum(vortex_mask)/2.0, c= 'k', ls='-.',lw = 1)
                axhline(y=mean(goodness_list_for_plot), c= 'k', ls=':',lw = 1)
                
                xlim(-0.05*len(goodness_list_for_plot),1.05*len(goodness_list_for_plot))
                
                print "Number of possible configurations = %s"%len(goodness_list_for_plot)
                print "Number of pixels being compared: %s (within %s pixels of each of %s vortices)"%(sum(vortex_mask),vortex_roi,n_vortices)
                print "Best match: %s pixels" %max(goodness_list_for_plot)
                print "Average match: %s pixels (%s percent of all pixels under consideration)"%(mean(goodness_list_for_plot),100*mean(goodness_list_for_plot)/sum(vortex_mask))

                
                
            ### For the remaining calculations, we're going to recalculate the vortex locations in units
            ### of the trap radius, centred on the trap:
            vortex_locations = []
            for point in points_of_interest:
                # We want the vortex locations relative to the centre of the cloud, in units of the cloud radius
                px = (point[0][0] - ROIs['roi0']['x1'] - wx/2)/circ_radius
                py = (point[0][1] - ROIs['roi0']['y1'] - wy/2)/circ_radius
                vortex_locations.append([px,py])
            vortex_locations = array(vortex_locations)
            
            ###### Cluster/Dipole Classification ######
            ### For details, please refer to "Einstein-Bose condensation of Onsager vortices"
            ### R. N. Valani, A. J. Groszek, and T. P. Simula, New Journal of Physics 20, 53038 (2018)
            # First, get nearest opposites
            nearest_opposites = []
            for i in range(n_vortices):
                d = 9000
                nv = None
                for j in range(n_vortices):
                    # Only looking for cases where this vortex isn't in a dipole, and the other vortex has the opposite sign
                    if vortex_signs[i] == -vortex_signs[j]:
                        dn = norm(vortex_locations[i]-vortex_locations[j])
                        if dn < d:
                            d = dn
                            nv = j
                # in this case, we care about the distance, not who it is
                nearest_opposites.append([d,nv])
                    
            # Second, get list of vortices closer than the nearest opposite
            cluster_candidates = []
            nearest_vortex_distance = []
            for i in range(n_vortices):
                local_cluster_candidates = []
                current_closest_found = 9000
                for j in range(n_vortices):
                    if any(vortex_locations[i] != vortex_locations[j]): # not the same vortex
                        dn = norm(vortex_locations[i]-vortex_locations[j])
                        if dn < nearest_opposites[i][0]:
                            local_cluster_candidates.append(j)
                        
                        if dn < current_closest_found:
                            current_closest_found = dn
                cluster_candidates.append(local_cluster_candidates)
                nearest_vortex_distance.append(current_closest_found)
            
            nearest_vortex_distance = array(nearest_vortex_distance)
            
            # Now, go through and check: 1) If no cluster candidates, is the nearest opposite mutual, with no cluster candidates? if So, it's a dipole.
            # 2) If there are cluster candidates, are they mutual? If so, it's a cluster.
            
            dipoles = []
            clusters = []
            for i in range(n_vortices):
                if cluster_candidates[i] == []:
                    # there are no cluster candidates, so we need to see if this is part of a dipole
                    dipole_candidate = nearest_opposites[i][1]
                    if cluster_candidates[dipole_candidate] == [] and nearest_opposites[dipole_candidate][1] == i:
                        # We've got a dipole!
                        dipoles.append([i,dipole_candidate])
                        # Now mark i's nearest_opposite as None to prevent it being added a second time when we get up to testing dipole_candidate
                        nearest_opposites[i][1] = None
                else:
                    # there are cluster candidates, are they mutual?
                    for j in cluster_candidates[i]:
                        if i in cluster_candidates[j]:
                            clusters.append((i,j,{'weight': norm(vortex_locations[i]-vortex_locations[j])}))
            
            ## Now that we know which vortices are in dipoles and clusters, we can also calculate the nearest neighbour distance for each class of vortex
            dipole_distances = []
            for pair in dipoles:
                # Only need to worry about adding one of the 2, as they are mutual nearest anyway
                dipole_distances.append(nearest_vortex_distance[pair[0]])
              
            # Make a graph out of the edges between the clusters:
            cluster_graph = NX.Graph()
            cluster_graph.add_edges_from(clusters)
                
            number_of_clusters = number_connected_components(cluster_graph)
            n_clustered = NX.number_of_nodes(cluster_graph)
            if n_clustered:
                mean_vortices_per_cluster = n_clustered/number_of_clusters
            else:
                mean_vortices_per_cluster = 0
            
            run.save_result('mean_vortices_per_cluster',mean_vortices_per_cluster)
            run.save_result('number_of_clusters',number_of_clusters)
            
            ## Calculate the cluster edges in case we want to use them later
            T = NX.minimum_spanning_tree(cluster_graph)
            cluster_plot_lines = NX.edges(T)
            

            
            # We may as well also calculate the nearest neighbour distance for each clustered vortex
            cluster_nearests = []
            for graph in NX.connected_component_subgraphs(cluster_graph):
                nodes =  NX.nodes(graph)
                node_coordinates = []
                for node in nodes:
                    node_coordinates.append(vortex_locations[node])
                    cluster_nearests.append(nearest_vortex_distance[node])
                
            
            ## So, now we have a list of dipole sizes and cluster sizes, we can save these results
            
            mean_dipole_length = mean(dipole_distances)
            mean_cluster_spacing = mean(cluster_nearests)
            mean_nearest_vortex = mean(nearest_vortex_distance)
            mean_vortex_spacing =  2.0/sqrt(n_vortices)
            
            run.save_result('mean_dipole_length',mean_dipole_length)
            run.save_result('mean_cluster_spacing',mean_cluster_spacing)
            run.save_result('mean_nearest_vortex',mean_nearest_vortex)
            run.save_result('mean_vortex_spacing',mean_vortex_spacing)
            
            
            run.save_result_array('dipoles',dipoles)
            run.save_result_array('cluster_edges',cluster_plot_lines)
            
            
            ###### Count each classification group ######
            # We have already counted the clusters in the code above, using a Graph
            # Dipole counting is easy, since we made sure that we didn't double count:
            n_dipole = 2 * len(dipoles)
            # what is left must be free vortices
            n_free = n_vortices - n_clustered - n_dipole

            ###### Calculate net polarisation ######
            polarisation = sum(vortex_signs)
            
            
            ## Now that everything is done, we can print and save the results
            
            print "Total number of vortices: %s" %n_vortices
            print "Number of vortices in dipoles: %s (%s %%)" %(n_dipole , 100.*n_dipole/n_vortices)
            print "Number of vortices in clusters: %s (%s %%)" %(n_clustered, 100.*n_clustered/n_vortices)
            print "Number of free vortices: %s (%s %%)" %(n_free, 100. * n_free/n_vortices)
            print "Polarisation: %s" %polarisation
            
            run.save_result('N_v',n_vortices)
            run.save_result('N_d',n_dipole)
            run.save_result('N_c',n_clustered)
            run.save_result('N_f',n_free)
            run.save_result('p',polarisation)
            
            run.save_result_array('calculated_flow',full_projected_field)
            run.save_result_array('calculated_flow_vector',full_vector_field)
            run.save_result_array('vortex_signs',vortex_signs)
            
            
            
            ### Now generate some figures to show the calculated Bragg signal for the best match, and the classified vortices that have been identified
            
            ## Plot vortex locations over OD ##
            figure("roi0")
            imshow(raw_OD, vmin = 0, vmax=2, cmap="gray",  interpolation="none")
            colorbar()
            title("roi0")
            xlim(ROIs['roi0']['x1'],ROIs['roi0']['x2'])
            ylim(ROIs['roi0']['y2'],ROIs['roi0']['y1'])
            circ = Circle((ROIs['roi0']['x1']+wx/2.0,ROIs['roi0']['y1']+wy/2.0), circ_radius, fill=False, ls = "--", ec="Red")
            ax = gca()
            ax.add_patch(circ)
            if len(points_of_interest):
                for point in points_of_interest:
                    plot(*point[0], color='b', marker='o', label = "POI", picker = True)
                
            ## Plot Bragg signal with vortex locations ##
            figure("Differential")
            imshow(differential_OD, vmin = -1, vmax = 1,cmap="RdBu",  interpolation="none")
            colorbar()
            title("Differential")
            circ = Circle((wx/2.0,wy/2.0), circ_radius, fill=False, ls = "--")
            ax = gca()
            ax.add_patch(circ)
            if len(points_of_interest):
                for point in points_of_interest:
                    px = point[0][0] - ROIs['roi0']['x1']
                    py = point[0][1] - ROIs['roi0']['y1']
                    plot(px,py, color='k', marker='x')
            
            xlim(0,wx)
            ylim(wy,0)
            
            ## Plot calculated flow field ##
            figure('Calculated flow field')
            imshow(full_projected_field,cmap="RdBu",vmin = -0.001, vmax = 0.001)
            title('Projected flow field of best match candidate')
            colorbar()
            
            ## Plot the classified vortices, with lines joining the clusters and dipoles
            figure("Clasified vortices")
            # first, plot all the vortices
            vortex_locations_array = array(vortex_locations)
            pos_neg_cmap = matplotlib.colors.ListedColormap(['g','b'], name='from_list', N=2)
            scatter(vortex_locations_array[:,0],vortex_locations_array[:,1],c=vortex_signs,marker = 'o', lw = 0, cmap = pos_neg_cmap)
            
            # then, add lines between dipoles
            for pair in dipoles:
                plot([vortex_locations_array[pair[0],0],vortex_locations_array[pair[1],0]],[vortex_locations_array[pair[0],1],vortex_locations_array[pair[1],1]],c = 'r')
            
                
            # and finally, join the dots for the clusters
            for edge in cluster_plot_lines:
                if vortex_signs[edge[0]] == 1:
                    col = 'b'
                else:
                    col = 'g'
                plot([vortex_locations_array[edge[0],0],vortex_locations_array[edge[1],0]],[vortex_locations_array[edge[0],1],vortex_locations_array[edge[1],1]],c = col)
            
            
            circ = Circle((0,0), 1, fill=False, ls = "--")        
            ax = gca()
            ax.add_patch(circ)
            xlim(-wx/(2*circ_radius),wx/(2*circ_radius))
            ylim(wy/(2*circ_radius),-wy/(2*circ_radius))
            colorbar()
            
            
            
            ##### Net dipole moment calculation #####
            summed_signed_locations = array([0.0,0.0])
            for i in range(len(vortex_signs)):
                summed_signed_locations += vortex_signs[i] * vortex_locations[i]
            
            dipole_moment = norm(summed_signed_locations)/n_vortices
            run.save_result('D',dipole_moment)
            
            
            ###### Now get first and second order correlation functions ######
            
            ## Need 1st & second nearest neighbours
            
            ## C_1 = 1/(Nv) SUM_i->Nv((c_i1)), where c_i1 = 1 for same sign nearest neighbour, or -1 for opposite sign
            ## C_2 = 1/(2*Nv) SUM_i->Nv(SUM_j=1,2(c_ij)), where c_ij = 1 if ith vortex is same sign as jth neares neighbour, else -1
            

            C_1 = 0
            C_2 = 0
            for i in range(n_vortices):
                distances = []
                for j in range(n_vortices):
                    vx = vortex_locations[i][0]
                    vy = vortex_locations[i][1]
                    
                    nx = vortex_locations[j][0]
                    ny = vortex_locations[j][1]
                    
                    dn = sqrt((vx-nx)**2 + (vy-ny)**2)
                    distances.append(dn)
                
                # Now, find the index of the closest and second closest. Note that there will also be a 0 distance vortex - itself!
                sorted_indices = argsort(distances)
                # nearest_neighbours.append(sorted_indices[1])
                # second_nearest.append(sorted_indices[2])
                
                ## First add the nearest neighbour term
                C_1 += vortex_signs[i] * vortex_signs[sorted_indices[1]]
                C_2 += vortex_signs[i] * vortex_signs[sorted_indices[1]]
                ## Then add the second-nearest term to C_2 -- only works for 3 or more vortices!
                if len(sorted_indices)>2:
                    C_2 += vortex_signs[i] * vortex_signs[sorted_indices[2]]
            ## Normalise
            C_1 = C_1/n_vortices
            C_2 = C_2/(2*n_vortices)
            
            run.save_result('C_1',C_1)
            run.save_result('C_2',C_2)

if __name__ == "__main__":
    run_code()