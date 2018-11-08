###=================generate_time_series_plots_per_grid.py==================###
###                                                                         ###
### Generates Figures 4, S3-S6 for "Evolution of large-scale flow from      ###
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
import xlsxwriter
from scipy.stats import sem
from labscript_utils.labconfig import LabConfig
import os
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerLine2D
from mpl_toolkits.axes_grid.inset_locator import InsetPosition
from vortex_thermometer import thermometer_cdN
from analysislib.johnstone_vortices_2018.parameters import *

## Check where to save output:
exp_config = LabConfig()
results_path = exp_config.get('paths', 'analysis_output_folder')

sequence_list = [   
                    '20171004T093812',  #4.2
                    '20171004T121555',  #6.0
                    '20171005T090729',  #7.9
                    '20171004T144649',  #9.7
                    '20171006T101525'   #11.5
                    ]

## Get the dataframe from lyse
full_df = lyse.data(timeout=10)

## Setup a spreadsheet for saving data to
workbook = xlsxwriter.Workbook(os.path.join(results_path,'time_data.xlsx'))
bold = workbook.add_format({'bold': True})
bold_center = workbook.add_format({'bold': True,'align':'center'})
boldshaded = workbook.add_format({'bold':True, 'bg_color': '#CCCCCC'})
shaded = workbook.add_format({'bg_color': '#CCCCCC','align':'center','left':1,'right':1})
leftborder = workbook.add_format({'left':1})
rightborder = workbook.add_format({'right':1})
boldleftborder = workbook.add_format({'bold': True,'left':1})
boldrightborder = workbook.add_format({'bold': True,'right':1})
## Add the headers
data_table_headings = ['Hold time',
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
                            'sem(E_k)',
                            'N_v',
                            'sem(N_v)',
                            'l_v',
                            'sem(l_v)']


def plot_time_series(specify_sequence):
    '''Generates a figure with graphs of data vs time for a specified sequence of experiments,
    corresponding to a fixed grid width.'''
    
    print "Working on sequence", specify_sequence
        
    ## Get a dataframe containing only the data for the current grid
    sequences = [ os.path.split(path)[1][0:15] for path in full_df['filepath'] ]
    df = full_df[ [(sequence == specify_sequence) for sequence in sequences] ]
    ## Throw out rejected shots
    df = df[df['review_vortex_locations','Good']==True]
    
    grid_radius = df['vortex_beam_radius'][0] * dmd_pixel_size
    grid_radius_micron = grid_radius*1e6
    ## Set up worksheet for data saving
    worksheet = workbook.add_worksheet('%s micron grid'%grid_radius_micron)
    worksheet.write(0,0,'Time series data for %s micron grid'%grid_radius_micron,bold)
    
    for head_i, heading in enumerate(data_table_headings):
        worksheet.write(1,head_i,heading,bold)
    
    worksheet_spectra = workbook.add_worksheet('%s micron grid Spectra'%grid_radius_micron)
    worksheet_spectra.write(0,0,'Kinetic energy spectra for %s micron grid'%grid_radius_micron,bold)
    worksheet_spectra.write(1, 0, 'Hold time:', boldshaded)
    worksheet_spectra.write(2,0,'k vector',bold)
    
    ### Set up the figure windows and axes for this grid ###
    fig = figure('dynamics %s'%specify_sequence,figsize=(4.75,2.8))
    
    ## We use gridspec to arrange the subfigures
    ## there are two with the same splitting, so that we can control padding better for the larger log-scale graph
    gs = GridSpec(5,3,width_ratios = [2,1,1])
    gs_log = GridSpec(5,3,width_ratios = [2,1,1])

    ## Create an axis for each graph
    class_ax = fig.add_subplot(gs[0, 0])
    corr_ax = fig.add_subplot(gs[1, 0])
    dm_ax = fig.add_subplot(gs[2, 0])
    temp_ax = fig.add_subplot(gs[3, 0])
    e_ax = fig.add_subplot(gs[4, 0])
    spec_ax = fig.add_subplot(gs_log[:3, 1:3])
    
    logN_ax = fig.add_subplot(gs_log[3:, 1])
    logL_ax = fig.add_subplot(gs_log[3:, 2])

    ## Set up an inset axis to place the colorbar on the energy spectrum plot in a nice position
    colorbar_location = [0.73,0.8,0.22,0.03]
    col_ax = axes([0,0,1,1],frameon = False)
    col_ip = InsetPosition(spec_ax,colorbar_location)
    col_ax.set_axes_locator(col_ip)

    ## Calculations for where to place the inset figure on the energy spectrum such that the x-axis lines up
    espec_x_min = 2e-2
    espec_x_max = 5
    espec_inset_x_min = 3.5e-2
    espec_inset_x_max = 1.25
    figw = log(espec_x_max)-log(espec_x_min)
    insetoffset = (log(espec_inset_x_min) - log(espec_x_min))/figw
    insetw = (log(espec_inset_x_max)-log(espec_inset_x_min))/figw
    delta_e_location = [insetoffset,.05,insetw,0.35]
    delta_e_ax = axes([0,0,1,1])
    delta_e_ip = InsetPosition(spec_ax,delta_e_location)
    delta_e_ax.set_axes_locator(delta_e_ip)
    
    ## An iterator to use to cycle through colours on the energy spectra
    e_spec_color=iter(cm.inferno_r(np.linspace(0.1,1,10)))
    
    
    ### Create empty lists in which we can later store the average value and the uncertainty of each measurement at each hold time
    
    ## Classified vortex numbers
    avg_N_c = []
    avg_u_N_c = []
    
    avg_N_d = []
    avg_u_N_d = []
    
    avg_N_f = []
    avg_u_N_f = []
    
    avg_N_v = []
    avg_u_N_v = []

    ## Dipole moment
    avg_d = []
    avg_u_d = []
    
    ## Total incompressible kinetic energy per vortex
    avg_en = []
    avg_u_en = []
    
    ## Correlation function
    avg_c1 = []
    avg_u_c1 = []
    
    ## Mean nearest-neighbor vortex distance
    avg_mean_nearest_vortex = []
    u_avg_mean_nearest_vortex = []
    
    ## and of course, we need a time vector
    t_vals = []
    
    ## Now, we group the dataframe by hold time, and for each time, add the time to t_vals,
    ## calculate the mean and sem of each measurement, and add them to their appropriate lists:
    spreadsheet_row = 0
    for hold_time, group in df.groupby('vortex_spoon_wait_time'):
        t_vals.append(hold_time)
        
        N_v = group[('vortex_signs_classification_and_statistics','N_v')]
        avg_N_v.append(mean(N_v))
        avg_u_N_v.append(sem(N_v))
        
        ## We want the fractional populations, so before averaging,
        ## each absolute number of clusters has to be divided by the
        ## total number of vortices in that shot.
        N_c = group[('vortex_signs_classification_and_statistics','N_c')]/N_v
        avg_N_c.append(mean(N_c))
        avg_u_N_c.append(sem(N_c))   
        
        N_d = group[('vortex_signs_classification_and_statistics','N_d')]/N_v
        avg_N_d.append(mean(N_d))
        avg_u_N_d.append(sem(N_d))
        
        N_f = group[('vortex_signs_classification_and_statistics','N_f')]/N_v
        avg_N_f.append(mean(N_f))
        avg_u_N_f.append(sem(N_f))
        
        
        D = group[('vortex_signs_classification_and_statistics','D')]
        avg_d.append(mean(D))
        avg_u_d.append(sem(D)) 
        
                
        eN = group[('energy_spectra','KE_per_vortex')]
        avg_en.append(mean(eN))
        avg_u_en.append(sem(eN))
        
                
        C_1 = group[('vortex_signs_classification_and_statistics','C_1')]
        avg_c1.append(mean(C_1))
        avg_u_c1.append(sem(C_1))
        

        mean_nearest_vortex = group[('vortex_signs_classification_and_statistics','mean_nearest_vortex')]
        avg_mean_nearest_vortex.append(mean(mean_nearest_vortex))
        u_avg_mean_nearest_vortex.append(sem(mean_nearest_vortex))
        
        
        
        ### Kinetic energy spectra ###
        
        ## We make a list containing arrays with the energy spectra of each shot, for this current hold time
        e_specs = []
        for path in group['filepath']:
            run = lyse.Run(path,no_write=True)
            e_spec = run.get_result_array('energy_spectra','KE_per_vortex_spec')
            k_vec = run.get_result_array('energy_spectra','k_vec')
            e_specs.append(e_spec)
            
        ## Now take that list of arrays and make it an array itself, so that we can take its average and sem across each shot at each k-vector
        e_specs = array(e_specs)
        av_e_spec = mean(e_specs,axis = 0)
        u_e_spec = sem(e_specs,axis = 0)
        
        ## We're going to plot this in the loop, so that we don't have to worry about storing these values for each time point independantly
        sca(spec_ax)
        ## Get the colour to use in this loop from the iterator
        col_now = next(e_spec_color)
        loglog(k_vec,av_e_spec, c = col_now,  label = r"$E(k)$, t = %s"%hold_time)
        
        ## We want to plot the difference from t=0.5 in the inset, so on the first loop, we will get the t=0.5 data,
        ## then subsequent loops will subtract this from the current time's data to get the difference.
        if hold_time == 0.5:
            initial_spec = av_e_spec
            delta_e_ax.plot(k_vec,zeros(shape(initial_spec)), c = col_now,  label = r"$E(k)$, t = %s"%hold_time)
        else:
            delta_spec = av_e_spec - initial_spec
            delta_e_ax.plot(k_vec,delta_spec*1e7, c = col_now,  label = r"$E(k)$, t = %s"%hold_time)
        
        
        ## Add energy spectra data for this time point to the spreadsheet:
        for k_vec_i in range(len(k_vec)):
            if hold_time == 0.5:
                worksheet_spectra.write(k_vec_i+3,0,k_vec[k_vec_i])
            worksheet_spectra.write(k_vec_i+3, 2*spreadsheet_row+1,av_e_spec[k_vec_i],leftborder)
            worksheet_spectra.write(k_vec_i+3, 2*spreadsheet_row+2,u_e_spec[k_vec_i],rightborder)
            
        worksheet_spectra.merge_range(1, 2*spreadsheet_row+1, 1, 2*spreadsheet_row+2,hold_time,shaded)
        worksheet_spectra.write(2,2*spreadsheet_row+1,'E(k)',boldleftborder)
        worksheet_spectra.write(2,2*spreadsheet_row+2,'sem(E(k))',boldrightborder)
        spreadsheet_row += 1
    ##### END OF LOOP OVER TIME #####
    
    ### Now, calculate the temperature at each hold time, based on the fraction of clusters and dipoles, and the total vortex number.
    beta_measured = thermometer_cdN(array(avg_N_c), array(avg_N_d), array(avg_N_v))
    beta_upper_error = thermometer_cdN(array(avg_N_c)+array(avg_u_N_c), array(avg_N_d)-array(avg_u_N_d), array(avg_N_v))
    beta_lower_error = thermometer_cdN(array(avg_N_c)-array(avg_u_N_c), array(avg_N_d)+array(avg_u_N_d), array(avg_N_v))

    ### Now we have everything that we want to plot! ###

    ## First save it to the spreadsheet
    for time_i, hold_time in enumerate(t_vals):
        worksheet.write(time_i+2,0,hold_time)
        
        worksheet.write(time_i+2,1,avg_N_c[time_i])
        worksheet.write(time_i+2,2,avg_u_N_c[time_i])
        
        worksheet.write(time_i+2,3,avg_N_d[time_i])
        worksheet.write(time_i+2,4,avg_u_N_d[time_i])
        
        worksheet.write(time_i+2,5,avg_N_f[time_i])
        worksheet.write(time_i+2,6,avg_u_N_f[time_i])
        
        worksheet.write(time_i+2,7,avg_c1[time_i])
        worksheet.write(time_i+2,8,avg_u_c1[time_i])
        
        worksheet.write(time_i+2,9,avg_d[time_i])
        worksheet.write(time_i+2,10,avg_u_d[time_i])
        
        worksheet.write(time_i+2,11,beta_measured[time_i])
        worksheet.write(time_i+2,12,beta_measured[time_i] - beta_upper_error[time_i])
        worksheet.write(time_i+2,13,beta_lower_error[time_i] - beta_measured[time_i])
        
        worksheet.write(time_i+2,14,avg_en[time_i])
        worksheet.write(time_i+2,15,avg_u_en[time_i])
        
        worksheet.write(time_i+2,16,avg_N_v[time_i])
        worksheet.write(time_i+2,17,avg_u_N_v[time_i])
        
        worksheet.write(time_i+2,18,avg_mean_nearest_vortex[time_i])
        worksheet.write(time_i+2,19,u_avg_mean_nearest_vortex[time_i])
        
        
    
    ## Temperature graph
    sca(temp_ax)
    temp_ax.errorbar(t_vals, beta_measured, yerr = array([beta_measured - beta_lower_error, beta_upper_error - beta_measured]), c='k',mfc = 'k',fmt="o",  ms=5, mew = 0, capsize = 1, capthick =1, ls = '-')
    temp_ax.axhline(y=0, c= 'k', ls=':',lw = 0.5)
    temp_ax.invert_yaxis()
    xlim(0,5.5)
    ylabel(r"$\beta$")


    
    ## Vortex number graph
    sca(logN_ax)
    logN_ax.errorbar(t_vals, avg_N_v, yerr=avg_u_N_v,c="k", fmt="o", ms=5,mew = 0, capsize = 1, capthick =1)
    logN_ax.set_xscale('log')
    logN_ax.set_yscale('log')
    ylabel("$N_v$")
    xlabel('Hold time (s)')
    ### Generate guide-to-eye lines to plot on the vortex number panel
    ### The values of these lines are chosen for each grid such that the lines are close to the data
    t_for_lines = linspace(0.4,6,12)
    ## 6.0
    if specify_sequence == '20171004T121555':
        ## 3.5 body
        tn25 = 15*t_for_lines**(-2./5)
        # 2 body
        tn1 = 20*t_for_lines**(-1.)
    ## 11.5
    elif specify_sequence == '20171006T101525':
        ## 3.5 body
        tn25 = 13*t_for_lines**(-2./5)
        # 2 body
        tn1 = 30*t_for_lines**(-1.)
    ##4.2
    elif specify_sequence == '20171004T093812':
        ## 3.5 body
        tn25 = 11*t_for_lines**(-2./5)
        # 2 body
        tn1 = 15*t_for_lines**(-1.)
    ## 9.7
    elif specify_sequence == '20171004T144649':
        ## 3.5 body
        tn25 = 13*t_for_lines**(-2./5)
        # 2 body
        tn1 = 25*t_for_lines**(-1.)
    ## 7.9
    else:
        ## 3.5 body
        tn25 = 14*t_for_lines**(-2./5)
        # 2 body
        tn1 = 25*t_for_lines**(-1.)
    
    loglog(t_for_lines, tn25, c='k', lw = 0.5, ls = '-',zorder = -1)
    loglog(t_for_lines, tn1, c='k', lw = 0.5, ls = ':',zorder = -1)

    xlim(0.4,6)
    ylim(3,25)



    ## Mean distance graph
    sca(logL_ax)
    logL_ax.errorbar(t_vals, avg_mean_nearest_vortex, yerr=u_avg_mean_nearest_vortex,c="k", fmt="o", ms=5,mew = 0, capsize = 1, capthick =1)
    logL_ax.set_xscale('log')
    logL_ax.set_yscale('log')
    ylabel("$l_v/R_\perp$")
    xlabel('Hold time (s)')
    ### Generate guide-to-eye lines to plot on the vortex spacing panel
    ### The values of these lines are chosen for each grid such that the lines are close to the data
    ## 6.0
    if specify_sequence == '20171004T121555':
        t12 = 0.22 * t_for_lines**(0.5)
        t15 = 0.23 * t_for_lines**(1./5)
    ## 11.5
    elif specify_sequence == '20171006T101525':
        t12 = 0.22 * t_for_lines**(0.5)
        t15 = 0.26 * t_for_lines**(1./5)
    ## 4.2
    elif specify_sequence == '20171004T093812':
        t12 = 0.23 * t_for_lines**(0.5)
        t15 = 0.27 * t_for_lines**(1./5)
    ## 9.7
    elif specify_sequence == '20171004T144649':
        t12 = 0.22 * t_for_lines**(0.5)
        t15 = 0.26 * t_for_lines**(1./5)
    ##7.9
    else:
        t12 = 0.2 * t_for_lines**(0.5)
        t15 = 0.25 * t_for_lines**(1./5)

    loglog(t_for_lines, t12, c='k', lw = 0.5, ls = ':',label = r'$t^{1/2}$',zorder = -1)
    loglog(t_for_lines, t15, c='k', lw = 0.5, ls='-', label = r'$t^{1/5}$',zorder = -1)

    xlim(0.4,6)
    ylim(0.18,0.55)

    
    ## Classified populations graph
    sca(class_ax)
    class_ax.errorbar(t_vals, avg_N_f, yerr=avg_u_N_f, c=colour_free, mec = colour_free, mfc = colour_free, fmt="o", ms=5, mew = 0, capsize = 1, capthick =1, ls = '-', label = "Free vortices")
    class_ax.errorbar(t_vals, avg_N_d, yerr=avg_u_N_d, c=colour_dipole, mec = colour_dipole, mfc = colour_dipole, fmt="o", ms=5, mew = 0, capsize = 1, capthick =1, ls = '-', label = "Dipole vortices")
    class_ax.errorbar(t_vals, avg_N_c, yerr=avg_u_N_c, c=colour_cluster, mec = colour_cluster, mfc = colour_cluster, fmt="o", ms=5, mew = 0, capsize = 1, capthick =1, ls = '-', label = "Clustered vortices")

    ylabel("$p_i$")
    xlim(0,5.5)
    
    ## Add a legend - this is done manually so that we can squash it down to the smallest possible size
    leg_cluster = mlines.Line2D([], [], color=colour_cluster, marker='o',
                          ms=5, mew = 0, ls = '-', label='$p_c$')
    leg_dipole = mlines.Line2D([], [], color=colour_dipole, marker='o',
                              ms=5, mew = 0, ls = '-', label='$p_d$')
    leg_free = mlines.Line2D([], [], color=colour_free, marker='o',
                              ms=5, mew = 0, ls = '-', label='$p_f$')

    legend( handles=[leg_cluster,leg_dipole,leg_free],
            bbox_to_anchor = (0.52,0.99),
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
    
    ## Total incompressible kinetic energy graph
    sca(e_ax)
    ylabel(r"$E_k/N_v$")
    e_ax.errorbar(t_vals, avg_en, yerr=avg_u_en, c="k", fmt="o", ls = '-', ms=5,mew = 0, capsize = 1, capthick =1)
    xlim(0,5.5)
    
    
    ## Dipole moment graph
    sca(dm_ax)
    ylabel('$d/R_\perp$')
    errorbar(t_vals, avg_d, yerr = avg_u_d, c="k", fmt="o", ls = '-', ms=5,mew = 0, capsize = 1, capthick =1)
    xlim(0,5.5)
    
    ## Correlation function graph
    sca(corr_ax)
    corr_ax.axhline(y=0, c= 'k', ls=':',lw=0.5)

    ylabel(r"$C_1$")
    errorbar(t_vals, avg_c1, yerr=avg_u_c1, c="k", fmt="o", label = "$C_1$", ls = '-', ms=5, mew = 0, capsize = 1, capthick =1)
    xlim(0,5.5)
    

    #### Done plotting the data, now tidy up the figure formatting ####
    
    ## Turn off tick labels on the x-axis of the stacked graphs
    class_ax.xaxis.set_ticklabels([])
    temp_ax.xaxis.set_ticklabels([])
    corr_ax.xaxis.set_ticklabels([])
    dm_ax.xaxis.set_ticklabels([])

    class_ax.tick_params(axis='x',which='minor',bottom='off', top="off")
    temp_ax.tick_params(axis='x',which='minor',bottom='off', top="off")
    corr_ax.tick_params(axis='x',which='minor',bottom='off', top="off")
    dm_ax.tick_params(axis='x',which='minor',bottom='off', top="off")
    e_ax.tick_params(axis='x',which='minor',bottom='off', top="off")

    ## label positioning common values
    figlabelx = -0.19
    figlabely = 0.68
    labelx = -0.15
    labely = -0.38

    ## Fix classification graph formatting
    sca(class_ax)
    ## Add the label
    class_ax.yaxis.set_label_coords(labelx, 0.5)
    class_ax.xaxis.set_label_coords(0.5, labely)
    class_ax.set_title('A', fontdict = blackfont,x=figlabelx,y=figlabely)
    ## Set the y-axis limits and tick marks depending on the grid (values chosen to best display data and provide sensible tick marks)
    ## 6.0
    if specify_sequence == '20171004T121555':
        yticks([0.2,0.4,0.6])
        ylim(0,0.8)
    ## 11.5
    elif specify_sequence == '20171006T101525':
        yticks([0.2,0.4,0.6])
        ylim(0,0.8)
    ## 4.2
    elif specify_sequence == '20171004T093812':
        yticks([0.2,0.4,0.6])
        ylim(0,0.8)
    ## 9.7
    elif specify_sequence == '20171004T144649':
        yticks([0.2,0.4,0.6])
        ylim(0,0.8)
    ## 7.9
    else:
        yticks([0.2,0.4,0.6])
        ylim(0,0.8)
    
    
    ## Fix correlation function graph formatting
    sca(corr_ax)
    ## Add the label
    corr_ax.yaxis.set_label_coords(labelx, 0.5)
    corr_ax.xaxis.set_label_coords(0.5, labely)
    corr_ax.set_title('B', fontdict = blackfont,x=figlabelx,y=figlabely)
    # since this graph has another on top of it, don't plot the top border, it will rely on the bottom border of the next graph up
    corr_ax.spines["top"].set_visible(False)
    
    ## Set the y-axis limits and tick marks depending on the grid (values chosen to best display data and provide sensible tick marks)
    ## 9.7
    if specify_sequence == '20171004T144649':
        yticks([0,0.2,0.4])
    ## 7.9    
    elif specify_sequence == '20171005T090729':
        yticks([-0.2,0,0.2,0.4])
    ## 11.5
    elif specify_sequence =='20171006T101525':
        yticks([0.1,0.3, 0.5])
    ## 4.2
    elif specify_sequence == '20171004T093812':
        yticks([-0.6,-0.4,-0.2,0,0.2])
    ## 6.0
    else:
        ylim(-0.3,0.35)
        yticks([-0.2,0,0.2])
    
    ## Fix dipole moment graph formatting
    sca(dm_ax)
    ## Add the label
    dm_ax.yaxis.set_label_coords(labelx, 0.5)
    dm_ax.xaxis.set_label_coords(0.5, labely)
    dm_ax.set_title('C', fontdict = blackfont,x=figlabelx,y=figlabely)
    # since this graph has another on top of it, don't plot the top border, it will rely on the bottom border of the next graph up
    dm_ax.spines["top"].set_visible(False)
    
    ## Set the y-axis limits and tick marks depending on the grid (values chosen to best display data and provide sensible tick marks)
    # 6.0
    if specify_sequence == '20171004T121555':
        yticks([0.1,0.2,0.3])
        ylim(0.05,0.35)
    ## 11.5
    elif specify_sequence == '20171006T101525':
        ylim(0.18,0.38)
        yticks([0.2,0.25,0.3,0.35])
    # 4.2
    elif specify_sequence == '20171004T093812':
        yticks([0.1,0.2,0.3])
    ## 9.7
    elif specify_sequence == '20171004T144649':
        yticks([0.22,0.26,0.3])
    # 7.9
    else:
        yticks([0.15,0.2,0.25,0.3])
    
    
    ## Fix temperature graph formatting
    sca(temp_ax)
    ## Add the label
    temp_ax.yaxis.set_label_coords(labelx, 0.5)
    temp_ax.xaxis.set_label_coords(0.5, labely)
    temp_ax.set_title('D', fontdict = blackfont,x=figlabelx,y=figlabely)
    # since this graph has another on top of it, don't plot the top border, it will rely on the bottom border of the next graph up
    temp_ax.spines["top"].set_visible(False)

    ## Set the y-axis limits and tick marks depending on the grid (values chosen to best display data and provide sensible tick marks)
    # 7.9 micron
    if specify_sequence == '20171005T090729':
        yticks([0,-0.2,-0.4,-0.6])
    # 11.5 micron
    elif specify_sequence =='20171006T101525':
        ylim(-0.18,-0.78)
        yticks([-0.3,-0.5,-0.7])
    # 4.2 micron
    elif specify_sequence == '20171004T093812':
        yticks([0.2,0,-0.2,-0.4,-0.6])
        ylim(0.4,-0.7)
    # 9.7 micron
    elif specify_sequence =='20171004T144649':
        ylim(-0.1,-0.75)
        yticks([-0.2,-0.4,-0.6])
    # 6.0 micron
    else:
        ylim(0.2,-0.65)
        yticks([0.1,-0.1,-0.3,-0.5])
        
        

    ## Fix energy graph formatting
    sca(e_ax)
    ## Add the label
    e_ax.yaxis.set_label_coords(labelx, 0.5)
    e_ax.xaxis.set_label_coords(0.5, labely)
    e_ax.set_title('E', fontdict = blackfont,x=figlabelx,y=figlabely)
    # since this graph has another on top of it, don't plot the top border, it will rely on the bottom border of the next graph up
    e_ax.spines["top"].set_visible(False)
    # also set the x-label for this graph, since it's on the bottom of the stack
    xlabel("Hold time (s)")
    ## Set the y-axis limits and tick marks depending on the grid (values chosen to best display data and provide sensible tick marks)
    ## 6.0
    if specify_sequence == '20171004T121555':
        ylim(2.25,3.35)
        yticks([2.4,2.6,2.8,3,3.2])
    ## 7.9
    elif specify_sequence == '20171005T090729':
        yticks([2.8,3.0,3.2])
    ## 4.2
    elif specify_sequence == '20171004T093812':
        yticks([2.4,2.8,3.2])
    ## 11.5
    elif specify_sequence =='20171006T101525':
        yticks([3.0,3.2,3.4,3.6])
    ## 9.7
    else:
        yticks([3,3.2,3.4])

    
    #### Now fix the spines on the stacked left-hand graphs
    spine_left_value = -0.8
    class_ax.spines['top'].set_bounds(spine_left_value, class_ax.viewLim.intervalx[1])
    class_ax.spines['bottom'].set_bounds(spine_left_value, class_ax.viewLim.intervalx[1])
    corr_ax.spines['bottom'].set_bounds(spine_left_value, corr_ax.viewLim.intervalx[1])
    dm_ax.spines['bottom'].set_bounds(spine_left_value, dm_ax.viewLim.intervalx[1])
    temp_ax.spines['bottom'].set_bounds(spine_left_value, temp_ax.viewLim.intervalx[1])
    e_ax.spines['bottom'].set_bounds(spine_left_value, e_ax.viewLim.intervalx[1])

    
    ## Fix energy spectrum graph formatting
    sca(spec_ax)
    ylabel(r'$E(k)/N_v$')
    xlabel(r'$k\xi$')

    ## Add vertical lines:
    ## Note, k_vec is in units of xi
    # k*xi = 1
    plt.axvline(x = 1, c = '0.5', ls = ':',zorder = -3)
    
    # system size
    k_sys = xi*2*pi/(circ_radius*pixel_size*2)
    plt.axvline(x = k_sys, c = '0.5',ls = '-',zorder = -3)
    
    # grid size
    k_grid_radius = xi*pi/grid_radius
    plt.axvline(x = k_grid_radius, c = '0.5', ls='--',zorder = -3)

    ## Generate guides to the eye at power-law scalings
    minus1 = 2e-7 * k_vec**(-1)
    minus3 =2e-7 *k_vec**(-3)
    minus53= 1.55e-7 *k_vec**(-5./3)

    plt.loglog(k_vec,minus1, c='k', ls = '-', lw=0.5, zorder = -3)
    plt.loglog(k_vec,minus3, c='k', ls = ':', lw=0.5, zorder = -3)
    plt.loglog(k_vec,minus53, c='k', ls = '--', lw=0.5, zorder = -3) 

    ## Generate a temporary figure environment for using to generate a colorbar for the energy spectra.
    fakefig = figure('Colorbar')
    cbar_show = imshow(array([[0,5]]),cmap="inferno_r")

    ## Go back to the spectrum axis (creating the figure above will have changed the current axis) and create the colorbar, based on the above imshow
    sca(spec_ax)
    cb = colorbar(cbar_show, cax=col_ax,ticks=[0,5],orientation="horizontal")
    cb.outline.set_linewidth(0.5)
    cb.ax.set_xticklabels(cb.ax.get_xticklabels(), fontsize=5, y=1)
    
    plt.text(colorbar_location[0]+0.5*colorbar_location[2],colorbar_location[1]+2*colorbar_location[3],"Hold time (s)", fontsize = 5,ha ='center',transform=spec_ax.transAxes)
    
    ## Now format the figure itself
    xlim(espec_x_min,espec_x_max)
    ylim(1e-9,4e-6)
    yticks([1e-8,1e-7,1e-6])

    spec_ax.yaxis.set_label_coords(-0.1, 0.5)
    spec_ax.xaxis.set_label_coords(0.5, -0.1)
    spec_ax.set_title('F', fontdict = blackfont,x=-0.14,y=0.88)

    spec_ax.tick_params(axis = 'y', direction='in', pad=1)
    spec_ax.tick_params(axis = 'x', direction='in', pad=1.5)

    ## Fix vortex number graph formatting
    sca(logN_ax)
    logN_ax.yaxis.set_label_coords(-0.22, 0.5)
    logN_ax.xaxis.set_label_coords(0.5, -0.24)
    logN_ax.set_title('G', fontdict = blackfont,x=-0.32,y=0.8)

    logN_ax.xaxis.set_major_formatter(ScalarFormatter())
    logN_ax.yaxis.set_major_formatter(ScalarFormatter())
    yticks([4,6,10,20])
    xticks([0.5,1,2,5])
    
    ## Fix vortex spacing graph formatting
    sca(logL_ax)
    logL_ax.yaxis.set_label_coords(-0.22, 0.5)
    logL_ax.xaxis.set_label_coords(0.5, -0.24)
    logL_ax.set_title('H', fontdict = blackfont,x=-0.32,y=0.8)

    logL_ax.xaxis.set_major_formatter(ScalarFormatter())
    logL_ax.yaxis.set_major_formatter(ScalarFormatter())
    yticks([0.2,0.3,0.4,0.5])
    xticks([0.5,1,2,5])
    
    ## Fix the inset in the energy spectrum graph's formatting
    sca(delta_e_ax)
    title(r'$\Delta E(k)/N_v \times 10^7$', fontdict = tinyfont, x=0.2,y=0.9)
    
    xlim(espec_inset_x_min,espec_inset_x_max)
    delta_e_ax.set_xscale("log", nonposx='clip')
    locator_params(axis='y',nbins=4)
    
    delta_e_ax.tick_params(axis = 'both', length = 2)
    delta_e_ax.tick_params(axis = 'x', which = 'minor', length = 1)
    
    ## Make everything a bit smaller for this inset!
    for label in (delta_e_ax.get_xticklabels() + delta_e_ax.get_yticklabels()):
        label.set_fontsize(tinyfont['size'])
    
    [i.set_linewidth(0.5) for i in delta_e_ax.spines.itervalues()]
    
    delta_e_ax.xaxis.set_major_formatter(NullFormatter())
    delta_e_ax.patch.set_alpha(0)
    delta_e_ax.yaxis.offsetText.set_fontsize(tinyfont['size'])
    delta_e_ax.tick_params(axis = 'y', direction='in', pad=1)
    
    ## Finally, adjust the spacing and padding on the GridSpecs to make the figure look right
    gs.update(left=0.085, right=0.85, top = 0.995, bottom = 0.1, wspace=0, hspace = 0.02)
    gs_log.update(left=0.12, right=0.995, top = 0.995, bottom = 0.1, wspace=0.35, hspace = 1.5)
    ## Make sure we've got the correct figure, before saving the pdf
    figure('dynamics %s'%specify_sequence)
    
    savefig(os.path.join(results_path,'individual_run_%s.pdf'%int(grid_radius*1e6)),transparent = True)


    
for specify_sequence in sequence_list:
    plot_time_series(specify_sequence)
# Close the spreadsheet
workbook.close()

## If the script is being run from outside lyse, show the graphs:
if not lyse.spinning_top: show()