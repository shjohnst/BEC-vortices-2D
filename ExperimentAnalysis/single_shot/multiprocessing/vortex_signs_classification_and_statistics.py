###======multiprocessing/vortex_signs_classification_and_statistics.py======###
###                                                                         ###
### This script is a wrapper function for the single-shot analysis routine  ###
### of the same name, that splits the vortex configuration space across     ###
### multiple CPU cores to be run in parallel.                               ###
###                                                                         ###
### Copyright (C) 2018 Philip Starkey and Shaun Johnstone                   ###
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



### NOTE: You must have 2^n_vortices >= num_processes
### Adjust n_vortices_single_process in the base "vortex_signs_classification_and_statistics.py"
### where appropriate to find best balance between single- and multi-process analysis for your
### CPU speed and number of cores.

from lyse import *
import time
import itertools
import ast
import sys, os
from pylab import *
from psutil import Popen, cpu_count
from subprocess import PIPE, CREATE_NEW_CONSOLE
# from labscript_utils.labconfig import LabConfig
num_processes = cpu_count()
print "Using %s processes"%num_processes

import time

from analysislib.single_shot.vortex_signs_classification_and_statistics import run_code, n_vortices_single_process

multiprocess = True
with h5py.File(path,'r') as h5_file:
    if "results/vortex_signs_classification_and_statistics" in h5_file:
        N = h5_file["results"]["vortex_signs_classification_and_statistics"].attrs.get('N_v',None)
        if N is not None and N <= n_vortices_single_process:
            multiprocess = False
            print N
            print "Skipping, has been analysed with single process"
if multiprocess:
    start_time = time.time()
    p_list = []
    p_launch_times = []
    
    ## Get analysislib to find the path to set as the current working directory
    # exp_config = LabConfig()
    # analysislibpath = exp_config.get('paths', 'analysislib')
    # scriptdir = os.path.join(analysislibpath,'single_shot')
    scriptdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    for i, config in enumerate(map(list, itertools.product([-1, 1], repeat=int(log2(num_processes))))):
        p_launch_times.append(time.time())
        p = Popen(['python.exe', 'vortex_signs_classification_and_statistics.py', str(path), str(config)], cwd=scriptdir, stdout=PIPE, stderr=PIPE)#, creationflags=CREATE_NEW_CONSOLE)
        p_list.append(p)
        print "process %d PID=%d"%(i, p.pid)
        
    result = []
    for i, p in enumerate(p_list):
        stdout, stderr = p.communicate()
        print 'subprocess %d output (note all stdout is printed first, followed by all stderr)'%i 
        print 'launched at: ', p_launch_times[i]
        print stdout.replace('\r', '')
        sys.stderr.write( stderr.replace('\r', ''))
        sys.stderr.flush()
        print '----'

        
        try:
            result.append(ast.literal_eval(stdout.split('\r\n')[-2]))
        except:
            sys.stderr.write('Failed to process result from subprocess %d'%i)
            sys.stderr.flush()
        
    end_time = time.time()
        
    print 'Best configurations that were returned:'
    for r in result:
        print r
    # print result
    print "subprocess time: ",end_time-start_time
    print "Now determining the best configuration from the candidate configurations returned by each process..."
    run = Run(path)
    run_code(result, run)