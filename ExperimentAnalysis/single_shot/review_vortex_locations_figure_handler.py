###===============review_vortex_locations_figure_handler.py=================###
###                                                                         ###
### This script is not designed to be run itself, but have the class        ###
### VortexPoints imported by review_vortex_locations.py to provide tools    ###
### for adding and subtracting vortices interactively via a figure in lyse. ###
###                                                                         ###
### Copyright (C) 2018  Shaun Johnstone and Philip Starkey                  ###
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

from lyse import *
from qtutils.qt import QtWidgets, QtGui


class VortexPoints(Plot):
    def __init__(self, *args, **kwargs):
        Plot.__init__(self, *args, **kwargs)
        self.points_of_interest_action = self.navigation_toolbar.addAction(
            QtGui.QIcon(':qtutils/fugue/target--plus'),
           'Add points of interest',self.on_add_points_of_interest_triggered)
        self.points_of_interest_action.setToolTip('Add points of interest')
        self.points_of_interest_action.setCheckable(True)
        
        self.actions_enabled = True
        self.canvas.mpl_connect('button_press_event',self.canvasClicked)

        self.remove_points_of_interest_action = self.navigation_toolbar.addAction(
            QtGui.QIcon(':qtutils/fugue/target--minus'),
           'Remove points of interest',self.on_remove_points_of_interest_triggered)
        self.remove_points_of_interest_action.setToolTip('Remove points of interest')
        self.remove_points_of_interest_action.setCheckable(True)
        
        self.canvas.mpl_connect('pick_event', self.onPick)
        
        
        self.clear_points_of_interest_action = self.navigation_toolbar.addAction(
            QtGui.QIcon(':qtutils/fugue/target--exclamation'),
           'Clear all points of interest',self.on_clear_points_of_interest_triggered)
        self.clear_points_of_interest_action.setToolTip('Cear all points of interest')
        
        
        self.done_action = self.navigation_toolbar.addAction(QtWidgets.QIcon(':qtutils/fugue/disk--arrow'), 'Done', self.save_points_of_interest)
        self.done_action.setToolTip('Save vortex locations')
        
        self.reject_action = self.navigation_toolbar.addAction(QtWidgets.QIcon(':qtutils/fugue/disk--minus'), 'Reject', self.on_reject)
        self.reject_action.setToolTip('Reject shot')
        
        self.cancel_action = self.navigation_toolbar.addAction(QtWidgets.QIcon(':qtutils/fugue/disk--exclamation'), 'Cancel', self.on_reject)
        self.cancel_action.setToolTip('Cancel (do not save any results for this shot)')

    def on_add_points_of_interest_triggered(self):
        if self.points_of_interest_action.isChecked():            
            if self.navigation_toolbar._active == "PAN":
                self.navigation_toolbar.pan()
            elif self.navigation_toolbar._active == "ZOOM":
                self.navigation_toolbar.zoom()
            if self.remove_points_of_interest_action.isChecked():
                self.remove_points_of_interest_action.setChecked(False)
            self.navigation_toolbar.set_message('Add points of interest')
    
    def on_remove_points_of_interest_triggered(self):
        if self.remove_points_of_interest_action.isChecked():            
            if self.navigation_toolbar._active == "PAN":
                self.navigation_toolbar.pan()
            elif self.navigation_toolbar._active == "ZOOM":
                self.navigation_toolbar.zoom()
            if self.points_of_interest_action.isChecked():
                self.points_of_interest_action.setChecked(False)
            self.navigation_toolbar.set_message('Remove points of interest')
            
    def canvasClicked( self, event ):
        if self.points_of_interest_action.isChecked() and self.actions_enabled:
            button=event.button
            x=event.xdata
            y=event.ydata
            ax = event.inaxes
            ax.plot(x,y,'ro',label = "POI", picker = True)
            self.draw()
            
    def onPick(self,event):
        if self.remove_points_of_interest_action.isChecked() and self.actions_enabled:
            if event.artist.get_label() == "POI":
                event.artist.remove()
                self.draw()
    def on_clear_points_of_interest_triggered(self):
        for i, ax in enumerate(self.figure.axes):
            artists = ax.get_children()
            for child in artists:
                if child.get_label() == "POI":
                    try:
                        child.remove()
                    except:
                        pass
        self.draw()
    
    def disable_actions(self):
        self.actions_enabled = False
        self.points_of_interest_action.setEnabled(False)
        self.remove_points_of_interest_action.setEnabled(False)
        self.clear_points_of_interest_action.setEnabled(False)
        self.done_action.setEnabled(False)
        self.reject_action.setEnabled(False)
        self.cancel_action.setEnabled(False)
        
    def enable_actions(self):
        self.actions_enabled = True
        self.points_of_interest_action.setEnabled(True)
        self.remove_points_of_interest_action.setEnabled(True)
        self.clear_points_of_interest_action.setEnabled(True)
        self.done_action.setEnabled(True)
        self.reject_action.setEnabled(True)
        self.cancel_action.setEnabled(True)
    
    def save_points_of_interest(self):
        self.disable_actions()
        points_of_interest = {}
        for i, ax in enumerate(self.figure.axes):
            points = []
            artists = ax.get_children()
            for child in artists:
                if child.get_label() == "POI":
                    try:
                        points.append(child.get_xydata())
                    except:
                        pass       
            points_of_interest[i] = points
        points_of_interest = points_of_interest[0]
        routine_storage.run.save_result_array('vortex_points',points_of_interest)
        routine_storage.run.save_result('N_v',len(points_of_interest))
        routine_storage.run.save_result('Good',True)
        delay_event.set()
   
        
        
    def on_reject(self):
        self.disable_actions()
        routine_storage.run.save_result_array('vortex_points',[])
        routine_storage.run.save_result('N_v',False)
        routine_storage.run.save_result('Good',False)
        delay_event.set()
        
    def on_cancel(self):
        self.disable_actions()
        delay_event.set()
    
    def analysis_complete(self, *args, **kwargs):
        self.enable_actions()