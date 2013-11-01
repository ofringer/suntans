#!/usr/bin/python
"""
SUNTANS NetCDF plotting GUI
"""
import os
import wx

# The recommended way to use wx with mpl is with the WXAgg
# backend. 
#
import matplotlib
matplotlib.use('WXAgg')
from matplotlib.figure import Figure
from matplotlib.collections import PolyCollection
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar

from sunpy import Spatial
from datetime import datetime
import numpy as np

import pdb

class SunPlotPy(wx.Frame, Spatial):
    """ 
    The main frame of the application
    """
    title = 'sunplot(py)'

    # Plotting options
    autoclim=True
    showedges=False
    bgcolor='k'
    textcolor='w'
    
    def __init__(self):
        wx.Frame.__init__(self, None, -1, self.title)
        
        self.create_menu()
        self.create_status_bar()
        self.create_main_panel()
        
        #self.draw_figure()

    def create_menu(self):
        self.menubar = wx.MenuBar()
        
        menu_file = wx.Menu()
        m_expt = menu_file.Append(-1, "&Open file\tCtrl-O", "Open netcdf file")
        self.Bind(wx.EVT_MENU, self.on_open_file, m_expt)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)
        
        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&About\tF1", "About the demo")
        self.Bind(wx.EVT_MENU, self.on_about, m_about)
        
        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_help, "&Help")
        self.SetMenuBar(self.menubar)

    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas 
             * mpl navigation toolbar
             * Control panel for interaction
        """
        self.panel = wx.Panel(self)
        
        # Create the mpl Figure and FigCanvas objects. 
        # 5x4 inches, 100 dots-per-inch
        #
        self.dpi = 100
        self.fig = Figure((7.0, 6.0), dpi=self.dpi,facecolor=self.bgcolor)
        self.canvas = FigCanvas(self.panel, -1, self.fig)
        
        # Since we have only one plot, we can use add_axes 
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.axes = self.fig.add_subplot(111)
        SetAxColor(self.axes,self.textcolor,self.bgcolor)
        
        # Bind the 'pick' event for clicking on one of the bars
        #
        #self.canvas.mpl_connect('pick_event', self.on_pick)
        
        ########
        # Create widgets
        ########
        self.variable_list = wx.ComboBox(
            self.panel, 
            size=(200,-1),
            choices=['Select a variable...'],
            style=wx.CB_READONLY)
        self.variable_list.Bind(wx.EVT_COMBOBOX, self.on_select_variable)
        
        self.time_list = wx.ComboBox(
            self.panel, 
            size=(200,-1),
            choices=['Select a time step...'],
            style=wx.CB_READONLY)
        self.time_list.Bind(wx.EVT_COMBOBOX, self.on_select_time)

        self.depthlayer_list = wx.ComboBox(
            self.panel, 
            size=(200,-1),
            choices=['Select a vertical layer...'],
            style=wx.CB_READONLY)
        self.depthlayer_list.Bind(wx.EVT_COMBOBOX, self.on_select_depth)

        self.show_edge_check = wx.CheckBox(self.panel, -1, 
            "Show Edges",
            style=wx.ALIGN_RIGHT)
        self.show_edge_check.Bind(wx.EVT_CHECKBOX, self.on_show_edges)

        cmaps = matplotlib.cm.datad.keys()
        cmaps.sort()
        self.colormap_list = wx.ComboBox(
            self.panel, 
            size=(100,-1),
            choices=cmaps,
            style=wx.CB_READONLY)
        self.colormap_list.Bind(wx.EVT_COMBOBOX, self.on_select_cmap)
        self.colormap_label = wx.StaticText(self.panel, -1,"Colormap:")

        self.clim_check = wx.CheckBox(self.panel, -1, 
            "Manual color limits ",
            style=wx.ALIGN_RIGHT)
        self.clim_check.Bind(wx.EVT_CHECKBOX, self.on_clim_check)

        self.climlow = wx.TextCtrl(
            self.panel, 
            size=(100,-1),
            style=wx.TE_PROCESS_ENTER)
        self.climlow.Bind(wx.EVT_TEXT_ENTER, self.on_climlow)
        
        self.climhigh = wx.TextCtrl(
            self.panel, 
            size=(100,-1),
            style=wx.TE_PROCESS_ENTER)
        self.climhigh.Bind(wx.EVT_TEXT_ENTER, self.on_climhigh)
 


        # Labels
        self.variable_label = wx.StaticText(self.panel, -1,"Variable:",size=(200,-1))
        self.time_label = wx.StaticText(self.panel, -1,"Time step:",size=(200,-1))
        self.depth_label = wx.StaticText(self.panel, -1,"Vertical level:",size=(200,-1))


        # Create the navigation toolbar, tied to the canvas
        #
        self.toolbar = NavigationToolbar(self.canvas)
        
        #########
        # Layout with box sizers
        #########
        
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)

        self.vbox.AddSpacer(10)
        #self.vbox.Add((-1,25))

        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL

        self.hbox0 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox0.Add(self.show_edge_check, 0, border=10, flag=flags)
        self.hbox0.Add(self.colormap_label, 0, border=10, flag=flags)
        self.hbox0.Add(self.colormap_list, 0, border=10, flag=flags)
        self.hbox0.Add(self.clim_check, 0, border=10, flag=flags)
        self.hbox0.Add(self.climlow, 0, border=10, flag=flags)
        self.hbox0.Add(self.climhigh, 0, border=10, flag=flags)

        self.vbox.AddSpacer(5)
        self.hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox1.Add(self.variable_label, 0, border=10, flag=flags)
        self.hbox1.Add(self.time_label, 0, border=10, flag=flags)
        self.hbox1.Add(self.depth_label, 0, border=10, flag=flags)

        self.vbox.AddSpacer(5)
        self.hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox2.Add(self.variable_list, 0, border=10, flag=flags)
        self.hbox2.Add(self.time_list, 0, border=10, flag=flags)
        self.hbox2.Add(self.depthlayer_list, 0, border=10, flag=flags)
       
        self.vbox.Add(self.hbox1, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        self.vbox.Add(self.hbox2, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        self.vbox.Add(self.hbox0, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        
        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)
    
    ##########
    # Event functions
    ##########

    def create_figure(self):
        """ 
        Creates the figure
        """
        # Find the colorbar limits if unspecified
        if self.autoclim:
            self.clim = [self.data.min(),self.data.max()]
            self.climlow.SetValue('%3.1f'%self.clim[0])
            self.climhigh.SetValue('%3.1f'%self.clim[1])
         
        if self.__dict__.has_key('collection'):
            self.collection.remove()

        self.collection = PolyCollection(self.xy,)
        self.collection.set_array(np.array(self.data[:]))
        self.collection.set_clim(vmin=self.clim[0],vmax=self.clim[1])
        if not self.showedges:
            self.collection.set_edgecolors(self.collection.to_rgba(np.array((self.data[:])))) 
        self.axes.add_collection(self.collection)    
        self.axes.set_aspect('equal')
        self.axes.set_xlim(self.xlims)
        self.axes.set_ylim(self.ylims)
        title=self.axes.set_title(self.genTitle(),color=self.textcolor)
        self.axes.set_xlabel('Easting [m]')
        self.axes.set_ylabel('Northing [m]')

        # create a colorbar

        if not self.__dict__.has_key('cbar'):
            self.cbar = self.fig.colorbar(self.collection)
            SetAxColor(self.cbar.ax.axes,self.textcolor,self.bgcolor)
        else:
            pass
            #self.cbar.check_update(self.collection)
            #self.cbar.on_mappable_changed(self.collection)

        self.canvas.draw()
   
    def update_figure(self):
        if self.autoclim:
            self.clim = [self.data.min(),self.data.max()]
            self.climlow.SetValue('%3.1f'%self.clim[0])
            self.climhigh.SetValue('%3.1f'%self.clim[1])
 
        self.collection.set_array(np.array(self.data[:]))
        self.collection.set_clim(vmin=self.clim[0],vmax=self.clim[1])
        if not self.showedges:
            self.collection.set_edgecolors(self.collection.to_rgba(np.array((self.data[:])))) 
        else:
            self.collection.set_edgecolors('k')

        title=self.axes.set_title(self.genTitle(),color=self.textcolor)
        self.canvas.draw()
    
    def on_pick(self, event):
        # The event received here is of the type
        # matplotlib.backend_bases.PickEvent
        #
        # It carries lots of information, of which we're using
        # only a small amount here.
        # 
        box_points = event.artist.get_bbox().get_points()
        msg = "You've clicked on a bar with coords:\n %s" % box_points
        
        dlg = wx.MessageDialog(
            self, 
            msg, 
            "Click!",
            wx.OK | wx.ICON_INFORMATION)

        dlg.ShowModal() 
        dlg.Destroy()        
    
    def on_select_variable(self, event):
        vname = event.GetString()
        self.flash_status_message("Selecting variable: %s"%vname)
        # update the spatial object and load the data
        self.variable = vname
        self.loadData()

        # Check if the variable has a depth coordinate
        depthstr = ['']
        # If so populate the vertical layer box
        if self.hasDim('Nk'):
            depthstr = ['%3.1f'%self.z_r[k] for k in range(self.Nkmax)]
        elif self.hasDim('Nkw'):
            depthstr = ['%3.1f'%self.z_w[k] for k in range(self.Nkmax+1)]

        self.depthlayer_list.SetItems(depthstr)

        # Update the plot
        self.update_figure()


    def on_select_time(self, event):
        tindex = event.GetSelection()
        # Update the object time index and reload the data
        if not self.tstep==tindex:
            self.tstep=tindex
            self.loadData()
            self.flash_status_message("Selecting variable: %s..."%event.GetString())

            # Update the plot
            self.update_figure()


    def on_select_depth(self, event):
        kindex = event.GetSelection()
        if not self.klayer[0]==kindex:
            self.klayer = [kindex]
            self.loadData()       
            self.flash_status_message("Selecting depth: %s..."%event.GetString())

            # Update the plot
            self.update_figure()

    def on_open_file(self, event):
        file_choices = "NetCDF (*.nc)|*.nc|All Files (*.*)|*.*"
        
        dlg = wx.FileDialog(
            self, 
            message="Save plot as...",
            defaultDir=os.getcwd(),
            defaultFile="",
            wildcard=file_choices,
            style= wx.FD_MULTIPLE)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPaths()
            self.flash_status_message("Opening file: %s" % path)
            
            # Initialise the class
            Spatial.__init__(self,path)
            
            # Populate the drop down menus
            vnames = self.listCoordVars()
            self.variable_list.SetItems(vnames)
            
            # Update the time drop down list
            self.timestr = [datetime.strftime(tt,'%d-%b-%Y %H:%M:%S') for tt in self.time]
            self.time_list.SetItems(self.timestr)

            # Draw the depth
            if 'dv' in vnames:
                self.variable='dv'
                self.loadData()
                self.create_figure()
        
    def on_show_edges(self,event):
        sender=event.GetEventObject()
        self.showedges = sender.GetValue()

        # Update the figure
        self.update_figure()

    def on_clim_check(self,event):
        sender=event.GetEventObject()
        if sender.GetValue() == True:
            self.autoclim=False
            self.update_figure()
        else:
            self.autoclim=True
       

    def on_climlow(self,event):
        self.clim[0] = event.GetString()
        #self.update_figure()

    def on_climhigh(self,event):
        self.clim[1] = event.GetString()
        #self.update_figure()

    def on_select_cmap(self,event):
        cmap=event.GetString()
        self.collection.set_cmap(cmap)

        # Update the figure
        self.update_figure()



    def on_exit(self, event):
        self.Destroy()
        
    def on_about(self, event):
        msg = """ SUNTANS NetCDF visualization tool
        
            *Author: Matt Rayson
            *Institution: Stanford University
            *Created: October 2013
        """
        dlg = wx.MessageDialog(self, msg, "About", wx.OK)
        dlg.ShowModal()
        dlg.Destroy()
   
    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()

    def flash_status_message(self, msg, flash_len_ms=1500):
        self.statusbar.SetStatusText(msg)
        self.timeroff = wx.Timer(self)
        self.Bind(
            wx.EVT_TIMER, 
            self.on_flash_status_off, 
            self.timeroff)
        self.timeroff.Start(flash_len_ms, oneShot=True)
    
    def on_flash_status_off(self, event):
        self.statusbar.SetStatusText('')


def SetAxColor(ax,color,bgcolor):
    ax.set_axis_bgcolor(bgcolor)
    
    ax.yaxis.set_tick_params(color=color,labelcolor=color)
    ax.xaxis.set_tick_params(color=color,labelcolor=color)
    ax.yaxis.label.set_color(color)
    ax.xaxis.label.set_color(color)
    ax.spines['top'].set_color(color)
    ax.spines['bottom'].set_color(color)
    ax.spines['left'].set_color(color)
    ax.spines['right'].set_color(color)
    return ax
 

if __name__ == '__main__':
    app = wx.PySimpleApp()
    app.frame = SunPlotPy()
    app.frame.Show()
    app.MainLoop()

