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
# Set some default parameters
matplotlib.rcParams['text.color']='white'
matplotlib.rcParams['savefig.facecolor']='black'
matplotlib.rcParams['savefig.edgecolor']='black'
matplotlib.rcParams['figure.facecolor']='black'
matplotlib.rcParams['figure.edgecolor']='black'
matplotlib.rcParams['axes.facecolor']='black'
matplotlib.rcParams['axes.edgecolor']='white'
matplotlib.rcParams['axes.labelcolor']='white'
matplotlib.rcParams['xtick.color']='white'
matplotlib.rcParams['ytick.color']='white'
#matplotlib.rcParams['font.family']='serif'
#matplotlib.rcParams['font.sans-serif']=['Verdana']

from matplotlib.figure import Figure
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar
import matplotlib.animation as animation

from sunpy import Spatial, Grid
from untrim_tools import untrim_gridvars, untrim_griddims, UNTRIMSpatial
from ptm_tools import PtmBin
from datetime import datetime
import numpy as np

import pdb

class SunPlotPy(wx.Frame, Spatial, Grid ):
    """ 
    The main frame of the application
    """
    title = 'sunplot(py)'

    # Plotting options
    autoclim=True
    showedges=False
    bgcolor='k'
    textcolor='w'
    cmap='jet'

    # other flags
    collectiontype='cells'
    oldcollectiontype='cells'

    tindex=0 
    
    def __init__(self):
        wx.Frame.__init__(self, None, -1, self.title)
        
        self.create_menu()
        self.create_status_bar()
        self.create_main_panel()
        
        #self.draw_figure()

    def create_menu(self):
        self.menubar = wx.MenuBar()
        
        ###
        # File Menu
        ###
        menu_file = wx.Menu()
        # Load a hydro output file
        m_expt = menu_file.Append(-1, "&Open file\tCtrl-O", "Open netcdf file")
        self.Bind(wx.EVT_MENU, self.on_open_file, m_expt)

        # Load a grid file
        m_grid = menu_file.Append(-1, "&Load grid\tCtrl-G", "Load SUNTANS grid from folder")
        self.Bind(wx.EVT_MENU, self.on_load_grid, m_grid)

        # Load a particle file
        m_part = menu_file.Append(-1, "&Load PTM file\tCtrl-Shift-P", "Load a PTM file")
        self.Bind(wx.EVT_MENU, self.on_load_ptm, m_part)

        # Save current scene as an animation
        m_anim = menu_file.Append(-1,"&Save animation of current scene\tCtrl-S","Save animation")
        self.Bind(wx.EVT_MENU, self.on_save_anim, m_anim)

        # Save the current figure
        m_prin = menu_file.Append(-1,"&Print current scene\tCtrl-P","Save figure")
        self.Bind(wx.EVT_MENU, self.on_save_fig, m_prin)



        menu_file.AppendSeparator()
        # Exit
        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)

        ###
        # Tools menu
        ###
        menu_tools = wx.Menu()
        m_gridstat = menu_tools.Append(-1, "&Plot grid size statistics", "SUNTANS grid size")
        self.Bind(wx.EVT_MENU, self.on_plot_gridstat, m_gridstat)

        
        ###
        # Help Menu
        ###
        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&About\tF1", "About the demo")
        self.Bind(wx.EVT_MENU, self.on_about, m_about)
        
        
        # Add all of the menu bars
        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_tools, "&Tools")
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
        #self.fig = Figure((7.0, 6.0), dpi=self.dpi,facecolor=self.bgcolor)
        self.fig = Figure((7.0, 6.0), dpi=self.dpi)
        self.canvas = FigCanvas(self.panel, -1, self.fig)
        
        
        # Since we have only one plot, we can use add_axes 
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.axes = self.fig.add_subplot(111)
        #SetAxColor(self.axes,self.textcolor,self.bgcolor)
        
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
        #self.toolbar.toolitems[8][3]='my_save_fig'

        #def my_save_fig(self,*args):
        #    print 'saving figure'
        #    return "break"

        
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
            #self.collection.remove()
            self.axes.collections.remove(self.collection)
        else:
            # First call - set the axes limits
            self.axes.set_aspect('equal')
            self.axes.set_xlim(self.xlims)
            self.axes.set_ylim(self.ylims)
 

        if self.collectiontype=='cells':
            self.collection = PolyCollection(self.xy,cmap=self.cmap)
            self.collection.set_array(np.array(self.data[:]))
            if not self.showedges:
                self.collection.set_edgecolors(self.collection.to_rgba(np.array((self.data[:])))) 
        elif self.collectiontype=='edges':
            xylines = [self.xp[self.edges],self.yp[self.edges]]
            linesc = [zip(xylines[0][ii,:],xylines[1][ii,:]) for ii in range(self.Ne)]
            self.collection = LineCollection(linesc,array=np.array(self.data[:]),cmap=self.cmap)

        self.collection.set_clim(vmin=self.clim[0],vmax=self.clim[1])

        self.axes.add_collection(self.collection)    
        self.title=self.axes.set_title(self.genTitle(),color=self.textcolor)
        self.axes.set_xlabel('Easting [m]')
        self.axes.set_ylabel('Northing [m]')

        # create a colorbar

        if not self.__dict__.has_key('cbar'):
            self.cbar = self.fig.colorbar(self.collection)
            #SetAxColor(self.cbar.ax.axes,self.textcolor,self.bgcolor)
        else:
            #pass
            print 'Updating colorbar...'
            #self.cbar.check_update(self.collection)
            self.cbar.on_mappable_changed(self.collection)

        self.canvas.draw()
   
    def update_figure(self):
        if self.autoclim:
            self.clim = [self.data.min(),self.data.max()]
            self.climlow.SetValue('%3.1f'%self.clim[0])
            self.climhigh.SetValue('%3.1f'%self.clim[1])
        else:
            self.clim = [float(self.climlow.GetValue()),\
                float(self.climhigh.GetValue())]
 
        # check whether it is cell or edge type
        if self.hasDim(self.variable,self.griddims['Ne']):
            self.collectiontype='edges'
        elif self.hasDim(self.variable,self.griddims['Nc']):
            self.collectiontype='cells'

        # Create a new figure if the variable has gone from cell to edge of vice
        # versa
        if not self.collectiontype==self.oldcollectiontype:
            self.create_figure()
            self.oldcollectiontype=self.collectiontype

        self.collection.set_array(np.array(self.data[:]))
        self.collection.set_clim(vmin=self.clim[0],vmax=self.clim[1])

        # Cells only
        if self.collectiontype=='cells':
            if not self.showedges:
                self.collection.set_edgecolors(self.collection.to_rgba(np.array((self.data[:])))) 
            else:
                self.collection.set_edgecolors('k')
                self.collection.set_linewidths(0.2)

        # Update the title
        self.title=self.axes.set_title(self.genTitle(),color=self.textcolor)

        #Update the colorbar
        self.cbar.update_normal(self.collection)

        # redraw the figure
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
        self.loadData(variable=self.variable)

        # Check if the variable has a depth coordinate
        depthstr = ['']
        # If so populate the vertical layer box
        if self.hasDim(self.variable,self.griddims['Nk']):
            depthstr = ['%3.1f'%self.z_r[k] for k in range(self.Nkmax)]
            depthstr += ['surface','seabed']

        elif self.hasDim(self.variable,'Nkw'):
            depthstr = ['%3.1f'%self.z_w[k] for k in range(self.Nkmax+1)]

        self.depthlayer_list.SetItems(depthstr)

        # Update the plot
        self.update_figure()



    def on_select_time(self, event):
        self.tindex = event.GetSelection()
        # Update the object time index and reload the data
        if self.plot_type=='hydro':
            if not self.tstep==self.tindex:
                self.tstep=self.tindex
                self.loadData()
                self.flash_status_message("Selecting variable: %s..."%event.GetString())

                # Update the plot
                self.update_figure()
        elif self.plot_type=='particles':
            self.PTM.plot(self.tindex,ax=self.axes,\
                xlims=self.axes.get_xlim(),ylims=self.axes.get_ylim())
        
            self.canvas.draw()


    def on_select_depth(self, event):
        kindex = event.GetSelection()
        if not self.klayer[0]==kindex:
            # Check if its the seabed or surface value
            if kindex>=self.Nkmax:
                kindex=event.GetString()
            self.klayer = [kindex]
            self.loadData()       
            self.flash_status_message("Selecting depth: %s..."%event.GetString())

            # Update the plot
            self.update_figure()

    def on_open_file(self, event):
        file_choices = "SUNTANS NetCDF (*.nc)|*.nc*|UnTRIM NetCDF (*.nc)|*.nc*|All Files (*.*)|*.*"
        
        dlg = wx.FileDialog(
            self, 
            message="Open SUNTANS file...",
            defaultDir=os.getcwd(),
            defaultFile="",
            wildcard=file_choices,
            style= wx.FD_MULTIPLE)
        
        if dlg.ShowModal() == wx.ID_OK:
            self.plot_type='hydro'

            path = dlg.GetPaths()

            # Initialise the class
            if dlg.GetFilterIndex() == 0: #SUNTANS
                self.flash_status_message("Opening SUNTANS file: %s" % path)
                Spatial.__init__(self,path)
                startvar='dv'
            if dlg.GetFilterIndex()==1: #UnTRIM
                self.flash_status_message("Opening UnTRIMS file: %s" % path)
                #Spatial.__init__(self,path,gridvars=untrim_gridvars,griddims=untrim_griddims)
                UNTRIMSpatial.__init__(self,path)
                startvar='Mesh2_face_depth'
            
            # Populate the drop down menus
            vnames = self.listCoordVars()
            self.variable_list.SetItems(vnames)
            
            # Update the time drop down list
            self.timestr = [datetime.strftime(tt,'%d-%b-%Y %H:%M:%S') for tt in self.time]
            self.time_list.SetItems(self.timestr)

            # Draw the depth
            if startvar in vnames:
                self.variable=startvar
                self.loadData()
                self.create_figure()

    def on_load_grid(self, event):
        
        dlg = wx.DirDialog(
            self, 
            message="Open SUNTANS grid from folder...",
            defaultPath=os.getcwd(),
            style= wx.DD_DEFAULT_STYLE)
        
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()

            # Initialise the class
            self.flash_status_message("Opening SUNTANS grid from folder: %s" % path)
            Grid.__init__(self,path)

            # Plot the Grid
            self.axes,self.collection = self.plotmesh(ax=self.axes,edgecolors='y')

            # redraw the figure
            self.canvas.draw()

    def on_load_ptm(self, event):
        file_choices = "PTM Binary (*_bin.out)|*_bin.out|All Files (*.*)|*.*"
        
        dlg = wx.FileDialog(
            self, 
            message="Open PTM file...",
            defaultDir=os.getcwd(),
            defaultFile="",
            wildcard=file_choices,
            style= wx.FD_MULTIPLE)
        
        if dlg.ShowModal() == wx.ID_OK:
            self.plot_type = 'particles'
            path = dlg.GetPath()

            # Initialise the class
            self.flash_status_message("Opening PTM binary file: %s" % path)
            self.PTM = PtmBin(path)
            
            # Update the time drop down list
            self.timestr = [datetime.strftime(tt,'%d-%b-%Y %H:%M:%S') for tt in self.PTM.time]
            self.time_list.SetItems(self.timestr)

            # Plot the first time step
            if self.__dict__.has_key('xlims'):
                self.PTM.plot(self.PTM.nt-1,ax=self.axes,xlims=self.xlims,\
                ylims=self.ylims,fontcolor='w')
            else:
                self.PTM.plot(self.PTM.nt-1,ax=self.axes,fontcolor='w')
            # redraw the figure
            self.canvas.draw()

        
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
        self.cmap=event.GetString()
        self.collection.set_cmap(self.cmap)

        # Update the figure
        self.update_figure()

    def on_save_fig(self,event):
        """
        Save a figure of the current scene to a file
        """
        file_choices = " (*.png)|*.png| (*.pdf)|*.pdf |(*.jpg)|*.jpg |(*.eps)|*eps "
        filters=['.png','.pdf','.png','.png']

        
        dlg = wx.FileDialog(
            self, 
            message="Save figure to file...",
            defaultDir=os.getcwd(),
            defaultFile="",
            wildcard=file_choices,
            style= wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)

        if dlg.ShowModal() == wx.ID_OK:

            path = dlg.GetPath()
            ext = filters[dlg.GetFilterIndex()]
            if ext in path:
                outfile=path
            else:
                outfile = path+ext

            self.fig.savefig(outfile)

            


    def on_save_anim(self,event):
        """
        Save an animation of the current scene to a file
        """
        file_choices = "Quicktime (*.mov)|*.mov| (*.gif)|*.gif| (*.avi)|*.avi |(*.mp4)|*.mp4 "
        filters=['.mov','.gif','.avi','.mp4']

        
        dlg = wx.FileDialog(
            self, 
            message="Output animation file...",
            defaultDir=os.getcwd(),
            defaultFile="",
            wildcard=file_choices,
            style= wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)

        if dlg.ShowModal() == wx.ID_OK:

            path = dlg.GetPath()
            ext = filters[dlg.GetFilterIndex()]
            if ext in path:
                outfile=path
            else:
                outfile = path+ext
            self.flash_status_message("Saving figure to file: %s" %outfile)
            self.flash_status_message("Saving animation to file: %s" %outfile)

            # Create the animation
            #self.tstep = range(self.Nt) # Use all time steps for animation
            #self.animate(cbar=self.cbar,cmap=self.cmap,\
            #    xlims=self.axes.get_xlim(),ylims=self.axes.get_ylim())
            def updateScalar(i):
                self.tstep=[i]
                self.loadData()
                self.update_figure()
                return (self.title,self.collection)

            self.anim = animation.FuncAnimation(self.fig, updateScalar, frames=self.Nt, interval=50, blit=True)

            if ext=='.gif':
                self.anim.save(outfile,writer='imagemagick',fps=6)
            else:
                self.anim.save(outfile,fps=6,bitrate=3600)

            # Return the figure back to its status
            del self.anim
            self.tstep=self.tindex
            self.loadData()
            self.update_figure()

            # Bring up a dialog box
            dlg2= wx.MessageDialog(self, 'Animation complete.', "Done", wx.OK)
            dlg2.ShowModal()
            dlg2.Destroy()

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
   
    def on_plot_gridstat(self, event):
        """
        Plot the grid size histogram in a new figure
        """
        matplotlib.pyplot.figure()
        self.plothist()
        matplotlib.pyplot.show()


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

