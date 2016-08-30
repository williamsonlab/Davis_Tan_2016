"""
@author: Joey Davis <joeydavis@gmail.com>
"""

import os
import wx
import qMS
import qMSDefs
import string
import csv
import pandas as pd
import numpy
import sys
import matplotlib.gridspec as gridspec
import vizLib
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar
import mrmTools
import pylab

LIGHTSTRING = 'light '
HEAVYSTRING = 'heavy '
COLORMAP = pylab.cm.gist_rainbow
TITLE = 'Masse MRM'
ZOOMMAX = 5
FILENAMEHEADER = 'File Name'
class MasseFrame(wx.Frame):
    """ The main frame of the application
    """
    def __init__(self, df, fn, fsize=1.0, size=1.0):
        wx.Frame.__init__(self, None, -1, TITLE)                
        
        self.create_main_panel(df, fn, fsize=fsize, size=size)
        
        self.lowCheckNum.SetValue(False)
        self.highCheckNum.SetValue(True)
        self.lowCheckDen.SetValue(True)
        self.highCheckDen.SetValue(True)
        self.PPMtranLHRangeBypass.SetValue('10')
        self.PPMtranTRANALLRangeBypass.SetValue('15')
        self.PPMdsLRangeBypass.SetValue('-25 25')
        self.PPMdsHRangeBypass.SetValue('-25 25')       
        self.RTdsLHRangeBypass.SetValue('0.5')
        self.RTdsTRANALLRangeBypass.SetValue('0.5')
        self.RTpepTRANALLRangeBypass.SetValue('5')
        self.proteinZoomRangeBypass.SetValue('1 30')
        
        self.minIntensityBypass.SetValue('100 100')
        self.ratioLimBypass.SetValue('0.01 100')
        self.idpRangeBypass.SetValue('0.5')
        self.ldpRangeBypass.SetValue('0.5')
        self.lisSaveRangeBypass.SetValue('100 0')
        
        self.recalcAndDrawAll(setZero=True)

    def create_main_panel(self, df, fn, fsize=1.0, size=1.0):
        """ Creates the main panel with all the controls on it:
             * mpl canvas 
             * mpl navigation toolbar
             * Control panel for interaction
        """
        self.currentDirectory = os.getcwd()
        
        self.fsize = fsize
        self.size = size
        self.figdim = 6.5*fsize
        self.datafile = fn
        self.SetFont(wx.Font(10*size, wx.SWISS, wx.NORMAL, wx.NORMAL, False,'MS Shell Dlg 2'))
        
        self.calcNum = [LIGHTSTRING]
        self.calcDen = [LIGHTSTRING, HEAVYSTRING]

        self.currentHist = "PPMtranLH"
        self.dataFrame = df
        self.currentRow = self.dataFrame.iloc[0]
        self.currentISOFile = self.currentRow['UID']
        self.currentDataset = self.currentRow[FILENAMEHEADER]
        self.currentProtein = self.currentRow['Protein Name']
        self.proteinView = self.dataFrame[self.dataFrame['Protein Name'] == self.currentProtein]
        self.datasetView = self.dataFrame[self.dataFrame[FILENAMEHEADER] == self.currentDataset]

        self.positionLabelsDataset = qMS.sort_nicely(sorted(set(self.dataFrame['Protein Name'].values)))
        self.positionLabelsProtein = [i.split('.')[0] for i in qMS.sort_nicely(sorted(set(self.dataFrame[FILENAMEHEADER].values)))]

        self.panel = wx.Panel(self)
        '''Create the figure panels'''
        self.figLeft = Figure((self.figdim*1.1, self.figdim))
        self.canvasLeft = FigCanvas(self.panel, wx.ID_ANY, self.figLeft)
        self.figRight = Figure((self.figdim*1.1, self.figdim))
        self.canvasRight = FigCanvas(self.panel, wx.ID_ANY, self.figRight)
        
        gsLeft = gridspec.GridSpec(5,1)
        gsRight = gridspec.GridSpec(5,1)
        self.byDatasetPlot = self.figLeft.add_subplot(gsLeft[:4, :])
        self.byProteinPlot = self.figRight.add_subplot(gsRight[:4, :])
        
        self.histPlotDataset = self.figLeft.add_subplot(gsLeft[4, :])
        self.histPlotProtein = self.figRight.add_subplot(gsRight[4, :])
        self.figLeft.tight_layout()
        self.figRight.tight_layout()
        '''Create the list boxes'''
        self.datasetsList = wx.ListBox(self.panel, id=26, choices=[], style=wx.LB_SINGLE, name='Dataset', size=(100*size,100*size))        
        self.savedList = wx.ListBox(self.panel, id=26, choices=[], style=wx.LB_SINGLE, name='Saved', size=(100*size,100*size))
        self.filteredList = wx.ListBox(self.panel, id=26, choices=[], style=wx.LB_SINGLE, name='Filtered', size=(100*size,100*size))
        
        '''Create the buttons'''
        self.toolbarLeft = NavigationToolbar(self.canvasLeft)
        self.toolbarRight = NavigationToolbar(self.canvasRight)
        
        self.cb_grid = wx.CheckBox(self.panel, wx.ID_ANY, label="Grid", style=wx.ALIGN_RIGHT)        
        self.hideCheck = wx.CheckBox(self.panel, wx.ID_ANY, label="Hide", style=wx.ALIGN_RIGHT)        
        self.exportButton = wx.Button(self.panel, wx.ID_ANY, "Export")
        
        self.calcButton = wx.Button(self.panel, wx.ID_ANY, "Calc", size=(57*size,40))
        
        self.lowCheckNum = wx.CheckBox(self.panel, wx.ID_ANY, label="low")
        self.highCheckNum = wx.CheckBox(self.panel, wx.ID_ANY, label="high")
        
        self.lowCheckDen = wx.CheckBox(self.panel, wx.ID_ANY, label="low")
        self.highCheckDen = wx.CheckBox(self.panel, wx.ID_ANY, label="high")
        
        self.proteinZoomOn = wx.CheckBox(self.panel, wx.ID_ANY, label="pZoom", style=wx.ALIGN_RIGHT)
        self.zoomCheck = wx.CheckBox(self.panel, wx.ID_ANY, label="yZoom", style=wx.ALIGN_RIGHT)
        self.proteinZoomRangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)

        self.PPMtranLH_range_button = wx.Button(self.panel, wx.ID_ANY, "PPMLtoH",size=(75*size,-1))
        self.PPMtranTRANALL_range_button = wx.Button(self.panel, wx.ID_ANY, "PPMTtoP",size=(75*size,-1))
        self.PPMdsL_range_button = wx.Button(self.panel, wx.ID_ANY, "PPMdsL", size=(75*size,-1))
        self.PPMdsH_range_button = wx.Button(self.panel, wx.ID_ANY, "PPMdsH", size=(75*size,-1))
        self.RTdsLH_range_button = wx.Button(self.panel, wx.ID_ANY, "RTdsLH", size=(65*size,-1))
        self.RTdsTRANALL_range_button = wx.Button(self.panel, wx.ID_ANY, "RTdsTP", size=(65*size,-1))
        self.RTpepTRANALL_range_button = wx.Button(self.panel, wx.ID_ANY, "RTDrift", size=(65*size,-1))
        self.minIntensity_range_button = wx.Button(self.panel, wx.ID_ANY, "min inten.", size=(80*size,-1))
        self.ratio_range_button = wx.Button(self.panel, wx.ID_ANY, "amp[U/L]", size=(80*size,-1))
        self.idp_range_button = wx.Button(self.panel, wx.ID_ANY, "IDP", size=(80*size,-1))
        self.ldp_range_button = wx.Button(self.panel, wx.ID_ANY, "LDP", size=(80*size,-1))
        self.lis_range_button = wx.Button(self.panel, wx.ID_ANY, "LIS", size=(80*size,-1))
        
        self.PPMtranLHRangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.PPMtranTRANALLRangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.PPMdsLRangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.PPMdsHRangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.RTdsLHRangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.RTdsTRANALLRangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.RTpepTRANALLRangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.minIntensityBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.ratioLimBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.idpRangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.ldpRangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.lisSaveRangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        
        
        self.PPMtranLHOn = wx.CheckBox(self.panel, wx.ID_ANY, label="PPMLtoH")
        self.PPMtranTRANALLOn = wx.CheckBox(self.panel, wx.ID_ANY, label="PPMTtoP")
        self.PPMdsLOn = wx.CheckBox(self.panel, wx.ID_ANY, label="PPMdsL")
        self.PPMdsHOn = wx.CheckBox(self.panel, wx.ID_ANY, label="PPMdsH")
        self.RTdsLHOn = wx.CheckBox(self.panel, wx.ID_ANY, label="RTLtoH")
        self.RTdsTRANALLOn = wx.CheckBox(self.panel, wx.ID_ANY, label="RTTtoP")
        self.RTpepTRANALLOn = wx.CheckBox(self.panel, wx.ID_ANY, label="RTdrift")
        self.minIntensityOn = wx.CheckBox(self.panel, wx.ID_ANY, label="min amp")
        self.ratioLimOn = wx.CheckBox(self.panel, wx.ID_ANY, label="amp[U/L]")
        self.handSaveOn = wx.CheckBox(self.panel, wx.ID_ANY, label="save")
        self.handDeleteOn = wx.CheckBox(self.panel, wx.ID_ANY, label="delete")
        self.priorFilterOn = wx.CheckBox(self.panel, wx.ID_ANY, label="oldFilt")
        self.idpOn = wx.CheckBox(self.panel, wx.ID_ANY, label="IDP")
        self.ldpOn = wx.CheckBox(self.panel, wx.ID_ANY, label="LDP")
        self.lisOn = wx.CheckBox(self.panel, wx.ID_ANY, label="LISave")
        
        '''lay out the buttons'''
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.figbox = wx.BoxSizer(wx.HORIZONTAL)
        self.figbox.Add(self.canvasLeft, 0)
        self.listBox = wx.BoxSizer(wx.VERTICAL)
        self.datasetsListBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,self.datafile), wx.HORIZONTAL)
        self.datasetsListBox.Add(self.datasetsList, 1, flag=wx.GROW)    
        
        self.savedListBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'saved fits'), wx.HORIZONTAL)
        self.savedListBox.Add(self.savedList, 1, flag=wx.GROW)
        
        self.filteredListBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'filtered fits'), wx.HORIZONTAL)
        self.filteredListBox.Add(self.filteredList, 1, flag=wx.GROW)
        
        self.listBox.Add(self.datasetsListBox, 1, flag=wx.GROW)                
        self.listBox.Add(self.savedListBox, 1, flag=wx.GROW)        
        self.listBox.Add(self.filteredListBox, 1, flag=wx.GROW)
        self.figbox.Add(self.listBox, 1, flag=wx.GROW)
        self.figbox.Add(self.canvasRight, 0)
        self.vbox.Add(self.figbox, 0, flag = wx.GROW)

        #new horizontal box    
        self.toolNumBox = wx.BoxSizer(wx.HORIZONTAL)
        
        #add toolbar, draw, grid and exportto hbox
        self.imageToolsBoxL = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'image tools'), wx.HORIZONTAL)
        
        self.toolbarSubunits = wx.BoxSizer(wx.VERTICAL)
        self.toolbarSubunits.Add(self.toolbarLeft, 0, flag=wx.ALIGN_LEFT)
        self.imageToolsBoxL.Add(self.toolbarSubunits, 0, flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER)
        
        self.imageCheckToolsBox = wx.BoxSizer(wx.VERTICAL)        
        self.imageCheckToolsBox.Add(self.cb_grid, 0, flag=wx.ALIGN_LEFT)
        self.imageCheckToolsBox.Add(self.hideCheck, 0, flag=wx.ALIGN_LEFT)
        self.imageCheckToolsBox.Add(self.zoomCheck, 0, flag=wx.ALIGN_LEFT)
        self.imageToolsBoxL.Add(self.imageCheckToolsBox, 0, flag=wx.ALIGN_LEFT)
        self.toolNumBox.Add(self.imageToolsBoxL, 0, flag=wx.ALIGN_LEFT)

        #add open and export buttons
        self.fileToolsBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'file tools'), wx.VERTICAL)
        self.fileToolsBox.Add(self.calcButton, flag=wx.ALIGN_CENTER)
        self.fileToolsBox.Add(self.exportButton, 0, flag=wx.ALIGN_CENTER)       
        self.toolNumBox.Add(self.fileToolsBox, 0, flag=wx.ALIGN_LEFT)

        #add numerator box to hbox
        self.numBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'num'), wx.VERTICAL)
        self.numBox.Add(self.lowCheckNum, 0, flag=wx.ALIGN_LEFT)
        self.numBox.Add(self.highCheckNum, 0,flag=wx.ALIGN_LEFT)
        self.toolNumBox.Add(self.numBox, 0, flag=wx.ALIGN_LEFT)
        #add denominator box to hbox
        self.denBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'den'), wx.VERTICAL)
        self.denBox.Add(self.lowCheckDen, 0, flag=wx.ALIGN_LEFT)
        self.denBox.Add(self.highCheckDen, 0,flag=wx.ALIGN_LEFT)
        self.toolNumBox.Add(self.denBox, 0, flag=wx.ALIGN_LEFT)
        
        self.imageToolsBoxR = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'image tools'), wx.HORIZONTAL)
        self.imageToolsBoxR.Add(self.toolbarRight, 1, flag=wx.ALIGN_RIGHT)
        self.zoomBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'zoom tools'), wx.VERTICAL)
        self.zoomBox.Add(self.zoomCheck, 0, flag=wx.ALIGN_LEFT)
        self.zoomBox.Add(self.proteinZoomOn, 0, flag=wx.ALIGN_LEFT)
        self.zoomBox.Add(self.proteinZoomRangeBypass, 0, flag=wx.ALIGN_LEFT)
        
        self.filtersBox = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY, 'filters'), wx.HORIZONTAL)
        self.filts1 = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY), wx.VERTICAL)
        self.filts1.Add(self.PPMtranLHOn, 0, flag=wx.ALIGN_LEFT)
        self.filts1.Add(self.PPMtranTRANALLOn, 0, flag=wx.ALIGN_LEFT)
        self.filts1.Add(self.PPMdsLOn, 0, flag=wx.ALIGN_LEFT)
        self.filts2 = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY), wx.VERTICAL)
        self.filts2.Add(self.PPMdsHOn, 0, flag=wx.ALIGN_LEFT)
        self.filts2.Add(self.RTdsLHOn, 0, flag=wx.ALIGN_LEFT)
        self.filts2.Add(self.RTdsTRANALLOn, 0, flag=wx.ALIGN_LEFT)
        self.filts3 = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY), wx.VERTICAL)
        self.filts3.Add(self.RTpepTRANALLOn, 0, flag=wx.ALIGN_LEFT)
        self.filts3.Add(self.minIntensityOn, 0, flag=wx.ALIGN_LEFT)
        self.filts3.Add(self.handSaveOn, 0, flag=wx.ALIGN_LEFT)
        self.filts4 = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY), wx.VERTICAL)
        self.filts4.Add(self.handDeleteOn, 0, flag=wx.ALIGN_LEFT)
        self.filts4.Add(self.ratioLimOn, 0, flag=wx.ALIGN_LEFT)
        self.filts4.Add(self.idpOn, 0, flag=wx.ALIGN_LEFT)
        self.filts5 = wx.StaticBoxSizer(wx.StaticBox(self.panel, wx.ID_ANY), wx.VERTICAL)
        self.filts5.Add(self.ldpOn, 0, flag=wx.ALIGN_LEFT)
        self.filts5.Add(self.priorFilterOn, 0, flag=wx.ALIGN_LEFT)
        self.filts5.Add(self.lisOn, 0, flag=wx.ALIGN_LEFT)

        self.filtersBox.Add(self.filts1, 0, flag=wx.ALIGN_LEFT)
        self.filtersBox.Add(self.filts2, 0, flag=wx.ALIGN_LEFT)
        self.filtersBox.Add(self.filts3, 0, flag=wx.ALIGN_LEFT)
        self.filtersBox.Add(self.filts4, 0, flag=wx.ALIGN_LEFT)
        self.filtersBox.Add(self.filts5, 0, flag=wx.ALIGN_LEFT)
    
        self.toolNumBox.Add(self.zoomBox, 0, flag=wx.ALIGN_LEFT)
        self.toolNumBox.Add(self.filtersBox, 0, flag=wx.ALIGN_LEFT)
        self.expandBox = wx.BoxSizer(wx.HORIZONTAL)
        self.toolNumBox.Add(self.expandBox, 1, flag=wx.GROW)

        self.toolNumBox.Add(self.imageToolsBoxR, 0, flag=wx.ALIGN_RIGHT)
        
        #add hbox to vbox
        self.vbox.Add(self.toolNumBox, 0, flag = wx.GROW)
        
        # Sliders for setting the various cutoffs
        self.filterBox1 = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'filter cuttoff settings'), wx.HORIZONTAL)
        self.filterBox2 = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY), wx.HORIZONTAL)
        flags = wx.ALIGN_LEFT | wx.ALIGN_CENTER
        self.filterBox1.Add(self.PPMtranLH_range_button, 0, flag=flags)
        self.filterBox1.Add(self.PPMtranLHRangeBypass, 0, flag=flags)
        self.filterBox1.Add(self.PPMtranTRANALL_range_button, 0, flag=flags)
        self.filterBox1.Add(self.PPMtranTRANALLRangeBypass, 0, flag=flags)
        self.filterBox1.Add(self.PPMdsL_range_button, 0, flag=flags)
        self.filterBox1.Add(self.PPMdsLRangeBypass, 0, flag=flags)
        self.filterBox1.Add(self.PPMdsH_range_button, 0, flag=flags)
        self.filterBox1.Add(self.PPMdsHRangeBypass, 0, flag = flags)
        self.filterBox1.Add(self.RTdsLH_range_button, 0, flag=flags)
        self.filterBox1.Add(self.RTdsLHRangeBypass, 0, flag = flags)
        self.filterBox2.Add(self.RTdsTRANALL_range_button, 0, flag=flags)
        self.filterBox2.Add(self.RTdsTRANALLRangeBypass, 0, flag = flags)
        self.filterBox2.Add(self.RTpepTRANALL_range_button, 0, flag=flags)
        self.filterBox2.Add(self.RTpepTRANALLRangeBypass, 0, flag = flags)
        self.filterBox2.Add(self.minIntensity_range_button, 0, flag=flags)
        self.filterBox2.Add(self.minIntensityBypass, 0, flag = flags)
        self.filterBox2.Add(self.ratio_range_button, 0, flag=flags)
        self.filterBox2.Add(self.ratioLimBypass, 0, flag = flags)
        self.filterBox2.Add(self.idp_range_button, 0, flag=flags)
        self.filterBox2.Add(self.idpRangeBypass, 0, flag = flags)
        self.filterBox2.Add(self.ldp_range_button, 0, flag=flags)
        self.filterBox2.Add(self.ldpRangeBypass, 0, flag = flags)
        self.filterBox2.Add(self.lis_range_button, 0, flag=flags)
        self.filterBox2.Add(self.lisSaveRangeBypass, 0, flag = flags)
        self.vbox.Add(self.filterBox1, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        self.vbox.Add(self.filterBox2, 0, flag = wx.ALIGN_LEFT | wx.TOP)

        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)
        
        '''bind events for buttons'''
        self.canvasLeft.mpl_connect('pick_event', self.pickScatterPoint)
        self.canvasRight.mpl_connect('pick_event', self.pickScatterPointProtein)

        self.exportButton.Bind(wx.EVT_BUTTON, self.on_export_button)
        self.calcButton.Bind(wx.EVT_BUTTON, self.recalcAndDrawAll)
        self.cb_grid.Bind(wx.EVT_CHECKBOX, self.redrawAll)
        self.hideCheck.Bind(wx.EVT_CHECKBOX, self.redrawAll)
        self.zoomCheck.Bind(wx.EVT_CHECKBOX, self.redrawAll)

        self.PPMtranLH_range_button.Bind(wx.EVT_BUTTON, self.on_PPMtranLH_range_button)
        self.PPMtranTRANALL_range_button.Bind(wx.EVT_BUTTON, self.on_PPMtranTRANALL_range_button)
        self.PPMdsL_range_button.Bind(wx.EVT_BUTTON, self.on_PPMdsL_range_button)
        self.PPMdsH_range_button.Bind(wx.EVT_BUTTON, self.on_PPMdsH_range_button)
        self.RTdsLH_range_button.Bind(wx.EVT_BUTTON, self.on_RTdsLH_range_button)
        self.RTdsTRANALL_range_button.Bind(wx.EVT_BUTTON, self.on_RTdsTRANALL_range_button)
        self.RTpepTRANALL_range_button.Bind(wx.EVT_BUTTON, self.on_RTpepTRANALL_range_button)
        self.minIntensity_range_button.Bind(wx.EVT_BUTTON, self.on_minIntensity_range_button)
        self.idp_range_button.Bind(wx.EVT_BUTTON, self.on_idp_range_button)
        self.ldp_range_button.Bind(wx.EVT_BUTTON, self.on_ldp_range_button)
        self.ratio_range_button.Bind(wx.EVT_BUTTON, self.on_ratio_range_button)
        
        self.PPMtranLHOn.Bind(wx.EVT_CHECKBOX, self.on_checkbox)
        self.PPMtranTRANALLOn.Bind(wx.EVT_CHECKBOX, self.on_checkbox)
        self.PPMdsLOn.Bind(wx.EVT_CHECKBOX, self.on_checkbox)
        self.PPMdsHOn.Bind(wx.EVT_CHECKBOX, self.on_checkbox)
        self.RTdsLHOn.Bind(wx.EVT_CHECKBOX, self.on_checkbox)
        self.RTdsTRANALLOn.Bind(wx.EVT_CHECKBOX, self.on_checkbox)
        self.RTpepTRANALLOn.Bind(wx.EVT_CHECKBOX, self.on_checkbox)
        self.minIntensityOn.Bind(wx.EVT_CHECKBOX, self.on_checkbox)
        self.handSaveOn.Bind(wx.EVT_CHECKBOX, self.on_checkbox)
        self.handDeleteOn.Bind(wx.EVT_CHECKBOX, self.on_checkbox)
        self.priorFilterOn.Bind(wx.EVT_CHECKBOX, self.on_checkbox)
        self.idpOn.Bind(wx.EVT_CHECKBOX, self.on_checkbox)
        self.ldpOn.Bind(wx.EVT_CHECKBOX, self.on_checkbox)
        self.lisOn.Bind(wx.EVT_CHECKBOX, self.on_checkbox)
        self.ratioLimOn.Bind(wx.EVT_CHECKBOX, self.on_checkbox)
        
        self.datasetsList.Bind(wx.EVT_LISTBOX, self.on_datasetsBoxClick)
        self.savedList.Bind(wx.EVT_LISTBOX, self.on_savedBoxClick)
        self.filteredList.Bind(wx.EVT_LISTBOX, self.on_filteredBoxClick)
        
        self.panel.Bind(wx.EVT_KEY_UP, self.on_key_press)

    def on_savedBoxClick(self, event):
        if not (self.savedList.GetStringSelection() is u''):
            self.currentISOFile = self.savedList.GetStringSelection()
            self.savedList.SetStringSelection(self.currentISOFile)
            self.filteredList.SetStringSelection('')
            self.newPepSelection()
        
    def on_filteredBoxClick(self, event):
        if not (self.filteredList.GetStringSelection() is u''):
            self.currentISOFile = self.filteredList.GetStringSelection()
            self.filteredList.SetStringSelection(self.currentISOFile)
            self.savedList.SetStringSelection('')
            self.newPepSelection()
    
    def on_datasetsBoxClick(self, event):
        if not (self.datasetsList.GetStringSelection() is u''):
            self.currentDataset = self.datasetsList.GetStringSelection()
            self.datasetView = self.dataFrame[self.dataFrame[FILENAMEHEADER] == self.currentDataset]    
            self.currentRow = self.datasetView.iloc[0]
            self.currentISOFile = self.currentRow['UID']
            self.currentProtein = self.currentRow['Protein Name']
            self.calc_hist()
            self.redrawAll()
   
    def pickScatterPoint(self, event):
        ind = event.ind
        thisline = event.artist
        ydata = thisline.get_ydata()
        xdata = thisline.get_xdata()
        self.currentRow = self.datasetView[(self.datasetView['currentCalc'] == ydata[ind][0]) & \
                                                ((self.datasetView['currentPosDataset']+self.datasetView['colOff']) == xdata[ind][0])].iloc[0]
        self.currentISOFile = self.currentRow['UID']
        self.currentProtein = self.currentRow['Protein Name']
        self.proteinView = self.dataFrame[self.dataFrame['Protein Name'] == self.currentProtein]
        
        self.selectedPointDataset.set_data(self.currentRow['currentPosDataset']+self.currentRow['colOff'], 
                                           self.currentRow['currentCalc'])
        
        self.savedList.SetStringSelection(self.currentISOFile)
        self.filteredList.SetStringSelection(self.currentISOFile)
        
        self.updateTextDataset()
        self.calc_hist()
        self.canvasLeft.draw()
        self.calc_figureRight()
        self.canvasRight.draw()
    
    def pickScatterPointProtein(self, event):
        ind = event.ind
        thisline = event.artist
        ydata = thisline.get_ydata()
        xdata = thisline.get_xdata()
        self.currentRow = self.proteinView[(self.proteinView['currentCalc'] == ydata[ind][0]) & \
                                                ((self.proteinView['currentPosProtein']+self.proteinView['colOff']) == xdata[ind][0])].iloc[0]
        self.currentISOFile = self.currentRow['UID']
        self.currentDataset = self.currentRow[FILENAMEHEADER]
        
        theRows = self.proteinView[self.proteinView['TID'] == self.currentRow['TID']]
        self.selectedPointsProtein.set_data(theRows['currentPosProtein']+theRows['colOff'], theRows['currentCalc'])

        self.updateTextProtein()

        self.calc_figureLeft()
        self.calc_lists()
        self.calc_hist()
        self.canvasLeft.draw()
        self.canvasRight.draw()        
        
    def newPepSelection(self):
        self.currentRow = self.dataFrame[self.dataFrame['UID'] == self.currentISOFile].iloc[0]
        self.currentISOFile = self.currentRow['UID']
        self.currentProtein = self.currentRow['Protein Name']
        self.proteinView = self.dataFrame[self.dataFrame['Protein Name'] == self.currentProtein]
        self.updateTextDataset()
        self.calc_hist()
        self.calc_figureRight()
        self.selectedPointDataset.set_data(self.currentRow['currentPosDataset']+self.currentRow['colOff'], 
                                           self.currentRow['currentCalc'])
        self.canvasLeft.draw()
        self.canvasRight.draw()
        
    def on_key_press(self, event):
        event.Skip()
        c = event.GetKeyCode()
        if c is 68: #got a d keystroke
            self.dataFrame.loc[self.dataFrame['UID'] == self.currentISOFile,'handDelete'] = True
            self.dataFrame.loc[self.dataFrame['UID'] == self.currentISOFile,'handSave'] = False
            try:            
                myIndex = list(self.savedListItems).index(self.currentISOFile)
            except ValueError:
                myIndex = 0
            try:
                nextItem = self.savedListItems[myIndex+1]
            except IndexError:
                nextItem = self.savedListItems[myIndex-1]
            self.currentISOFile = nextItem
            self.recalcAndDrawAll()
            
        elif c is 83: # got a s keystroke
            self.dataFrame.loc[self.dataFrame['UID'] == self.currentISOFile,'handDelete'] = False
            self.dataFrame.loc[self.dataFrame['UID'] == self.currentISOFile,'handSave'] = True
            self.recalcAndDrawAll()
            
        elif c is 69:
            self.dataFrame.loc[self.dataFrame['TID'] == self.currentRow['TID'], 'handDelete'] = True
            self.dataFrame.loc[self.dataFrame['TID'] == self.currentRow['TID'], 'handSave'] = False
            try:
                myIndex = list(self.savedListItems).index(self.currentISOFile)
            except ValueError:
                myIndex = 0
            try:
                nextItem = self.savedListItems[myIndex+1]
            except IndexError:
                nextItem = self.savedListItems[myIndex-1]
            self.currentISOFile = nextItem
            self.recalcAndDrawAll()
            
        elif c is 87:
            self.dataFrame.loc[self.dataFrame['TID'] == self.currentRow['TID'], 'handDelete'] = False
            self.dataFrame.loc[self.dataFrame['TID'] == self.currentRow['TID'], 'handSave'] = True
            self.recalcAndDrawAll()

    def calc_lists(self):
        self.savedListItems = numpy.array(qMS.sort_nicely(list(self.datasetView[self.datasetView['allClear'] == True]['UID'].values)))
        self.filteredListItems = numpy.array(qMS.sort_nicely(list(self.datasetView[self.datasetView['allClear'] == False]['UID'].values)))
        self.savedList.Set(self.savedListItems)
        self.filteredList.Set(self.filteredListItems)
        self.savedList.SetStringSelection(self.currentISOFile)
        self.filteredList.SetStringSelection(self.currentISOFile)
        
        self.datasetsListItems = qMS.sort_nicely(list(self.dataFrame[FILENAMEHEADER].unique()))
        self.datasetsList.Set(self.datasetsListItems)
        self.datasetsList.SetStringSelection(self.currentDataset)

    def recalcAndDrawAll(self, setZero=False):
        self.calc_data()
        self.calc_hist()
        self.calc_figureLeft()
        self.calc_figureRight()
        self.calc_lists()
        self.draw_all()
    def redrawAll(self, event=None):
        self.calc_figureLeft()
        self.calc_figureRight()
        self.calc_lists()
        self.draw_all()

    def calc_data(self):
        self.calcNum = self.getChecksNum()
        self.calcDen = self.getChecksDen()
        self.dataFrame['currentCalc'] = mrmTools.calcValue(self.dataFrame, self.calcNum, self.calcDen)
        self.filterData()
        self.currentRow = self.dataFrame[self.dataFrame['UID'] == self.currentISOFile].iloc[0]
        self.currentISOFile = self.currentRow['UID']
        self.currentProtein = self.currentRow['Protein Name']        
        self.currentDataset = self.currentRow[FILENAMEHEADER]
        self.datasetView = self.dataFrame[self.dataFrame[FILENAMEHEADER] == self.currentDataset]
        self.proteinView = self.dataFrame[self.dataFrame['Protein Name'] == self.currentRow['Protein Name']]
               
    def filterData(self):
        self.findRanges()
        self.dataFrame['allClear'] = self.getPass() | self.dataFrame['handSave']
    
    def findRanges(self):    
        self.PPMtranLH = float(self.PPMtranLHRangeBypass.GetValue())
        self.PPMtranTRANALL = float(self.PPMtranTRANALLRangeBypass.GetValue())
        (self.PPMdsL_low, self.PPMdsL_high) = map(float, self.PPMdsLRangeBypass.GetValue().split(' '))
        (self.PPMdsH_low, self.PPMdsH_high) = map(float, self.PPMdsHRangeBypass.GetValue().split(' '))
        self.RTdsLH = float(self.RTdsLHRangeBypass.GetValue())
        self.RTdsTRANALL = float(self.RTdsTRANALLRangeBypass.GetValue())
        self.RTpepTRANALL = float(self.RTpepTRANALLRangeBypass.GetValue())
        (self.ratio_low, self.ratio_high) = map(float, self.ratioLimBypass.GetValue().split(' '))
        self.minIntensity = map(float, self.minIntensityBypass.GetValue().split(' '))
        self.idp = float(self.idpRangeBypass.GetValue())
        self.ldp = float(self.ldpRangeBypass.GetValue())
        (self.lisMinIntensityL, self.lisMinIntensityH) = map(float, self.lisSaveRangeBypass.GetValue().split(' '))

    def getPass(self):
        self.findRanges()
        filt = self.dataFrame['Missed Cleavages'] >= -1
        if self.priorFilterOn.IsChecked():
            filt = filt & (self.dataFrame['priorFilter'] == True)
        if self.PPMtranLHOn.IsChecked():
            filt =  filt & (self.dataFrame['PPMtranLH'] <= self.PPMtranLH)
        if self.PPMtranTRANALLOn.IsChecked():
            filt =  filt & (self.dataFrame['PPMtranTRANALL'] <= self.PPMtranTRANALL)
        if self.PPMdsLOn.IsChecked():
            filt = filt & (self.dataFrame[LIGHTSTRING+'Mass Error PPM'] >= self.PPMdsL_low) & (self.dataFrame[LIGHTSTRING+'Mass Error PPM'] <= self.PPMdsL_high)
        if self.PPMdsHOn.IsChecked():
            filt = filt & (self.dataFrame[HEAVYSTRING+'Mass Error PPM'] >= self.PPMdsH_low) & (self.dataFrame[HEAVYSTRING+'Mass Error PPM'] <= self.PPMdsH_high)
        if self.RTdsLHOn.IsChecked():
            filt = filt & (self.dataFrame['RTdsLH'] <= self.RTdsLH)
        if self.RTdsTRANALLOn.IsChecked():
            filt = filt & (self.dataFrame['RTdsTRANALL'] <= self.RTdsTRANALL)
        if self.RTpepTRANALLOn.IsChecked():
           filt = filt & (self.dataFrame['RTpepTRANALL'] <= self.RTpepTRANALL)
        if self.ratioLimOn.IsChecked():
            filt = filt & (self.dataFrame['ratio'] >= self.ratio_low) & (self.dataFrame['ratio'] <= self.ratio_high)
        if self.minIntensityOn.IsChecked():
            filt = filt & (self.dataFrame[LIGHTSTRING+'AdjArea'] >= self.minIntensity[0]) & (self.dataFrame[HEAVYSTRING+'AdjArea'] >= self.minIntensity[1])
        if self.idpOn.IsChecked():
            filt = filt & (self.dataFrame[LIGHTSTRING+'Isotope Dot Product'] >= self.idp) & (self.dataFrame[HEAVYSTRING+'Isotope Dot Product'] >= self.idp)
        if self.ldpOn.IsChecked():
            filt = filt & (self.dataFrame[LIGHTSTRING+'Library Dot Product'] >= self.ldp) & (self.dataFrame[HEAVYSTRING+'Library Dot Product'] >= self.ldp)
        if self.lisOn.IsChecked():
            filt = filt | (self.dataFrame[LIGHTSTRING+'AdjArea'] <= self.lisMinIntensityL) | (self.dataFrame[HEAVYSTRING+'AdjArea'] <= self.lisMinIntensityH)
        if self.handSaveOn.IsChecked():
            filt = filt | (self.dataFrame['handSave'] == True)
        if self.handDeleteOn.IsChecked():
            filt = filt & (self.dataFrame['handDelete'] == False)
        return filt
        
    def calc_figureLeft(self):
        s = ':'
        self.byDatasetPlot.clear()
        self.datasetView = self.dataFrame[self.dataFrame[FILENAMEHEADER] == self.currentDataset]
        self.passDatasetView = self.datasetView[self.datasetView['allClear'] == True]
        self.failDatasetView = self.datasetView[self.datasetView['allClear'] == False]
        
        for c in list(self.passDatasetView['colOff'].unique()):
            color = COLORMAP(c)
            toPlot = self.passDatasetView[self.passDatasetView['colOff'] == c]
            toPlotPrecursors = toPlot[toPlot['Fragment Ion Type'] == 'precursor']
            self.byDatasetPlot.plot(toPlot['currentPosDataset']+toPlot['colOff'], toPlot['currentCalc'], 'o', 
                                       color=color, picker=5, mec='none', alpha=1, ms=5*self.size)
            self.byDatasetPlot.plot(toPlotPrecursors['currentPosDataset']+toPlotPrecursors['colOff'], toPlotPrecursors['currentCalc'], 'o', 
                                       mfc='none', picker=5, mec='black', alpha=1, ms=5*self.size)
            
        self.byDatasetPlot.plot([-1], [-1], 'o', color='red', label='Saved: ' + str(self.passDatasetView.shape[0]), visible=False)
        self.byDatasetPlot.grid(self.cb_grid.IsChecked())
        self.selectedPointDataset, = self.byDatasetPlot.plot(self.currentRow['currentPosDataset'] + self.currentRow['colOff'], self.currentRow['currentCalc'], 
                                              'o', ms=30*self.size, alpha=0.5, color='yellow')
        
        meds = self.determineMedians(self.passDatasetView, 'currentPosDataset')
        self.byDatasetPlot.plot(meds[0], meds[1], 'g-', lw=2, label="Median : " + str(round(numpy.median(meds[1]),1)))
        
        if not self.hideCheck.IsChecked():
            mec = 'grey'
            mark = 'x'
            self.byDatasetPlot.plot(self.failDatasetView['currentPosDataset']+self.failDatasetView['colOff'], self.failDatasetView['currentCalc'], mark, 
                                    mec=mec, picker=5, label="Filtered : " + str(self.failDatasetView.shape[0]), color = 'grey')

        self.byDatasetPlot.set_xticks(range(1,int(self.datasetView['currentPosDataset'].max())+1))
        self.byDatasetPlot.set_xticklabels(self.positionLabelsDataset, rotation=90, size='small')
        if self.proteinZoomOn.IsChecked():
            (x_low, x_high) = map(float, self.proteinZoomRangeBypass.GetValue().split(' '))
            self.byDatasetPlot.set_xlim([x_low,x_high])
        else:
            self.byDatasetPlot.set_xlim([0.5,self.datasetView['currentPosDataset'].max()+1.5])

        self.byDatasetPlot.set_title(self.currentDataset+ " : " + setCurrentFrac(self.calcNum, self.calcDen))
        if not self.zoomCheck.IsChecked():
            self.byDatasetPlot.set_ylim([0,max(self.passDatasetView['currentCalc'].max(),1)])
        else:
            self.byDatasetPlot.set_ylim([0,ZOOMMAX])
            #self.byDatasetPlot.set_xlim([0,10])
        self.byDatasetPlot.legend()
        self.calc_tran_string()
        self.byDatasetPlot.text(0.02, 0.98,self.tranLabelString, horizontalalignment='left',verticalalignment='top',
                                transform = self.byDatasetPlot.transAxes,color = 'black', size=9*self.fsize)
        self.lTextPep = self.byDatasetPlot.text(0.125, 0.98,self.tranLightColumn,horizontalalignment='left',verticalalignment='top',
                                transform = self.byDatasetPlot.transAxes,color = '#0000FF', size=9*self.fsize)
        self.hTextPep = self.byDatasetPlot.text(0.225, 0.98,self.tranHeavyColumn,horizontalalignment='left',verticalalignment='top',
                                transform = self.byDatasetPlot.transAxes,color = '#FF0000', size=9*self.fsize)
        self.dsTextPep = self.byDatasetPlot.text(0.325, 0.98,self.tranDSColumn,horizontalalignment='left',verticalalignment='top',
                                transform = self.byDatasetPlot.transAxes,color = '#000000', size=9*self.fsize)
        vizLib.cleanAxis(self.byDatasetPlot, ticksOnly=True)

    def calc_figureRight(self):
        self.byProteinPlot.clear()
        self.proteinView = self.dataFrame[self.dataFrame['Protein Name'] == self.currentProtein]
        self.passProteinView = self.proteinView[self.proteinView['allClear'] == True]
        self.failProteinView = self.proteinView[self.proteinView['allClear'] == False]
        
        for c in list(self.passProteinView['colOff'].unique()):
            color = COLORMAP(c)
            toPlot = self.passProteinView[self.passProteinView['colOff'] == c]
            toPlotPrecursors = toPlot[toPlot['Fragment Ion Type'] == 'precursor']

            self.byProteinPlot.plot(toPlot['currentPosProtein']+toPlot['colOff'], toPlot['currentCalc'], 'o', 
                                       color=color, picker=5, mec='none', alpha=1, ms=5*self.size)
            self.byProteinPlot.plot(toPlotPrecursors['currentPosProtein']+toPlotPrecursors['colOff'], toPlotPrecursors['currentCalc'], 'o', 
                                       mfc='none', picker=5, mec='black', alpha=1, ms=5*self.size)            


        theRows = self.proteinView[self.proteinView['TID'] == self.currentRow['TID']]
        self.selectedPointsProtein, = self.byProteinPlot.plot(theRows['currentPosProtein'] + theRows['colOff'], theRows['currentCalc'], 
                                              'o', ms=30*self.size, alpha=0.25, color='red')
        self.byProteinPlot.plot([-1], [-1], 'o', color='green', label='Saved: ' + str(self.passProteinView.shape[0]), visible=False)
        
        self.byProteinPlot.grid(self.cb_grid.IsChecked())

        meds = self.determineMedians(self.passProteinView, 'currentPosProtein')
        self.byProteinPlot.plot(meds[0], meds[1], 'g-', lw=2, label="Median : " + str(round(numpy.median(meds[1]),1)))
        
        if not (self.hideCheck.IsChecked()):
            mec = 'grey'
            mark = 'x'
            self.byProteinPlot.plot(self.failProteinView['currentPosProtein']+self.failProteinView['colOff'], self.failProteinView['currentCalc'], mark, 
                                    mec=mec, picker=5, label="Filtered : " + str(self.failProteinView.shape[0]), color = 'grey')
                                    
        self.byProteinPlot.set_xticks(range(1,int(self.proteinView['currentPosProtein'].max())+1))
        self.byProteinPlot.set_xticklabels(self.positionLabelsProtein, rotation=90, size='small')
        self.byProteinPlot.set_title(self.currentProtein+" : "+self.currentRow['Peptide Modified Sequence']+' : '+self.currentRow['TID']+\
                                    ':'+self.currentRow['Library Name'][0:8]+':'+str(self.currentRow['Library Score1']))
        self.byProteinPlot.legend()
        self.byProteinPlot.set_xlim(0, int(self.proteinView['currentPosProtein'].max())+2)
        if not self.zoomCheck.IsChecked():
            self.byProteinPlot.set_ylim(0, max(self.passProteinView['currentCalc'].max(),1))
        else:
            self.byProteinPlot.set_ylim([0, ZOOMMAX])
            #self.byDatasetPlot.set_xlim([0,10])
        
        
        self.calc_prot_string()
        self.byProteinPlot.text(0.02, 0.98,self.protLabelString, horizontalalignment='left',verticalalignment='top',
                                transform = self.byProteinPlot.transAxes,color = 'black', size=9*self.fsize)
        self.lTextProt = self.byProteinPlot.text(0.125, 0.98,self.pepLightColumn,horizontalalignment='left',verticalalignment='top',
                                transform = self.byProteinPlot.transAxes,color = '#0000FF', size=9*self.fsize)
        self.hTextProt = self.byProteinPlot.text(0.225, 0.98,self.pepHeavyColumn,horizontalalignment='left',verticalalignment='top',
                                transform = self.byProteinPlot.transAxes,color = '#FF0000', size=9*self.fsize)              
        vizLib.cleanAxis(self.byProteinPlot, ticksOnly=True)

    def draw_all(self):
        self.canvasLeft.draw()
        self.canvasRight.draw()   
    
    def determineMedians(self, view, field):
        xs = list(set(view[field]))
        ys = [view[view[field]==k]['currentCalc'].median() for k in xs]
        xs = [float(i)+0.5 for i in xs]
        return [xs,ys]
        
    def getChecksNum(self):
        toReturn = []        
        if self.lowCheckNum.IsChecked():
            toReturn.append(LIGHTSTRING)
        if self.highCheckNum.IsChecked():
            toReturn.append(HEAVYSTRING)
        return toReturn
        
    def getChecksDen(self):
        toReturn = []        
        if self.lowCheckDen.IsChecked():
            toReturn.append(LIGHTSTRING)
        if self.highCheckDen.IsChecked():
            toReturn.append(HEAVYSTRING)
        return toReturn
    
    def updateTextDataset(self):
        self.calc_tran_string()
        self.lTextPep.set_text(self.tranLightColumn)
        self.hTextPep.set_text(self.tranHeavyColumn)
        self.dsTextPep.set_text(self.tranDSColumn)

    def updateTextProtein(self):
        self.calc_prot_string()
        self.lTextProt.set_text(self.pepLightColumn)
        self.hTextProt.set_text(self.pepHeavyColumn)

    def calc_tran_string(self):
        lPPM = self.currentRow[LIGHTSTRING+'Mass Error PPM']
        tlPPM = self.currentRow[LIGHTSTRING+'Average Mass Error PPM']
        lRT = round(self.currentRow[LIGHTSTRING+'Retention Time'],1)
        tlRT = round(self.currentRow[LIGHTSTRING+'Best Retention Time'],1)
        lFWHM = self.currentRow[LIGHTSTRING+'Fwhm']
        tlFWHM = round(self.passDatasetView[self.passDatasetView['PID'] == self.currentRow['PID']][LIGHTSTRING+'Fwhm'].median(),2)
        lHeight = self.currentRow[LIGHTSTRING+'Height']
        tlHeight = self.currentRow[LIGHTSTRING+'Max Height']
        lArea = self.currentRow[LIGHTSTRING+'Area']
        lBack = self.currentRow[LIGHTSTRING+'Background']
        tlAreaAdj = self.currentRow[LIGHTSTRING+'Total Area'] - self.currentRow[LIGHTSTRING+'Total Background']
        lMS2MS1 = round(self.currentRow[LIGHTSTRING+'Total Area Fragment']/self.currentRow[LIGHTSTRING+'Total Area MS1'],2)
        lIDP = round(self.currentRow[LIGHTSTRING+'Isotope Dot Product'],2)
        lLDP = round(self.currentRow[LIGHTSTRING+'Library Dot Product'],2)
        
        hPPM = self.currentRow[HEAVYSTRING+'Mass Error PPM']
        thPPM = self.currentRow[HEAVYSTRING+'Average Mass Error PPM']
        hRT = round(self.currentRow[HEAVYSTRING+'Retention Time'],1)
        thRT = round(self.currentRow[HEAVYSTRING+'Best Retention Time'],1)
        hFWHM = self.currentRow[HEAVYSTRING+'Fwhm']
        thFWHM = round(self.passDatasetView[self.passDatasetView['PID'] == self.currentRow['PID']][HEAVYSTRING+'Fwhm'].median(),2)
        hHeight = self.currentRow[HEAVYSTRING+'Height']
        thHeight = self.currentRow[HEAVYSTRING+'Max Height']
        hArea = self.currentRow[HEAVYSTRING+'Area']
        hBack = self.currentRow[HEAVYSTRING+'Background']
        thAreaAdj = self.currentRow[HEAVYSTRING+'Total Area'] - self.currentRow[HEAVYSTRING+'Total Background']
        hMS2MS1 = round(self.currentRow[HEAVYSTRING+'Total Area Fragment']/self.currentRow[HEAVYSTRING+'Total Area MS1'],2)
        hIDP = round(self.currentRow[HEAVYSTRING+'Isotope Dot Product'],2)
        hLDP = round(self.currentRow[HEAVYSTRING+'Library Dot Product'],2)
        hlDP = round(self.currentRow[HEAVYSTRING+'Ratio Dot Product'],2)
        
        ldsPPM = self.passDatasetView[LIGHTSTRING + 'Mass Error PPM'].median()
        hdsPPM = self.passDatasetView[HEAVYSTRING+ 'Mass Error PPM'].median()
        ldsFWHM = self.passDatasetView[LIGHTSTRING + 'Fwhm'].median()
        hdsFWHM = self.passDatasetView[HEAVYSTRING + 'Fwhm'].median()
        
        rowsLabels = ['','PPM/', 'IDP/LDP', 'RT/', 'FWHM/', 'Ht', 'MaxHt', 'Amp', 'Noise', 'totAmp', 'ms21:hlDP']        
        lTranString = [lPPM, tlPPM, lIDP, lLDP, lRT, tlRT, lFWHM, tlFWHM, lHeight, tlHeight, lArea, lBack, tlAreaAdj, lMS2MS1, hlDP]
        hTranString = [hPPM, thPPM, hIDP, hLDP, hRT, thRT, hFWHM, thFWHM, hHeight, thHeight, hArea, hBack, thAreaAdj, hMS2MS1, hlDP]
        dsTranString = [ldsPPM, hdsPPM, ldsFWHM, hdsFWHM]

        fString='{0}/{1}\n{2}/{3}\n{4}/{5}\n{6}/{7}\n{8:.1e}\n{9:.1e}\n{10:.1e}\n{11:.1e}\n{12:.1e}\n{13}:{14}'
        fdsString = '{0}/{1}\n\n\n{2}/{3}\n'
        
        self.tranLabelString = '\n'.join(rowsLabels)
        self.tranLightColumn = 'light/T\n'+fString.format(*lTranString)
        self.tranHeavyColumn = 'heavy/T\n'+fString.format(*hTranString)
        self.tranDSColumn = 'dsL/dsH\n'+fdsString.format(*dsTranString)

        return [lTranString, hTranString, dsTranString]

    def calc_prot_string(self):
        transitionSet = self.proteinView[self.proteinView['TID'] == self.currentRow['TID']]
        lPPM = transitionSet[LIGHTSTRING + 'Mass Error PPM'].median()
        hPPM = transitionSet[HEAVYSTRING + 'Mass Error PPM'].median()
        lIDP = round(transitionSet[LIGHTSTRING+'Isotope Dot Product'].median(),2)
        lLDP = round(transitionSet[LIGHTSTRING+'Library Dot Product'].median(),2)
        hIDP = round(transitionSet[HEAVYSTRING+'Isotope Dot Product'].median(),2)
        hLDP = round(transitionSet[HEAVYSTRING+'Library Dot Product'].median(),2)
        lRT = round(transitionSet[LIGHTSTRING+'Retention Time'].median(),1)
        hRT = round(transitionSet[HEAVYSTRING+'Retention Time'].median(),1)
        lRTt = round(transitionSet[LIGHTSTRING+'Best Retention Time'].median(),1)
        hRTt = round(transitionSet[HEAVYSTRING+'Best Retention Time'].median(),1)
        lFWHM = round(transitionSet[LIGHTSTRING+'Fwhm'].median(),2)
        hFWHM = round(transitionSet[HEAVYSTRING+'Fwhm'].median(),2)
        lArea = transitionSet[LIGHTSTRING+'Area'].median()
        hArea = transitionSet[HEAVYSTRING+'Area'].median()
        lBack = transitionSet[LIGHTSTRING+'Background'].median()
        hBack = transitionSet[HEAVYSTRING+'Background'].median()
        lMS2MS1 = round((transitionSet[LIGHTSTRING+'Total Area Fragment']/transitionSet[LIGHTSTRING+'Total Area MS1']).median(),2)
        hMS2MS1 = round((transitionSet[HEAVYSTRING+'Total Area Fragment']/transitionSet[HEAVYSTRING+'Total Area MS1']).median(),2)
        hlDP = round(transitionSet[HEAVYSTRING+'Ratio Dot Product'].median(),2)
        
        rowsLabels = ['','PPM', 'IDP/LDP', 'RT', 'FWHM', 'Amp', 'Noise', 'ms21:hlDP']
        lPepString = [lPPM, lIDP, lLDP, lRT, lRTt, lFWHM, lArea, lBack, lMS2MS1, hlDP]
        hPepString = [hPPM, hIDP, hLDP, hRT, hRTt, hFWHM, hArea, hBack, hMS2MS1, hlDP]
        self.protLabelString = '\n'.join(rowsLabels)
        fString='{0}\n{1}/{2}\n{3}/{4}\n{5}\n{6:.1e}\n{7:.1e}\n{8}:{9}'
        self.pepLightColumn = 'light/T\n'+fString.format(*lPepString)
        self.pepHeavyColumn = 'heavy/T\n'+fString.format(*hPepString)

    def calc_hist(self):
        name = self.currentHist
        dataClearDS=self.datasetView[self.datasetView['allClear'] == True][name].replace([numpy.inf, -numpy.inf], numpy.nan).dropna()
        dataClearProtein=self.proteinView[self.proteinView['allClear'] == True][name].replace([numpy.inf, -numpy.inf], numpy.nan).dropna()
        
        bin_num = min(30, len(list(set(dataClearDS))))
        if name == HEAVYSTRING+'AdjArea':
            bin_num = numpy.logspace(numpy.log10(dataClearDS.min()), numpy.log10(dataClearDS.max()), num=30)
        self.histPlotDataset.clear()
        self.histPlotDataset.hist(dataClearDS.values, bins = bin_num)
        self.histPlotDataset.text(0.05,0.75,name+'_DS', transform=self.histPlotDataset.transAxes)
        
        bin_num = min(30, len(list(set(dataClearProtein))))
        if name == HEAVYSTRING+'AdjArea':
            bin_num = numpy.logspace(numpy.log10(dataClearProtein.min()), numpy.log10(dataClearProtein.max()), num=30)
        self.histPlotProtein.clear()
        self.histPlotProtein.hist(dataClearProtein.values, bins = bin_num)
        self.histPlotProtein.text(0.05,0.75,name+'_Protein', transform=self.histPlotProtein.transAxes)
        if name == HEAVYSTRING+'AdjArea':
            self.histPlotProtein.set_xscale('log')
            self.histPlotDataset.set_xscale('log')
        vizLib.cleanAxis(self.histPlotDataset, ticksOnly=True)
        vizLib.cleanAxis(self.histPlotProtein, ticksOnly=True)

    def on_export_button(self, event):
        """
        Create and show the Save FileDialog
        """
        wildcard =  "Filtered Skyline CSV file (*_filt.csv)|*_filt.csv|"\
                    "All files (*.*)|*.*|"
        defFile = self.datafile[:-4]+'_filt.csv'
        dlg = wx.FileDialog(
            self, message="Save file as ...", 
            defaultDir=self.currentDirectory, 
            defaultFile=defFile, wildcard=wildcard, style=wx.SAVE
            )
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.recalcAndDrawAll()
            self.dataFrame['priorFilter'] = self.dataFrame['allClear']
            self.dataFrame.to_csv(path, index=False)
            #self.savePs(path[:-4]+'.filterParam')
        dlg.Destroy()
            
    def on_PPMtranLH_range_button(self, event):
        self.currentHist = "PPMtranLH"
        self.recalcAndDrawAll()
    def on_PPMtranTRANALL_range_button(self, event):
        self.currentHist = "PPMtranTRANALL"
        self.recalcAndDrawAll()
    def on_PPMdsL_range_button(self, event):
        self.currentHist = LIGHTSTRING + 'Mass Error PPM'
        self.recalcAndDrawAll()        
    def on_PPMdsH_range_button(self, event):
        self.currentHist = HEAVYSTRING + 'Mass Error PPM'
        self.recalcAndDrawAll()
    def on_RTdsLH_range_button(self, event):
        self.currentHist = "RTdsLH"
        self.recalcAndDrawAll()
    def on_RTdsTRANALL_range_button(self, event):
        self.currentHist = "RTdsTRANALL"
        self.recalcAndDrawAll()
    def on_RTpepTRANALL_range_button(self, event):
        self.currentHist = "RTpepTRANALL"
        self.recalcAndDrawAll()    
    def on_minIntensity_range_button(self, event):
        self.currentHist = HEAVYSTRING+'AdjArea'
        self.recalcAndDrawAll()
    def on_ratio_range_button(self, event):
        self.currentHist = "ratio"
        self.recalcAndDrawAll()
    def on_idp_range_button(self, event):
        self.currentHist = HEAVYSTRING + 'Isotope Dot Product'
        self.recalcAndDrawAll()
    def on_ldp_range_button(self, event):
        self.currentHist = HEAVYSTRING + 'Library Dot Product'
        self.recalcAndDrawAll()

    def on_checkbox(self, event):
        self.recalcAndDrawAll()
    def on_exit(self, event):
        self.Destroy()
###############################################
        #END OF CLASS DEFINTION
###############################################

        
def setCurrentFrac(calcNum, calcDen):
    num = [i.strip() for i in calcNum]
    den = [i.strip() for i in calcDen]
    return string.join(num, "+") + "/" + string.join(den, "+")
    
def openFile(fullpath):
    r = csv.reader(open(fullpath))
    header = r.next()
    
    if not (FILENAMEHEADER in header):
        print "incorrect csv file"
    else:
        dataFrame = mrmTools.readMRMCSV(fullpath, fileNameHeader = FILENAMEHEADER)
    return dataFrame

def startApp(dataFrame,fn, fsize=None, size=None):
    app = wx.App()
    screenSize = wx.DisplaySize()
    if size is None:
        size = 1.0*float(screenSize[0])/1600.0
        fsize = 1.0*float(screenSize[0])/1600.0
        print "detected screen resolution as " + str(screenSize[0]) + "; setting size=fsize="+str(round(size,2))+", these can be changed at bottom of masse.py if desired."
    vizLib.setRcs(scale=size*12, tickScale=0.75*size,  axisLineWidth=1*size)
    axisFont = {'titlesize' : size*12*1.2,
                'labelsize' : size*12*1.25}
    matplotlib.rc('axes', **axisFont)

    app.frame = MasseFrame(dataFrame, fn, fsize=fsize, size=size)
    app.frame.Show()
    app.MainLoop()    

def fileOpenStart(pathToFile=None):
        """
        Create and show the Open FileDialog
        """
        wildcard =  "All files (*.*)|*.*|"\
                    "Skyline output file (*.csv)|*.csv|"
        if pathToFile is None:            
            app = wx.App(False)  # Create a new app, don't redirect stdout/stderr to a window.
            frame = wx.Frame(None, wx.ID_ANY, "") # A Frame is a top-level window.
            dlg = wx.FileDialog(frame,
                                message="Choose a file",
                                defaultFile="",
                                wildcard=wildcard,
                                style=wx.OPEN | wx.CHANGE_DIR
                                )
            if dlg.ShowModal() == wx.ID_OK:
                pathToFile = dlg.GetPaths()[0]
        
        hold = pathToFile
        dp = '/'.join(hold.split('/')[:-1])+'/'
        fn = hold.split('/')[-1]
        print "opening file : " + pathToFile + "..."
        df = openFile(pathToFile)
        return [df, fn]
        frame.Destroy()
        app.Destroy()
        dlg.Destroy()

if __name__ == '__main__':
    if len(sys.argv) > 1:
        pathToFile=sys.argv[1]
    else:
        pathToFile=None
    #fsize=0.5
    #size=1.0
    fsize=None
    size=None
    #pathToFile = '/home/jhdavis/Dropbox/fromLinux/test_full_new.csv'
    
    [dfr, fn] = fileOpenStart(pathToFile)
    startApp(dfr, fn, fsize=fsize, size=size)
