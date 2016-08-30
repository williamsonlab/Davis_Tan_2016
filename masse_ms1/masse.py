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

class MasseFrame(wx.Frame):
    """ The main frame of the application
    """
    title = 'Masse filtering'
    
    def __init__(self, df, dp, fn, pulse=False, varLab=False, fsize=1.0, size=1.0):
        wx.Frame.__init__(self, None, -1, self.title)                
        
        self.create_main_panel(df, dp, fn, pulse=pulse, varLab=varLab, fsize=fsize, size=size)
        
        self.lowCheckNum.SetValue(True)
        self.lowCheckDen.SetValue(True)
        self.highCheckDen.SetValue(True)
        self.ppmDiffRangeBypass.SetValue('-100 100')
        self.N14RangeBypass.SetValue('-100 100')
        self.N15RangeBypass.SetValue('-100 100')
        self.missedRangeBypass.SetValue('0 5')        
        self.rtDiffRangeBypass.SetValue('-0.5 0.5')
        self.residRangeBypass.SetValue('0 5000')
        self.minIntensityBypass.SetValue('100')
        self.proteinZoomRangeBypass.SetValue('0 60')
        self.ratioLimBypass.SetValue('0.1 10')
        self.gwBypass.SetValue('0 1')

        self.ppmDiffOn.SetValue(False)
        self.N14On.SetValue(False)
        self.N15On.SetValue(False)
        self.missedOn.SetValue(False)
        self.rtOn.SetValue(False)
        self.residOn.SetValue(False)
        self.minIntensityOn.SetValue(False)
        self.ratioLimOn.SetValue(False)
        self.handSaveOn.SetValue(False)
        self.handDeleteOn.SetValue(False)
        self.gwOn.SetValue(False)
        
        if self.varLab:
            self.FRC_NXRangeBypass.SetValue('0 1')
            self.FRC_NXOn.SetValue(False)
        self.recalcAndDrawAll(setZero=True)

    def create_main_panel(self, df, dp, fn, pulse=False, varLab=False, fsize=1.0, size=1.0):
        """ Creates the main panel with all the controls on it:
             * mpl canvas 
             * mpl navigation toolbar
             * Control panel for interaction
        """
        self.fsize = fsize
        self.size = size
        self.SetFont(wx.Font(10*size, wx.SWISS, wx.NORMAL, wx.NORMAL, False,'MS Shell Dlg 2'))
        self.calcNum = ["AMP_U"]
        self.calcDen = ["AMP_U", "AMP_S"]
        self.currentHist = "ppmDiff"
        self.savedPoints=None
        self.filteredPoints=None
        #self.positionLabels = qMSDefs.positionLabels70S
        self.currentDirectory = os.getcwd()
        self.dataFrame = df
        self.datapath = dp
        self.datafile = fn
        self.pulse=pulse
        self.varLab=varLab
        self.figdim = 7.5*fsize
        self.positionLabels = qMS.sort_nicely(sorted(set(self.dataFrame['protein'].values)))
        
        self.panel = wx.Panel(self)
        
        '''Create the figure panels'''
        self.figLeft = Figure((self.figdim, self.figdim))
        self.canvasLeft = FigCanvas(self.panel, wx.ID_ANY, self.figLeft)
        self.figRight = Figure((self.figdim, self.figdim))
        self.canvasRight = FigCanvas(self.panel, wx.ID_ANY, self.figRight)
        
        gsLeft = gridspec.GridSpec(5,1)
        gsRight = gridspec.GridSpec(5,1)
        self.PLPlot = self.figLeft.add_subplot(gsLeft[:4, :])
        self.PNGPlot = self.figRight.add_subplot(gsRight[:3, :])
        self.CONPlot = self.figRight.add_subplot(gsRight[3:4, :])
        
        self.histPlotAll = self.figLeft.add_subplot(gsLeft[4, :])
        self.histPlotSelected = self.figRight.add_subplot(gsRight[4, :])
        self.figLeft.tight_layout()
        self.figRight.tight_layout()
        '''Create the list boxes'''
        self.savedList = wx.ListBox(self.panel, id=26, choices=[], style=wx.LB_SINGLE, name='Saved fits')
        self.filteredList = wx.ListBox(self.panel, id=26, choices=[], style=wx.LB_SINGLE, name='Filtered fits')
        
        '''Create the buttons'''
        self.toolbarLeft = NavigationToolbar(self.canvasLeft)
        self.toolbarRight = NavigationToolbar(self.canvasRight)
        
        self.cb_grid = wx.CheckBox(self.panel, wx.ID_ANY, label="Grid", style=wx.ALIGN_RIGHT)        
        self.hideCheck = wx.CheckBox(self.panel, wx.ID_ANY, label="Hide", style=wx.ALIGN_RIGHT)        
        self.zoomCheck = wx.CheckBox(self.panel, wx.ID_ANY, label="Zoom", style=wx.ALIGN_RIGHT)        
        self.exportButton = wx.Button(self.panel, wx.ID_ANY, "Export")
        self.openButton = wx.Button(self.panel, wx.ID_ANY, "Open")
        
        self.rother = wx.RadioButton(self.panel, label="other", style=wx.RB_GROUP)
        self.r70S = wx.RadioButton(self.panel, label="70S")
        self.r50S = wx.RadioButton(self.panel, label="50S")
        self.r30S = wx.RadioButton(self.panel, label="30S")
        
        self.savePButton = wx.Button(self.panel, wx.ID_ANY, "Save")
        self.loadPButton = wx.Button(self.panel, wx.ID_ANY, "Load")
        self.N1Button = wx.Button(self.panel, wx.ID_ANY, "N1", size=(29*size,-1))
        self.N2Button = wx.Button(self.panel, wx.ID_ANY, "1.5", size=(33*size,-1))
        self.N3Button = wx.Button(self.panel, wx.ID_ANY, "2", size=(24*size,-1))
        
        self.calcButton = wx.Button(self.panel, wx.ID_ANY, "Calculate", size=(57*size,40))
        
        self.lowCheckNum = wx.CheckBox(self.panel, wx.ID_ANY, label="low")
        if self.pulse:
            self.midCheckNum = wx.CheckBox(self.panel, wx.ID_ANY, label="mid")
        self.highCheckNum = wx.CheckBox(self.panel, wx.ID_ANY, label="high")
        
        self.lowCheckDen = wx.CheckBox(self.panel, wx.ID_ANY, label="low")
        if self.pulse:
            self.midCheckDen = wx.CheckBox(self.panel, wx.ID_ANY, label="mid")
        self.highCheckDen = wx.CheckBox(self.panel, wx.ID_ANY, label="high")
        
        self.proteinZoomOn = wx.CheckBox(self.panel, wx.ID_ANY, label="pZoom", style=wx.ALIGN_RIGHT)
        self.proteinZoomRangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)

        self.ppmDiff_range_button = wx.Button(self.panel, wx.ID_ANY, "PPM diff",size=(75*size,-1))
        self.N14_range_button = wx.Button(self.panel, wx.ID_ANY, "14N PPM", size=(75*size,-1))
        self.N15_range_button = wx.Button(self.panel, wx.ID_ANY, "15N PPM", size=(75*size,-1))
        self.missed_range_button = wx.Button(self.panel, wx.ID_ANY, "Missed cl.", size=(85*size,-1))
        self.rtDiff_range_button = wx.Button(self.panel, wx.ID_ANY, "RT diff", size=(65*size,-1))
        self.minIntensity_range_button = wx.Button(self.panel, wx.ID_ANY, "min inten.", size=(80*size,-1))
        self.resid_range_button = wx.Button(self.panel, wx.ID_ANY, "resid.", size=(70*size,-1))
        self.ratio_range_button = wx.Button(self.panel, wx.ID_ANY, "amp[U/L]", size=(80*size,-1))
        self.gw_range_button = wx.Button(self.panel, wx.ID_ANY, "GW", size=(50*size,-1))
        if self.varLab:
            self.FRC_NX_range_button = wx.Button(self.panel, wx.ID_ANY, "FRC_NX", size=(70*size,-1))
        
        self.ppmDiffRangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.N14RangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.N15RangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.missedRangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.rtDiffRangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.residRangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.minIntensityBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.ratioLimBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.gwBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        if self.varLab:
            self.FRC_NXRangeBypass = wx.TextCtrl(self.panel, size=(75,-1),style=wx.TE_PROCESS_ENTER)

        self.ppmDiffOn = wx.CheckBox(self.panel, wx.ID_ANY, label="ppmDiff")
        self.N14On = wx.CheckBox(self.panel, wx.ID_ANY, label="N14")
        self.N15On = wx.CheckBox(self.panel, wx.ID_ANY, label="N15")
        self.missedOn = wx.CheckBox(self.panel, wx.ID_ANY, label="missed")
        self.rtOn = wx.CheckBox(self.panel, wx.ID_ANY, label="RTDiff")
        self.residOn = wx.CheckBox(self.panel, wx.ID_ANY, label="residuals")
        self.minIntensityOn = wx.CheckBox(self.panel, wx.ID_ANY, label="min intensity")
        self.ratioLimOn = wx.CheckBox(self.panel, wx.ID_ANY, label="amp[U/L] ratio")
        if self.varLab:
            self.FRC_NXOn = wx.CheckBox(self.panel, wx.ID_ANY, label="FRC_NX")
        self.handSaveOn = wx.CheckBox(self.panel, wx.ID_ANY, label="curated_save")
        self.handDeleteOn = wx.CheckBox(self.panel, wx.ID_ANY, label="curated_remove")
        self.gwOn = wx.CheckBox(self.panel, wx.ID_ANY, label="GW")
        
        '''lay out the buttons'''
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.figbox = wx.BoxSizer(wx.HORIZONTAL)
        self.figbox.Add(self.canvasLeft, 0)
        self.listBox = wx.BoxSizer(wx.VERTICAL)
        self.savedListBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'saved fits'), wx.HORIZONTAL)
        self.savedListBox.Add(self.savedList, 1, flag=wx.GROW)
        
        self.filteredListBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'filtered fits'), wx.HORIZONTAL)
        self.filteredListBox.Add(self.filteredList, 1, flag=wx.GROW)
        
        self.listBox.Add(self.savedListBox, 1, flag=wx.GROW)        
        self.listBox.Add(self.filteredListBox, 1, flag=wx.GROW)
        self.figbox.Add(self.listBox, 1, flag=wx.GROW)
        self.figbox.Add(self.canvasRight, 0)
        self.vbox.Add(self.figbox, 0, flag = wx.GROW)

        #new horizontal box    
        self.toolNumBox = wx.BoxSizer(wx.HORIZONTAL)
        
        #add toolbar, draw, grid and exportto hbox
        self.imageToolsBoxL = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'image tools'), wx.HORIZONTAL)

        self.zoomBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'zoom tools'), wx.VERTICAL)
        self.zoomBox.Add(self.proteinZoomOn, 0, flag=wx.ALIGN_LEFT)
        self.zoomBox.Add(self.proteinZoomRangeBypass, 0, flag=wx.ALIGN_LEFT)       
        
        self.toolbarSubunits = wx.BoxSizer(wx.VERTICAL)
        self.toolbarSubunits.Add(self.toolbarLeft, 0, flag=wx.ALIGN_LEFT)
        self.subunitsBox = wx.BoxSizer(wx.HORIZONTAL)
        self.subunitsBox.Add(self.r70S, 0, flag=wx.ALIGN_LEFT)
        self.subunitsBox.Add(self.r50S, 0, flag=wx.ALIGN_LEFT)
        self.subunitsBox.Add(self.r30S, 0, flag=wx.ALIGN_LEFT)
        self.subunitsBox.Add(self.rother, 0, flag=wx.ALIGN_LEFT)
        self.toolbarSubunits.Add(self.subunitsBox, 0, flag=wx.ALIGN_CENTER | wx.ALIGN_CENTER)
        self.imageToolsBoxL.Add(self.toolbarSubunits, 0, flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER)
        
        self.imageCheckToolsBox = wx.BoxSizer(wx.VERTICAL)        
        self.imageCheckToolsBox.Add(self.cb_grid, 0, flag=wx.ALIGN_RIGHT)
        self.imageCheckToolsBox.Add(self.hideCheck, 0, flag=wx.ALIGN_RIGHT)
        self.imageCheckToolsBox.Add(self.zoomCheck, 0, flag=wx.ALIGN_RIGHT)
        self.imageToolsBoxL.Add(self.imageCheckToolsBox, 0, flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER)
        self.toolNumBox.Add(self.imageToolsBoxL, 0, flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER)

        #add open and export buttons
        self.fileToolsBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'file tools'), wx.VERTICAL)
        self.fileToolsBox.Add(self.openButton, 0, flag=wx.ALIGN_LEFT)
        self.fileToolsBox.Add(self.exportButton, 0, flag=wx.ALIGN_LEFT)        
        self.toolNumBox.Add(self.fileToolsBox, 0, flag=wx.ALIGN_LEFT | wx.GROW)
        
        #add loadP/SaveP buttons
        self.paramBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'filterParams'), wx.VERTICAL)
        self.paramBox.Add(self.loadPButton, 0, flag=wx.ALIGN_LEFT)
        self.paramBox.Add(self.savePButton, 0, flag=wx.ALIGN_LEFT)
        self.toolNumBox.Add(self.paramBox, 0, flag=wx.ALIGN_RIGHT | wx.GROW)
        
        
        #add calculate button
        self.calcNBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'calcDist'), wx.VERTICAL)
        self.calcNBox.Add(self.calcButton, flag=wx.ALIGN_TOP | wx.GROW)
        self.NBox = wx.BoxSizer(wx.HORIZONTAL)
        self.NBox.Add(self.N1Button, 0, flag=wx.ALIGN_LEFT)
        self.NBox.Add(self.N2Button, 0, flag=wx.ALIGN_LEFT)
        self.NBox.Add(self.N3Button, 0, flag=wx.ALIGN_LEFT)
        self.calcNBox.Add(self.NBox, 0, flag=wx.ALIGN_TOP)
        self.toolNumBox.Add(self.calcNBox, 0, flag=wx.ALIGN_RIGHT | wx.GROW)
        self.toolNumBox.Add(self.zoomBox, 0, flag=wx.ALIGN_LEFT)
        
        #add numerator box to hbox
        self.numBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'num'), wx.VERTICAL)
        self.numBox.Add(self.lowCheckNum, 0, flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER)
        if self.pulse:
            self.numBox.Add(self.midCheckNum, 0, flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER)
        self.numBox.Add(self.highCheckNum, 0,flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER)
        self.toolNumBox.Add(self.numBox, 0, flag=wx.ALIGN_RIGHT | wx.GROW)
        #add denominator box to hbox
        self.denBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'den'), wx.VERTICAL)
        self.denBox.Add(self.lowCheckDen, 0, flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER)
        if self.pulse:
            self.denBox.Add(self.midCheckDen, 0, flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER)
        self.denBox.Add(self.highCheckDen, 0,flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER)
        self.toolNumBox.Add(self.denBox, 0, flag=wx.ALIGN_RIGHT | wx.GROW)
        
        self.imageToolsBoxR = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'image tools'), wx.HORIZONTAL)
        self.imageToolsBoxR.Add(self.toolbarRight, 0, flag=wx.ALIGN_RIGHT)
        self.expandBox = wx.BoxSizer(wx.HORIZONTAL)
        self.toolNumBox.Add(self.expandBox, 1, flag=wx.GROW)
        self.toolNumBox.Add(self.imageToolsBoxR, 0, flag=wx.ALIGN_RIGHT | wx.GROW)
        
        #add hbox to vbox
        self.vbox.Add(self.toolNumBox, 0, flag = wx.GROW)
        
        # Sliders for setting the various cutoffs
        self.filterBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'filter cuttoff settings'), wx.HORIZONTAL)
        flags = wx.ALIGN_LEFT | wx.ALIGN_CENTER
        self.filterBox.Add(self.ppmDiff_range_button, 0, flag=flags)
        self.filterBox.Add(self.ppmDiffRangeBypass, 0, flag=flags)
        self.filterBox.Add(self.N14_range_button, 0, flag=flags)
        self.filterBox.Add(self.N14RangeBypass, 0, flag=flags)
        self.filterBox.Add(self.N15_range_button, 0, flag=flags)
        self.filterBox.Add(self.N15RangeBypass, 0, flag = flags)
        self.filterBox.Add(self.missed_range_button, 0, flag=flags)
        self.filterBox.Add(self.missedRangeBypass, 0, flag = flags)
        self.filterBox.Add(self.rtDiff_range_button, 0, flag=flags)
        self.filterBox.Add(self.rtDiffRangeBypass, 0, flag = flags)
        if self.varLab:
            self.filterBox.Add(self.FRC_NX_range_button, 0, flag=flags)
            self.filterBox.Add(self.FRC_NXRangeBypass, 0, flag = flags)
        self.filterBox.Add(self.resid_range_button, 0, flag=flags)
        self.filterBox.Add(self.residRangeBypass, 0, flag = flags)
        self.filterBox.Add(self.minIntensity_range_button, 0, flag=flags)
        self.filterBox.Add(self.minIntensityBypass, 0, flag = flags)
        self.filterBox.Add(self.ratio_range_button, 0, flag=flags)
        self.filterBox.Add(self.ratioLimBypass, 0, flag = flags)
        self.filterBox.Add(self.gw_range_button, 0, flag=flags)
        self.filterBox.Add(self.gwBypass, 0, flag = flags)
        self.vbox.Add(self.filterBox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        
        # Filters on/off checks
        self.controlFilters = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'filter on/off'), wx.HORIZONTAL)
        self.controlFilters.Add(self.ppmDiffOn, 0, flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER)
        self.controlFilters.Add(self.N14On, 0, flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER)
        self.controlFilters.Add(self.N15On, 0, flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER)
        self.controlFilters.Add(self.missedOn, 0, flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER)
        self.controlFilters.Add(self.rtOn, 0, flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER)
        if self.varLab:
            self.controlFilters.Add(self.FRC_NXOn, 0, flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER)
        self.controlFilters.Add(self.residOn, 0, flag = wx.ALIGN_RIGHT | wx.ALIGN_CENTER)
        self.controlFilters.Add(self.minIntensityOn, 0, flag = wx.ALIGN_RIGHT | wx.ALIGN_CENTER)
        self.controlFilters.Add(self.handSaveOn, 0, flag = wx.ALIGN_RIGHT | wx.ALIGN_CENTER)
        self.controlFilters.Add(self.handDeleteOn, 0, flag = wx.ALIGN_RIGHT | wx.ALIGN_CENTER)
        self.controlFilters.Add(self.ratioLimOn, 0, flag = wx.ALIGN_RIGHT | wx.ALIGN_CENTER)
        self.controlFilters.Add(self.gwOn, 0, flag = wx.ALIGN_RIGHT | wx.ALIGN_CENTER)
        self.vbox.Add(self.controlFilters, 0, flag = wx.ALIGN_LEFT | wx.TOP)

        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)
        
        '''bind events for buttons'''
        self.canvasLeft.mpl_connect('pick_event', self.pickScatterPoint)
        
        self.r70S.Bind(wx.EVT_RADIOBUTTON, self.onr70S_select)
        self.r50S.Bind(wx.EVT_RADIOBUTTON, self.onr50S_select)
        self.r30S.Bind(wx.EVT_RADIOBUTTON, self.onr30S_select)
        self.rother.Bind(wx.EVT_RADIOBUTTON, self.onrother_select)

        self.exportButton.Bind(wx.EVT_BUTTON, self.on_export_button)
        self.openButton.Bind(wx.EVT_BUTTON, self.on_open_button)
        self.calcButton.Bind(wx.EVT_BUTTON, self.on_calc_button)
        self.cb_grid.Bind(wx.EVT_CHECKBOX, self.on_cb_grid)
        self.hideCheck.Bind(wx.EVT_CHECKBOX, self.on_hideCheck)
        self.zoomCheck.Bind(wx.EVT_CHECKBOX, self.on_zoomCheck)
        self.proteinZoomOn.Bind(wx.EVT_CHECKBOX, self.on_pzoom_check)
        
        self.loadPButton.Bind(wx.EVT_BUTTON, self.on_loadP_button)
        self.savePButton.Bind(wx.EVT_BUTTON, self.on_saveP_button)
        self.N1Button.Bind(wx.EVT_BUTTON, self.on_N1_button)
        self.N2Button.Bind(wx.EVT_BUTTON, self.on_N2_button)
        self.N3Button.Bind(wx.EVT_BUTTON, self.on_N3_button)

        self.ppmDiff_range_button.Bind(wx.EVT_BUTTON, self.on_ppmDiff_range_button)
        self.N14_range_button.Bind(wx.EVT_BUTTON, self.on_N14_range_button)
        self.N15_range_button.Bind(wx.EVT_BUTTON, self.on_N15_range_button)
        self.missed_range_button.Bind(wx.EVT_BUTTON, self.on_missed_range_button)
        self.rtDiff_range_button.Bind(wx.EVT_BUTTON, self.on_rtDiff_range_button)
        self.resid_range_button.Bind(wx.EVT_BUTTON, self.on_resid_range_button)
        self.minIntensity_range_button.Bind(wx.EVT_BUTTON, self.on_minIntensity_range_button)
        self.ratio_range_button.Bind(wx.EVT_BUTTON, self.on_ratio_range_button)
        self.gw_range_button.Bind(wx.EVT_BUTTON, self.on_gw_range_button)
        if self.varLab:
            self.FRC_NX_range_button.Bind(wx.EVT_BUTTON, self.on_FRC_NX_range_button)
        
        self.ppmDiffOn.Bind(wx.EVT_CHECKBOX, self.on_ppmDiffOn)
        self.N14On.Bind(wx.EVT_CHECKBOX, self.on_N14On)
        self.N15On.Bind(wx.EVT_CHECKBOX, self.on_N15On)
        self.missedOn.Bind(wx.EVT_CHECKBOX, self.on_missedOn)
        self.rtOn.Bind(wx.EVT_CHECKBOX, self.on_rtOn)
        self.residOn.Bind(wx.EVT_CHECKBOX, self.on_residOn)
        self.minIntensityOn.Bind(wx.EVT_CHECKBOX, self.on_minIntensityOn)
        self.handSaveOn.Bind(wx.EVT_CHECKBOX, self.on_handSaveOn)
        self.handDeleteOn.Bind(wx.EVT_CHECKBOX, self.on_handDeleteOn)
        self.ratioLimOn.Bind(wx.EVT_CHECKBOX, self.on_ratioLimOn)
        self.gwOn.Bind(wx.EVT_CHECKBOX, self.on_gwOn)
        if self.varLab:
            self.FRC_NXOn.Bind(wx.EVT_CHECKBOX, self.on_FRC_NXOn)
        
        self.savedList.Bind(wx.EVT_LISTBOX, self.on_savedBoxClick)
        self.filteredList.Bind(wx.EVT_LISTBOX, self.on_filteredBoxClick)
        
        self.panel.Bind(wx.EVT_KEY_UP, self.on_key_press)

    def on_N1_button(self, event):
        self.setNormalValues(1.0)
    def on_N2_button(self, event):
        self.setNormalValues(1.5)
    def on_N3_button(self, event):
        self.setNormalValues(2.0)
    def setNormalValues(self, stds):
        muPPM = self.dataFrame['ppmDiff'].median()
        varPPM = self.dataFrame['ppmDiff'].std()
        
        mu14PPM = self.dataFrame['ppm_n14'].median()
        var14PPM = self.dataFrame['ppm_n14'].std()
        mu14PPM = self.dataFrame[(self.dataFrame['ppm_n14'] > mu14PPM-3*var14PPM) & (self.dataFrame['ppm_n14'] < mu14PPM+3*var14PPM)]['ppm_n14'].median()
        
        mu15PPM = self.dataFrame['ppm_n15'].median()
        var15PPM = self.dataFrame['ppm_n15'].std()
        mu15PPM = self.dataFrame[(self.dataFrame['ppm_n15'] > mu15PPM-3*var15PPM) & (self.dataFrame['ppm_n15'] < mu15PPM+3*var15PPM)]['ppm_n15'].median()
        
        if self.varLab:
            muFRCNX = self.dataFrame['FRC_NX'].median()
            varFRCNX = self.dataFrame['FRC_NX'].std()
        
        muRatio = self.dataFrame[(self.dataFrame['ratio'] > 0) & (self.dataFrame['ratio'] < 20)]['ratio'].median()
        varRatio = self.dataFrame[(self.dataFrame['ratio'] > 0) & (self.dataFrame['ratio'] < 20)]['ratio'].std()
        
        self.ppmDiffRangeBypass.SetValue(str(round(muPPM-stds*varPPM,2)) + ' ' + str(round(muPPM+stds*varPPM,2)))
        self.N14RangeBypass.SetValue(str(round(mu14PPM-stds*var14PPM,2)) + ' ' + str(round(mu14PPM+stds*var14PPM,2)))
        self.N15RangeBypass.SetValue(str(round(mu15PPM-stds*var15PPM,2)) + ' ' + str(round(mu15PPM+stds*var15PPM,2)))
        #self.ratioLimBypass.SetValue(str(round(muRatio-stds*varRatio,2)) + ' ' + str(round(muRatio+stds*varRatio,2)))
        
        if self.varLab:
            self.FRC_NXRangeBypass.SetValue(str(round(muFRCNX-stds*varFRCNX,2)) + ' ' + str(round(muFRCNX+stds*varFRCNX,2)))

    def on_loadP_button(self, event):
        f = open(self.datapath+'_last.filterParam', 'r')
        params = csv.reader(f)
        pdict = {}
        for l in params:
            pdict[l[0]] = l[1]
        f.close()
        self.lowCheckNum.SetValue(qMS.boolParse(pdict['lowNum']))
        if self.pulse:
            self.midCheckNum.SetValue(qMS.boolParse(pdict['midNum']))
        self.highCheckNum.SetValue(qMS.boolParse(pdict['highNum']))
        self.lowCheckDen.SetValue(qMS.boolParse(pdict['lowDen']))
        if self.pulse:
            self.midCheckDen.SetValue(qMS.boolParse(pdict['midDen']))
        self.highCheckDen.SetValue(qMS.boolParse(pdict['highDen']))
        
        self.cb_grid.SetValue(qMS.boolParse(pdict['gridChecked']))
        self.zoomCheck.SetValue(qMS.boolParse(pdict['zoomChecked']))
        self.hideCheck.SetValue(qMS.boolParse(pdict['hideChecked']))
        
        self.ppmDiffRangeBypass.SetValue(pdict['ppmDiff_low'] + ' ' + pdict['ppmDiff_high'])
        self.N14RangeBypass.SetValue(pdict['ppm_n14_low'] + ' ' + pdict['ppm_n14_high'])
        self.N15RangeBypass.SetValue(pdict['ppm_n15_low'] + ' ' + pdict['ppm_n15_high'])
        self.missedRangeBypass.SetValue(pdict['missed_low'] + ' ' + pdict['missed_high'])
        self.rtDiffRangeBypass.SetValue(pdict['rtDiff_low'] + ' ' + pdict['rtDiff_high'])
        if self.varLab:
            self.FRC_NXRangeBypass.SetValue(pdict['FRC_NX_low'] + ' ' + pdict['FRC_NX_high'])
            self.FRC_NXOn.SetValue(qMS.boolParse(pdict['FRC_NX']))
        self.residRangeBypass.SetValue(pdict['resid_low'] + ' ' + pdict['resid_high'])
        self.ratioLimBypass.SetValue(pdict['ratio_low'] + ' ' + pdict['ratio_high'])
        self.gwBypass.SetValue(pdict['GW_low'] + ' ' + pdict['GW_high'])
        self.minIntensityBypass.SetValue(pdict['minIntensity'])

        self.ppmDiffOn.SetValue(qMS.boolParse(pdict['ppmDiff']))
        self.N14On.SetValue(qMS.boolParse(pdict['n14']))
        self.N15On.SetValue(qMS.boolParse(pdict['n15']))
        self.missedOn.SetValue(qMS.boolParse(pdict['missed']))
        self.rtOn.SetValue(qMS.boolParse(pdict['rtDiff']))
        self.residOn.SetValue(qMS.boolParse(pdict['resid']))
        self.ratioLimOn.SetValue(qMS.boolParse(pdict['ratio']))
        self.gwOn.SetValue(qMS.boolParse(pdict['GW']))
        self.minIntensityOn.SetValue(qMS.boolParse(pdict['minIntensityOn']))
        self.handSaveOn.SetValue(qMS.boolParse(pdict['handSaveOn']))
        self.handDeleteOn.SetValue(qMS.boolParse(pdict['handDeleteOn']))
        f.close()
        self.recalcAndDrawAll()

    def on_saveP_button(self, event):
        self.savePs(self.datapath+'_last.filterParam')
    def savePs(self, p):
        (ppmDiff_low, ppmDiff_high) = map(str, self.ppmDiffRangeBypass.GetValue().split(' '))
        (ppm_n14_low, ppm_n14_high) = map(str, self.N14RangeBypass.GetValue().split(' '))
        (ppm_n15_low, ppm_n15_high) = map(str, self.N15RangeBypass.GetValue().split(' '))
        (missed_low, missed_high) = map(str, self.missedRangeBypass.GetValue().split(' '))
        (rtDiff_low, rtDiff_high) = map(str, self.rtDiffRangeBypass.GetValue().split(' '))
        if self.varLab:
            (FRC_NX_low, FRC_NX_high) = map(str, self.FRC_NXRangeBypass.GetValue().split(' '))
        (resid_low, resid_high) = map(str, self.residRangeBypass.GetValue().split(' '))
        (ratio_low, ratio_high) = map(str, self.ratioLimBypass.GetValue().split(' '))
        (gw_low, gw_high) = map(str, self.gwBypass.GetValue().split(' '))
        minI = float(self.minIntensityBypass.GetValue())

        outstr = 'param,' + 'value'
        
        outstr = outstr + '\nppmDiff_low,' + ppmDiff_low + '\nppmDiff_high,' + ppmDiff_high
        outstr = outstr + '\nppm_n14_low,' + ppm_n14_low + '\nppm_n14_high,' + ppm_n14_high
        outstr = outstr + '\nppm_n15_low,' + ppm_n15_low + '\nppm_n15_high,' + ppm_n15_high
        outstr = outstr + '\nmissed_low,' + missed_low + '\nmissed_high,' + missed_high
        outstr = outstr + '\nrtDiff_low,' + rtDiff_low + '\nrtDiff_high,' + rtDiff_high
        if self.varLab:
            outstr = outstr + '\nFRC_NX_low,' + FRC_NX_low + '\nFRC_NX_high,' + FRC_NX_high
            outstr = outstr + '\nFRC_NX,' + str(self.FRC_NXOn.IsChecked())
        outstr = outstr + '\nresid_low,' + resid_low + '\nresid_high,' + resid_high
        outstr = outstr + '\nratio_low,' + ratio_low + '\nratio_high,' + ratio_high
        outstr = outstr + '\nGW_low,' + gw_low + '\nGW_high,' + gw_high
        outstr = outstr + '\nminIntensity,' + str(minI)
        
        outstr = outstr + '\ngridChecked,' + str(self.cb_grid.IsChecked())
        outstr = outstr + '\nzoomChecked,' + str(self.zoomCheck.IsChecked())
        outstr = outstr + '\nhideChecked,' + str(self.hideCheck.IsChecked())

        outstr = outstr + '\nlowNum,' + str(self.lowCheckNum.IsChecked())
        if self.pulse:
            outstr = outstr + '\nmidNum,' + str(self.midCheckNum.IsChecked())
        outstr = outstr + '\nhighNum,' + str(self.highCheckNum.IsChecked())
        outstr = outstr + '\nlowDen,' + str(self.lowCheckDen.IsChecked())
        if self.pulse:
            outstr = outstr + '\nmidDen,' + str(self.midCheckDen.IsChecked())
        outstr = outstr + '\nhighDen,' + str(self.highCheckDen.IsChecked())
        
        outstr = outstr + '\nppmDiff,' + str(self.ppmDiffOn.IsChecked())
        outstr = outstr + '\nn14,' + str(self.N14On.IsChecked())
        outstr = outstr + '\nn15,' + str(self.N15On.IsChecked())
        outstr = outstr + '\nmissed,' + str(self.missedOn.IsChecked())
        outstr = outstr + '\nrtDiff,' + str(self.rtOn.IsChecked())
        outstr = outstr + '\nresid,' + str(self.residOn.IsChecked())
        outstr = outstr + '\nratio,' + str(self.ratioLimOn.IsChecked())
        outstr = outstr + '\nGW,' + str(self.gwOn.IsChecked())
        outstr = outstr + '\nminIntensityOn,' + str(self.minIntensityOn.IsChecked())
        outstr = outstr + '\nhandSaveOn,' + str(self.handSaveOn.IsChecked())
        outstr = outstr + '\nhandDeleteOn,' + str(self.handDeleteOn.IsChecked())

        f = open(p, 'w')
        f.write(outstr)
        f.close()
        
    def onr70S_select(self, event):
        self.dataFrame['currentPos'] = self.dataFrame['70Spos']
        self.positionLabels = qMSDefs.positionLabels70S
        self.recalcAndDrawAll()
        
    def onr50S_select(self, event):
        self.dataFrame['currentPos'] = self.dataFrame['50Spos']        
        self.positionLabels = qMSDefs.positionLabels50S
        self.recalcAndDrawAll()
        
    def onr30S_select(self, event):
        self.dataFrame['currentPos'] = self.dataFrame['30Spos']
        self.positionLabels = qMSDefs.positionLabels30S
        self.recalcAndDrawAll()
        
    def onrother_select(self, event):
        self.dataFrame['currentPos'] = self.dataFrame['otherpos']
        self.positionLabels = qMS.sort_nicely(sorted(set(self.dataFrame['protein'].values)))
        self.recalcAndDrawAll()

    def on_open_button(self, event):
        """
        Create and show the Open FileDialog
        """
        wildcard =  "All files (*.*)|*.*|"\
                    "Preprocessed _iso_res.csv file (*_iso_res.csv)|*_iso_res.csv|"\
                    "Massacre iso_csv file (*_iso.csv)|*_iso.csv|"
        dlg = wx.FileDialog(
            self, message="Choose a file",
            defaultDir=self.currentDirectory, 
            defaultFile="",
            wildcard=wildcard,
            style=wx.OPEN | wx.CHANGE_DIR
            )
        
        if dlg.ShowModal() == wx.ID_OK:
            fullname = dlg.GetPaths()[0].split('/')
            dpa = '/'.join(fullname[:-1]) + '/'
            self.currentDirectory = dpa
            fna = fullname[-1]
            [dfr, pul, vlab] = openFile(dpa+fna)
            startApp(dfr, dpa, fna, pul, vlab, fsize=self.fsize, size=self.size)

        dlg.Destroy()
    def on_export_button(self, event):
        """
        Create and show the Save FileDialog
        """
        wildcard =  "Filtered _iso_res_filt.csv file (*_iso_res_filt.csv)|*_iso_res_filt.csv|"\
                    "All files (*.*)|*.*|"
        defFile = self.datafile[:-4]+'_filt.csv'
        dlg = wx.FileDialog(
            self, message="Save file as ...", 
            defaultDir=self.currentDirectory, 
            defaultFile=defFile, wildcard=wildcard, style=wx.SAVE
            )
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.calcPoints()
            self.savedPoints.to_csv(path, index=False)
            print path
            self.savePs(path[:-4]+'.filterParam')
        dlg.Destroy()
        
    def on_cb_grid(self, event):
        self.calc_figure()
        self.canvasLeft.draw()
    def on_zoomCheck(self, event):
        self.calc_figure()
        self.canvasLeft.draw()
    def on_hideCheck(self, event):
        self.calc_figure()
        self.canvasLeft.draw()
    
    def on_pzoom_check(self, event):
        self.calc_figure()
        self.canvasLeft.draw()
        
    def on_ppmDiff_range_button(self, event):
        self.currentHist = "ppmDiff"
        self.recalcAndDrawAll()
    def on_N14_range_button(self, event):
        self.currentHist = "ppm_n14"
        self.recalcAndDrawAll()        
    def on_N15_range_button(self, event):
        self.currentHist = "ppm_n15"
        self.recalcAndDrawAll()
    def on_missed_range_button(self, event):
        self.currentHist = "missed"
        self.recalcAndDrawAll()
    def on_rtDiff_range_button(self, event):
        self.currentHist = "rtDiff"
        self.recalcAndDrawAll()
    def on_FRC_NX_range_button(self, event):
        self.currentHist = "FRC_NX"
        self.recalcAndDrawAll()
    def on_resid_range_button(self, event):
        self.currentHist = "resid"
        self.recalcAndDrawAll()
    def on_minIntensity_range_button(self, event):
        self.currentHist = "minIntensity"
        self.recalcAndDrawAll()
    def on_ratio_range_button(self, event):
        self.currentHist = "ratio"
        self.recalcAndDrawAll()
    def on_gw_range_button(self, event):
        self.currentHist = "GW"
        self.recalcAndDrawAll()
    def on_calc_button(self, event):
        self.recalcAndDrawAll()

    def on_ppmDiffOn(self, event):
        self.recalcAndDrawAll()
    def on_N14On(self, event):
        self.recalcAndDrawAll()
    def on_N15On(self, event):
        self.recalcAndDrawAll()
    def on_missedOn(self, event):
        self.recalcAndDrawAll()
    def on_rtOn(self, event):
        self.recalcAndDrawAll()
    def on_residOn(self, event):
        self.recalcAndDrawAll()
    def on_FRC_NXOn(self, event):
        self.recalcAndDrawAll()
    def on_minIntensityOn(self, event):
        self.recalcAndDrawAll()
    def on_ratioLimOn(self, event):
        self.recalcAndDrawAll()
    def on_handSaveOn(self, event):
        self.recalcAndDrawAll()
    def on_handDeleteOn(self, event):
        self.recalcAndDrawAll()
    def on_gwOn(self, event):
        self.recalcAndDrawAll()

    def on_exit(self, event):
        self.Destroy()

    def on_savedBoxClick(self, event):
        if not (self.savedList.GetStringSelection() is u''):
            self.currentISOFile = self.savedList.GetStringSelection()
            self.currentRow = self.dataFrame[self.dataFrame['isofile'] == self.currentISOFile]
            self.newSelection()
        
    def on_filteredBoxClick(self, event):
        if not (self.filteredList.GetStringSelection() is u''):
            self.currentISOFile = self.filteredList.GetStringSelection()
            self.currentRow = self.dataFrame[self.dataFrame['isofile'] == self.currentISOFile]
            self.newSelection()

    def on_key_press(self, event):
        event.Skip()
        c = event.GetKeyCode()
        if c is 68: #got a d keystroke
            rowIX = self.dataFrame[self.dataFrame['isofile']==self.currentISOFile].index.values[0]
            self.dataFrame.loc[rowIX, 'handDelete'] = True
            self.dataFrame.loc[rowIX, 'handSave'] = False

            #row = self.dataFrame[self.dataFrame['isofile'] == self.currentISOFile]
            #row['handDelete'] = True
            #row['handSave'] = False

            myIndex = list(self.savedListItems).index(self.currentISOFile)
            try:
                nextItem = self.savedListItems[myIndex+1]
            except IndexError:
                nextItem = self.savedListItems[myIndex-1]
            self.currentISOFile = nextItem
            self.currentRow = self.dataFrame[self.dataFrame['isofile'] == self.currentISOFile]

            self.newSelection()            
            #self.dataFrame.update(row)
            self.recalcAndDrawAll()
            
        elif c is 83: # got a s keystroke
            #row = self.dataFrame[self.dataFrame['isofile'] == self.currentISOFile]
            #row['handDelete'] = False
            #row['handSave'] = True
            rowIX = self.dataFrame[self.dataFrame['isofile']==self.currentISOFile].index.values[0]
            self.dataFrame.loc[rowIX, 'handDelete'] = False
            self.dataFrame.loc[rowIX, 'handSave'] = True      
        
            #self.dataFrame.update(row)
            self.recalcAndDrawAll()

        
    def newSelection(self):
        self.calc_fit()
        self.canvasRight.draw()
        self.selectedPoint.set_data(self.currentRow['currentPos'].values[0], self.currentRow['currentCalc'].values[0])
        self.canvasLeft.draw()
    
    def recalcAndDrawAll(self, setZero=False):
        self.calc_data()
        self.calcPoints()
        if setZero is True:
            self.currentRow = self.savedPoints[0:1]
            self.currentISOFile = self.currentRow['isofile'].values[0]
        else:
            self.currentRow = self.dataFrame[self.dataFrame['isofile'] == self.currentISOFile]
        self.updateLists()
        self.calc_figure()
        self.calc_hist()
        self.calc_fit()
        self.draw_all()
    
    def calc_data(self):
        self.UID_output_list = []  
        self.calcNum = self.getChecksNum()
        self.calcDen = self.getChecksDen()
        self.dataFrame['currentCalc'] = qMS.calcValue(self.dataFrame, self.calcNum, self.calcDen)
        
    def calcPoints(self):
        self.savedPoints = self.getPass(True)
        self.filteredPoints = self.getPass(False)

    def updateLists(self):
        self.savedListItems = self.savedPoints['isofile'].values
        sli = qMS.sort_nicely(list(self.savedListItems))
        self.savedListItems = numpy.array(sli)
        self.filteredListItems = self.filteredPoints['isofile'].values
        fli = qMS.sort_nicely(list(self.filteredListItems))
        self.filteredListItems = numpy.array(fli)
        self.savedList.Set(self.savedListItems)
        self.filteredList.Set(self.filteredListItems)
        self.savedList.SetStringSelection(self.currentISOFile)
        self.filteredList.SetStringSelection(self.currentISOFile)

    def determineMedians(self):
        xs = list(set(self.savedPoints['currentPos']))
        ys = [self.savedPoints[self.savedPoints['currentPos']==k]['currentCalc'].median() for k in xs]
        return [xs,ys]
        
    def calc_figure(self):
        self.PLPlot.clear()
        self.selectedPoint, = self.PLPlot.plot(self.currentRow['currentPos'].values[0], self.currentRow['currentCalc'].values[0], 
                                              'o', ms=20, alpha=0.5, color='yellow', visible=True)
        self.PLPlot.grid(self.cb_grid.IsChecked())
        self.PLPlot.plot(self.savedPoints['currentPos'], self.savedPoints['currentCalc'], 'ro', picker=5, label="Saved : " + str(len(self.savedPoints['currentCalc'].values)))
        meds = self.determineMedians()
        self.PLPlot.plot(meds[0], meds[1], 'g-', lw=2, label="Median : " + str(round(numpy.median(meds[1]),1)))
        if self.hideCheck.IsChecked():
            self.PLPlot.plot(self.filteredPoints['currentPos'], self.filteredPoints['currentCalc'], 'x', mec='grey', picker=5, label="Filtered : " + str(len(self.filteredPoints['currentCalc'].values)))
        else:
            self.PLPlot.plot(self.filteredPoints['currentPos'], self.filteredPoints['currentCalc'], 'bo', picker=5, label="Filtered : "  + str(len(self.filteredPoints['currentCalc'].values)))
        self.PLPlot.set_xticks(range(1,int(self.dataFrame['currentPos'].max())+1))
        if self.positionLabels[0].count("-ECOLI")>0:
            self.PLPlot.set_xticklabels([i.split('-')[0] for i in self.positionLabels], rotation=90, size=10)
        else:
            self.PLPlot.set_xticklabels(self.positionLabels, rotation=90, size=10)
        self.PLPlot.set_title(self.datafile + " : " + setCurrentFrac(self.calcNum, self.calcDen))
        self.PLPlot.set_xlim([0,self.dataFrame['currentPos'].max()+1])
        if not self.zoomCheck.IsChecked():
            self.PLPlot.set_ylim([0,max(self.savedPoints['currentCalc'].max(),1)])
        else:
            self.PLPlot.set_ylim([0,3])
        
        if self.proteinZoomOn.IsChecked():
            (x_low, x_high) = map(float, self.proteinZoomRangeBypass.GetValue().split(' '))
            self.PLPlot.set_xlim([x_low,x_high])
        else:
            self.PLPlot.set_xlim([0.5,self.dataFrame['currentPos'].max()+1.5])        
        self.PLPlot.legend()
        vizLib.cleanAxis(self.PLPlot, ticksOnly=True)
        
        """    
    def calc_fit(self):
        plotsfile = self.datapath+self.currentISOFile+'.plots'
        txtfile = self.datapath+self.currentISOFile+'.txt'
        df = pd.read_csv(plotsfile)
        df2 = pd.read_csv(txtfile, header=None, sep=' ')
        del df['residAdj']
        self.PNGPlot.clear()
        #self.PNGPlot.plot(df2[0][0:len(df['dat'])].values, df['dat'].values, 'o', markersize=3, label='data', markerfacecolor='None', markeredgecolor='red')
        self.PNGPlot.plot(df2[0][0:len(df['dat'])].values, df['dat'].values, 'o', markersize=6, markerfacecolor='None', markeredgecolor='red')
        self.PNGPlot.plot(df2[0][0:len(df['dat'])].values, df['dat'].values, 'r-', linewidth=2, label='data')
        self.PNGPlot.plot(df2[0][0:len(df['fit'])].values, df['fit'].values, 'b-', linewidth=2, label='fit')
        self.PNGPlot.plot(df2[0][0:len(df['resid'])].values, df['resid'].values, 'g-', linewidth=2, label='residual')
        self.PNGPlot.set_xlabel('m/z')
        self.PNGPlot.set_ylabel('intensity')
        self.PNGPlot.set_xlim(df2[0][0:len(df['dat'])].values.min(), df2[0][0:len(df['dat'])].values.max())
        self.PNGPlot.legend()
        row = self.dataFrame[self.dataFrame['isofile'] == self.currentISOFile]
        passing = self.testPassRow(row)
        stringColor = 'black'
        if row['handDelete'].values[0] is True:
            stringColor = 'red'
        elif row['handSave'].values[0] is True:
            stringColor = 'green'
        dataString =    "ppmDiff: " + str(round(row['ppmDiff'].values[0],1)) + " : " + str(passing['ppmDiff']) +\
                        "\nN14: " + str(round(row['ppm_n14'].values[0],1)) + " : " + str(passing['N14']) +\
                        "\nN15: " + str(round(row['ppm_n15'].values[0],1)) + " : " + str(passing['N15']) +\
                        "\nmissed: " + str(row['missed'].values[0]) + " : " + str(passing['missed']) +\
                        "\nrtDiff: " + str(round(row['rtDiff'].values[0],3)) + " : " + str(passing['rtDiff']) +\
                        "\nresid: " + str(round(row['resid'].values[0],1)) + " : " + str(passing['resid']) +\
                        "\nratio: " + str(round(row['ratio'].values[0],1)) + " : " + str(passing['ratio']) +\
                        "\ngw: " + str(round(row['GW'].values[0],6)) + " : " + str(passing['GW']) +\
                        "\nintensity: " + str(round(row['minIntensity'].values[0],1)) + " : " + str(passing['minIntensity']) +\
                        "\ncurrentCalc: " + str(round(row['currentCalc'].values[0],1))
                        
        if self.varLab:
            dataString = dataString + "\nFRC_NX: " + str(round(row['FRC_NX'].values[0],3)) + " : " + str(passing['FRC_NX'])
        self.PNGPlot.text(0.98, 0.45,dataString,
                          horizontalalignment='right',
                          verticalalignment='top',
                          transform = self.PNGPlot.transAxes,
                          color = stringColor)
        
        dataString = str(row['seqmod'].values[0]) + " : z=" + str(row['charge'].values[0])
        self.PNGPlot.text(0.02, 0.98,dataString,
                          horizontalalignment='left',
                          verticalalignment='top',
                          transform = self.PNGPlot.transAxes,
                          color = stringColor)
        self.PNGPlot.set_title(self.currentISOFile)
        vizLib.cleanAxis(self.PNGPlot, ticksOnly=True)
        """

    def calc_fit(self):
        plotsfile = self.datapath+self.currentISOFile+'.plots'
        txtfile = self.datapath+self.currentISOFile+'.txt'
        #confile = self.datapath+self.currentISOFile+'.newcon'
        
        df = pd.read_csv(plotsfile)
        df2 = pd.read_csv(txtfile, header=None, sep=' ')
        #df3 = pd.read_csv(confile, header=None, skiprows=6, sep=' ')

        del df['residAdj']
        self.PNGPlot.clear()
        #self.CONPlot.clear()
        #self.PNGPlot.plot(df2[0][0:len(df['dat'])].values, df['dat'].values, 'o', markersize=3, label='data', markerfacecolor='None', markeredgecolor='red')
        self.PNGPlot.plot(df2[0][0:len(df['dat'])].values, df['dat'].values, 'o', markersize=6, markerfacecolor='None', markeredgecolor='red')
        self.PNGPlot.plot(df2[0][0:len(df['dat'])].values, df['dat'].values, 'r-', linewidth=2, label='data')
        self.PNGPlot.plot(df2[0][0:len(df['fit'])].values, df['fit'].values, 'b-', linewidth=2, label='fit')
        self.PNGPlot.plot(df2[0][0:len(df['resid'])].values, df['resid'].values, 'g-', linewidth=2, label='residual')
        self.PNGPlot.set_xlabel('m/z')
        self.PNGPlot.set_ylabel('intensity')
        self.PNGPlot.set_xlim(df2[0][0:len(df['dat'])].values.min(), df2[0][0:len(df['dat'])].values.max())
        self.PNGPlot.legend()
        row = self.dataFrame[self.dataFrame['isofile'] == self.currentISOFile]
        passing = self.testPassRow(row)
        stringColor = 'black'
        if row['handDelete'].values[0] is True:
            stringColor = 'red'
        elif row['handSave'].values[0] is True:
            stringColor = 'green'
        dataString =    "ppmDiff: " + str(round(row['ppmDiff'].values[0],1)) + " : " + str(passing['ppmDiff']) +\
                        "\nN14: " + str(round(row['ppm_n14'].values[0],1)) + " : " + str(passing['N14']) +\
                        "\nN15: " + str(round(row['ppm_n15'].values[0],1)) + " : " + str(passing['N15']) +\
                        "\nmissed: " + str(row['missed'].values[0]) + " : " + str(passing['missed']) +\
                        "\nrtDiff: " + str(round(row['rtDiff'].values[0],3)) + " : " + str(passing['rtDiff']) +\
                        "\nresid: " + str(round(row['resid'].values[0],1)) + " : " + str(passing['resid']) +\
                        "\nratio: " + str(round(row['ratio'].values[0],1)) + " : " + str(passing['ratio']) +\
                        "\ngw: " + str(round(row['GW'].values[0],6)) + " : " + str(passing['GW']) +\
                        "\nintensity: " + str(round(row['minIntensity'].values[0],1)) + " : " + str(passing['minIntensity']) +\
                        "\ncurrentCalc: " + str(round(row['currentCalc'].values[0],1))
                        
        if self.varLab:
            dataString = dataString + "\nFRC_NX: " + str(round(row['FRC_NX'].values[0],3)) + " : " + str(passing['FRC_NX'])
        self.PNGPlot.text(0.98, 0.45,dataString,
                          horizontalalignment='right',
                          verticalalignment='top',
                          transform = self.PNGPlot.transAxes,
                          color = stringColor)
        
        dataString = str(row['seqmod'].values[0]) + " : z=" + str(row['charge'].values[0])
        self.PNGPlot.text(0.02, 0.98,dataString,
                          horizontalalignment='left',
                          verticalalignment='top',
                          transform = self.PNGPlot.transAxes,
                          color = stringColor)
        self.PNGPlot.set_title(self.currentISOFile)

        #conMax = df3.max().max()
        #for ypos in range(len(df3.index)):
        #    self.CONPlot.scatter(range(len(df3.columns)), [ypos]*len(df3.columns), s=20, alpha=0.5, c=df3.ix[ypos].values/conMax, cmap=pylab.cm.autumn_r, edgecolors='none')
        #self.CONPlot.set_xlim(0, len(df3.columns))
        #self.CONPlot.set_ylim(1, len(df3.index)-1)
        #self.CONPlot.set_xticks([])
        #self.CONPlot.set_yticks([])
        
        #vizLib.cleanAxis(self.CONPlot, complete=True)
        vizLib.cleanAxis(self.PNGPlot, ticksOnly=True)

    def calc_hist(self):
        name = self.currentHist
        dataAll = self.dataFrame[name]
        dataSelected=self.savedPoints[name]
        
        bin_num = min(30, len(list(set(dataAll))))
        self.histPlotAll.clear()
        self.histPlotAll.hist(dataAll, bins = bin_num)
        self.histPlotAll.text(0.05,0.75,name+'_All', transform=self.histPlotAll.transAxes)
        if self.currentHist == "minIntensity":
            self.histPlotAll.set_xscale('log')
        
        bin_num = min(30, len(list(set(dataSelected))))
        self.histPlotSelected.clear()
        self.histPlotSelected.hist(dataSelected.values, bins = bin_num)
        self.histPlotSelected.text(0.05,0.75,name+'_Selected', transform=self.histPlotSelected.transAxes)
        if self.currentHist == "minIntensity":
            self.histPlotSelected.set_xscale('log')
        vizLib.cleanAxis(self.histPlotSelected, ticksOnly=True)
        vizLib.cleanAxis(self.histPlotAll, ticksOnly=True)
    
    def draw_all(self):
        self.canvasLeft.draw()
        self.canvasRight.draw()
        
    def getChecksNum(self):
        toReturn = []        
        if self.lowCheckNum.IsChecked():
            toReturn.append('AMP_U')
        if self.pulse and self.midCheckNum.IsChecked():
            toReturn.append('AMP_L')
        if self.highCheckNum.IsChecked():
            toReturn.append('AMP_S')
        return toReturn
        
    def getChecksDen(self):
        toReturn = []        
        if self.lowCheckDen.IsChecked():
            toReturn.append('AMP_U')
        if self.pulse and self.midCheckDen.IsChecked():
            toReturn.append('AMP_L')
        if self.highCheckDen.IsChecked():
            toReturn.append('AMP_S')
        return toReturn

    def findRanges(self):    
        (self.ppmDiff_low, self.ppmDiff_high) = map(float, self.ppmDiffRangeBypass.GetValue().split(' '))
        (self.ppm_n14_low, self.ppm_n14_high) = map(float, self.N14RangeBypass.GetValue().split(' '))
        (self.ppm_n15_low, self.ppm_n15_high) = map(float, self.N15RangeBypass.GetValue().split(' '))
        (self.missed_low, self.missed_high) = map(float, self.missedRangeBypass.GetValue().split(' '))
        (self.rtDiff_low, self.rtDiff_high) = map(float, self.rtDiffRangeBypass.GetValue().split(' '))
        (self.resid_low, self.resid_high) = map(float, self.residRangeBypass.GetValue().split(' '))
        (self.ratio_low, self.ratio_high) = map(float, self.ratioLimBypass.GetValue().split(' '))
        (self.gw_low, self.gw_high) = map(float, self.gwBypass.GetValue().split(' '))
        self.minIntensity = float(self.minIntensityBypass.GetValue())
        if self.varLab:
            (self.FRC_NX_low, self.FRC_NX_high) = map(float, self.FRC_NXRangeBypass.GetValue().split(' '))

    def testPassRow(self, row):
        self.findRanges()
        passing = {'ppmDiff':True, 'N14':True, 'N15':True, 'missed':True, 'rtDiff':True, 'resid':True, 'varLab':True}
        
        passing['ppmDiff'] = (row['ppmDiff'].values[0] >= self.ppmDiff_low) and (row['ppmDiff'].values[0] <= self.ppmDiff_high) 
        passing['N14'] = (row['ppm_n14'].values[0] >= self.ppm_n14_low) and (row['ppm_n14'].values[0] <= self.ppm_n14_high)
        passing['N15'] = (row['ppm_n15'].values[0] >= self.ppm_n15_low) and (row['ppm_n15'].values[0] <= self.ppm_n15_high)
        passing['missed'] = (row['missed'].values[0] >= self.missed_low) and (row['missed'].values[0] <= self.missed_high)
        passing['rtDiff'] = (row['rtDiff'].values[0] >= self.rtDiff_low) and (row['rtDiff'].values[0] <= self.rtDiff_high)
        passing['resid'] = (row['resid'].values[0] >= self.resid_low) and (row['resid'].values[0] <= self.resid_high) 
        passing['minIntensity'] = (row['minIntensity'].values[0] >= self.minIntensity)
        passing['ratio'] = (row['ratio'].values[0] >= self.ratio_low) and (row['ratio'].values[0] <= self.ratio_high) 
        passing['GW'] = (row['GW'].values[0] >= self.gw_low) and (row['GW'].values[0] <= self.gw_high) 
        if self.varLab:
            passing['FRC_NX'] = (row['FRC_NX'].values[0] >= self.FRC_NX_low) and (row['FRC_NX'].values[0] <= self.FRC_NX_high)
        return passing

    def getPass(self, t):
        self.findRanges()
        filt = (self.dataFrame['missed'] >= -1)
        #filt = pd.DataFrame([True]*len(self.dataFrame.index))
        if self.ppmDiffOn.IsChecked():
            filt =  (self.dataFrame['ppmDiff'] >= self.ppmDiff_low) & (self.dataFrame['ppmDiff'] <= self.ppmDiff_high) 
        if self.N14On.IsChecked():
            filt = filt & (self.dataFrame['ppm_n14'] >= self.ppm_n14_low) & (self.dataFrame['ppm_n14'] <= self.ppm_n14_high)
        if self.N15On.IsChecked():
            filt = filt & (self.dataFrame['ppm_n15'] >= self.ppm_n15_low) & (self.dataFrame['ppm_n15'] <= self.ppm_n15_high)
        if self.missedOn.IsChecked():
            filt = filt & (self.dataFrame['missed'] >= self.missed_low) & (self.dataFrame['missed'] <= self.missed_high)
        if self.rtOn.IsChecked():
            filt = filt & (self.dataFrame['rtDiff'] >= self.rtDiff_low) & (self.dataFrame['rtDiff'] <= self.rtDiff_high)
        if self.residOn.IsChecked():
            filt = filt & (self.dataFrame['resid'] >= self.resid_low) & (self.dataFrame['resid'] <= self.resid_high)
        if self.ratioLimOn.IsChecked():
            filt = filt & (self.dataFrame['ratio'] >= self.ratio_low) & (self.dataFrame['ratio'] <= self.ratio_high)
        if self.gwOn.IsChecked():
            filt = filt & (self.dataFrame['GW'] >= self.gw_low) & (self.dataFrame['GW'] <= self.gw_high)
        if self.minIntensityOn.IsChecked():
            filt = filt & (self.dataFrame['minIntensity'] >= self.minIntensity)
        if self.varLab:
            if self.FRC_NXOn.IsChecked():
                filt = filt & (self.dataFrame['FRC_NX'] >= self.FRC_NX_low) & (self.dataFrame['FRC_NX'] <= self.FRC_NX_high)
        if self.handSaveOn.IsChecked():
            filt = filt | (self.dataFrame['handSave'] == True)
        if self.handDeleteOn.IsChecked():
            filt = filt & (self.dataFrame['handDelete'] == False)
        if t:
            return self.dataFrame[filt]
        else:
            return self.dataFrame[~filt]
            
    def pickScatterPoint(self, event):
        ind = event.ind
        thisline = event.artist
        ydata = thisline.get_ydata()
        xdata = thisline.get_xdata()
        self.currentRow = self.dataFrame[   (self.dataFrame['currentCalc'] == ydata[ind][0]) & \
                                                ((self.dataFrame['currentPos'] == xdata[ind][0]))]
        self.currentISOFile = self.currentRow['isofile'].values[0]
        self.selectedPoint.set_data(self.currentRow['currentPos'].values[0], self.currentRow['currentCalc'].values[0])
        self.canvasLeft.draw()
        self.calc_fit()
        self.canvasRight.draw()
        self.savedList.SetStringSelection(self.currentISOFile)
        self.filteredList.SetStringSelection(self.currentISOFile)
        
def setCurrentFrac(calcNum, calcDen):
    num = [i[-1:] for i in calcNum]
    den = [i[-1:] for i in calcDen]
    return string.join(num, "+") + " / " + string.join(den, "+")
    
def transformValue(inputData, transformData, toTransform, normalizeTo):
    inputNorm = float(inputData.loc[toTransform])/float(inputData.loc[normalizeTo])
    transformNorm = float(transformData.loc[toTransform])/float(transformData.loc[normalizeTo])
    inputNormOffset = inputNorm - transformNorm
    inputTransformed = inputNormOffset*float(inputData.loc[normalizeTo])
    return inputTransformed

def openFile(fullpath):
    r = csv.reader(open(fullpath))
    header = r.next()
    
    if not ('resid' in header and 'minIntensity' in header and 'ratio' in header and 'currentCalc' in header):
        print "preprocesing : " + fullpath + "..."
        dataFrame = qMS.preProcessIsoCSV(fullpath, True)
    else:
        dataFrame = qMS.readIsoCSV(fullpath, noProcess=False)

    puls = 'AMP_L' in dataFrame.columns
    vla = 'FRC_NX' in dataFrame.columns
    return [dataFrame, puls, vla]

def startApp(dataFrame, datapath, filename, pulse, varLab, fsize=None, size=None):
    app = wx.App()
    screenSize = wx.DisplaySize()
    if size is None:
        size = 1.0*float(screenSize[0])/1600.0
        fsize = 1.0*float(screenSize[0])/1600.0
        print "detected screen resolution as " + str(screenSize[0]) + "; setting size=fsize="+str(round(size,2))+", these can be changed at bottom of masse.py if desired."
    vizLib.setRcs(scale=size*12, tickScale=0.75*size, axisLineWidth=1)
    axisFont = {'titlesize' : size*12*1.2,
                'labelsize' : size*12*1.25}
    matplotlib.rc('axes', **axisFont)

    app.frame = MasseFrame(dataFrame, datapath, filename, pulse=pulse, varLab=varLab, fsize=fsize, size=size)
    app.frame.Show()
    app.MainLoop()    

def fileOpenStart(pathToFile=None):
        """
        Create and show the Open FileDialog
        """
        wildcard =  "All files (*.*)|*.*|"\
                    "Preprocessed _iso_res.csv file (*_iso_res.csv)|*_iso_res.csv|"\
                    "Massacre iso_csv file (*_iso.csv)|*_iso.csv|"
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
        [df, p, vl] = openFile(pathToFile)
        return [df, dp, fn, p, vl]
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
    [dfr, dpa, fna, pul, vlab] = fileOpenStart(pathToFile)
    
    startApp(dfr, dpa, fna, pul, vlab, fsize=fsize, size=size)