import os
import wx
import qMS
import qMSDefs
import string
import csv
import pandas
import numpy
import sys
import matplotlib.gridspec as gridspec
import vizLib
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar

FILENAMEHEADER = 'shortName'
UIDHEADER = 'UID'
ISOFILEHEADER = 'isofile'
PROTEINHEADER = 'protein'
FIGSIZE = 7.5

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
        self.priorFilters.SetValue(False)
        
        self.ppmDiffRangeBypass.SetValue('-100 100')
        self.N14RangeBypass.SetValue('-100 100')
        self.N15RangeBypass.SetValue('-100 100')
        self.residRangeBypass.SetValue('5000')
        self.ratioLimBypass.SetValue('0 1000')
        self.minIntLimBypass.SetValue('10')
        self.lastHistLimBypass.SetValue('-1 1')
        self.proteinZoomRangeBypass.SetValue('0 60')
        self.recalcAll()
        self.currentRow = self.dataFrame[self.dataFrame[UIDHEADER] == self.currentRow[UIDHEADER]].iloc[0]
        self.redrawAll(setZero=True)

    def create_main_panel(self, df, dp, fn, pulse=False, varLab=False, fsize=1.0, size=1.0):
        """ Creates the main panel with all the controls on it:
             * mpl canvas 
             * mpl navigation toolbar
             * Control panel for interaction
        """
        self.lastHistField = 'rtDiff'
        if varLab:
            self.lastHistField = 'FRC_NX'
        self.fsize = fsize
        self.size = size
        self.SetFont(wx.Font(10*size, wx.SWISS, wx.NORMAL, wx.NORMAL, False,'MS Shell Dlg 2'))
        self.calcNum = ["AMP_U"]
        self.calcDen = ["AMP_U", "AMP_S"]
        self.currentDirectory = os.getcwd()
        self.dataFrame = df
        self.dataFrame['minIntensity'].replace([-numpy.inf, numpy.inf], numpy.nan, inplace=True)
        self.dataFrame['minIntensity'].fillna(0, inplace=True)
        self.dataFrame['currentPosDataset'] = self.dataFrame['currentPos']
        frh = {'ppmDiff_low':-100, 'ppmDiff_high':100, 'ppm_n14_low':-100, 'ppm_n14_high':100,
               'ppm_n15_low':-100, 'ppm_n15_high':100, 'resid':100, 'ratio_low': 0, 'ratio_high':1000} 
        self.filterRangeHash = {i:frh for i in list(self.dataFrame[FILENAMEHEADER].unique())}
        positionFileDict = {key:int(value)+1 for value,key in enumerate(qMS.sort_nicely(sorted(self.dataFrame[FILENAMEHEADER].unique())))}
        positionFileLookup = pandas.Series(positionFileDict)
        self.dataFrame['currentPosProtein']=positionFileLookup[self.dataFrame[FILENAMEHEADER]].values

        self.positionLabelsDataset = [i[:3] for i in qMS.sort_nicely(sorted(set(self.dataFrame[PROTEINHEADER].values)))]
        self.positionLabelsProtein = [i.split('.')[0] for i in qMS.sort_nicely(sorted(set(self.dataFrame[FILENAMEHEADER].values)))]
        
        startDataset = qMS.sort_nicely(sorted(self.dataFrame[FILENAMEHEADER].unique()))[0]
        self.datasetView = self.dataFrame[self.dataFrame[FILENAMEHEADER] == startDataset]
        startUID = qMS.sort_nicely(list(self.datasetView[self.datasetView['allClear'] == True][UIDHEADER].values))[0]
        self.currentRow = self.datasetView[self.datasetView[UIDHEADER] == startUID].iloc[0]

        self.proteinView = self.dataFrame[self.dataFrame[PROTEINHEADER] == self.currentRow[PROTEINHEADER]]

        self.datapath = dp
        self.datafile = fn
        self.pulse=pulse
        self.varLab=varLab
        self.figdim = FIGSIZE*fsize
        
        self.panel = wx.Panel(self)
        
        '''Create the figure panels'''
        self.figLeft = Figure((self.figdim, self.figdim))
        self.canvasLeft = FigCanvas(self.panel, wx.ID_ANY, self.figLeft)
        self.figRight = Figure((self.figdim, self.figdim))
        self.canvasRight = FigCanvas(self.panel, wx.ID_ANY, self.figRight)
        
        gsLeft = gridspec.GridSpec(6,3)
        gsRight = gridspec.GridSpec(5,25)
        self.PLPlotDataset = self.figLeft.add_subplot(gsLeft[:4, :])
        self.hist1 = self.figLeft.add_subplot(gsLeft[4,0])
        self.hist2 = self.figLeft.add_subplot(gsLeft[4,1])
        self.hist3 = self.figLeft.add_subplot(gsLeft[5,0])
        self.hist4 = self.figLeft.add_subplot(gsLeft[5,1])
        self.hist5 = self.figLeft.add_subplot(gsLeft[4,2])
        self.hist6 = self.figLeft.add_subplot(gsLeft[5,2])
        self.PNGPlot = self.figRight.add_subplot(gsRight[:3, 1:])
        self.PLPlotProtein = self.figRight.add_subplot(gsRight[3:, :])
        
        '''Create the list boxes'''
        self.datasetsList = wx.ListBox(self.panel, id=26, choices=[], style=wx.LB_SINGLE, name='Dataset', size=(100*size,100*size))
        self.savedList = wx.ListBox(self.panel, id=26, choices=[], style=wx.LB_SINGLE, name='Saved', size=(100*size,100*size))
        self.filteredList = wx.ListBox(self.panel, id=26, choices=[], style=wx.LB_SINGLE, name='Filtered', size=(100*size,100*size))
        
        '''Create the buttons'''
        self.toolbarLeft = NavigationToolbar(self.canvasLeft)
        self.toolbarRight = NavigationToolbar(self.canvasRight)
        
        self.cb_grid = wx.CheckBox(self.panel, wx.ID_ANY, label="Grid", style=wx.ALIGN_RIGHT)        
        self.hideCheck = wx.CheckBox(self.panel, wx.ID_ANY, label="Hide", style=wx.ALIGN_RIGHT)        
        self.zoomCheck = wx.CheckBox(self.panel, wx.ID_ANY, label="Zoom", style=wx.ALIGN_RIGHT)        
        self.priorFilters = wx.CheckBox(self.panel, wx.ID_ANY, label="oFilt")        
        self.exportButton = wx.Button(self.panel, wx.ID_ANY, "Export")
        self.openButton = wx.Button(self.panel, wx.ID_ANY, "Open")
                
        self.calcButton = wx.Button(self.panel, wx.ID_ANY, "Calc", size=(57*size,40))
        
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
        
        self.handSaveOn = wx.CheckBox(self.panel, wx.ID_ANY, label="save")
        self.handDeleteOn = wx.CheckBox(self.panel, wx.ID_ANY, label="delete")

        
        self.ppmDiff_range_button = wx.Button(self.panel, wx.ID_ANY, "PPM diff",size=(75*size,-1))
        self.N14_range_button = wx.Button(self.panel, wx.ID_ANY, "14N PPM", size=(75*size,-1))
        self.N15_range_button = wx.Button(self.panel, wx.ID_ANY, "15N PPM", size=(75*size,-1))
        self.resid_range_button = wx.Button(self.panel, wx.ID_ANY, "resid.", size=(70*size,-1))
        self.ratio_range_button = wx.Button(self.panel, wx.ID_ANY, "cc", size=(80*size,-1))
        self.minInt_range_button = wx.Button(self.panel, wx.ID_ANY, "minInt", size=(80*size,-1))
        self.lastHist_range_button = wx.Button(self.panel, wx.ID_ANY, self.lastHistField, size=(80*size,-1))

        self.ppmDiffRangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.N14RangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.N15RangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.residRangeBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.ratioLimBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.minIntLimBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        self.lastHistLimBypass = wx.TextCtrl(self.panel, size=(75*size,-1),style=wx.TE_PROCESS_ENTER)
        
        '''lay out the buttons'''
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.figbox = wx.BoxSizer(wx.HORIZONTAL)
        self.figbox.Add(self.canvasLeft, 0)
        self.listBox = wx.BoxSizer(wx.VERTICAL)
        self.datasetsListBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'datasets'), wx.HORIZONTAL)
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
                
        #add calculate button
        self.calcNBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'calcDist'), wx.VERTICAL)
        self.calcNBox.Add(self.calcButton, flag=wx.ALIGN_TOP)
        self.calcNBox.Add(self.priorFilters, 0, flag=wx.ALIGN_LEFT)
        self.toolNumBox.Add(self.calcNBox, 0, flag=wx.ALIGN_RIGHT | wx.GROW)
        
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
        
        self.zoomBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'zoom tools'), wx.VERTICAL)
        self.zoomBox.Add(self.proteinZoomOn, 0, flag=wx.ALIGN_LEFT)
        self.zoomBox.Add(self.proteinZoomRangeBypass, 0, flag=wx.ALIGN_LEFT)

        self.handBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,'hand tools'), wx.VERTICAL)
        self.handBox.Add(self.handSaveOn, 0, flag=wx.ALIGN_LEFT)
        self.handBox.Add(self.handDeleteOn, 0, flag=wx.ALIGN_LEFT)
        

        self.expandBox = wx.BoxSizer(wx.HORIZONTAL)
        self.toolNumBox.Add(self.handBox, 0, flag=wx.ALIGN_LEFT)
        self.toolNumBox.Add(self.zoomBox, 0, flag=wx.ALIGN_LEFT)
        self.toolNumBox.Add(self.expandBox, 1, flag=wx.GROW)
        self.toolNumBox.Add(self.imageToolsBoxR, 0, flag=wx.ALIGN_RIGHT | wx.GROW)
        
        #add hbox to vbox
        self.vbox.Add(self.toolNumBox, 0, flag = wx.GROW)
        
        # Sliders for setting the various cutoffs
        self.filterBox = wx.StaticBoxSizer(wx.StaticBox(self.panel,wx.ID_ANY,), wx.HORIZONTAL)
        flags = wx.ALIGN_LEFT | wx.ALIGN_CENTER
        self.filterBox.Add(self.ppmDiff_range_button, 0, flag=flags)
        self.filterBox.Add(self.ppmDiffRangeBypass, 0, flag=flags)
        self.filterBox.Add(self.N14_range_button, 0, flag=flags)
        self.filterBox.Add(self.N14RangeBypass, 0, flag=flags)
        self.filterBox.Add(self.N15_range_button, 0, flag=flags)
        self.filterBox.Add(self.N15RangeBypass, 0, flag = flags)
        self.filterBox.Add(self.resid_range_button, 0, flag=flags)
        self.filterBox.Add(self.residRangeBypass, 0, flag = flags)
        self.filterBox.Add(self.ratio_range_button, 0, flag=flags)
        self.filterBox.Add(self.ratioLimBypass, 0, flag = flags)
        self.filterBox.Add(self.minInt_range_button, 0, flag=flags)
        self.filterBox.Add(self.minIntLimBypass, 0, flag = flags)
        self.filterBox.Add(self.lastHist_range_button, 0, flag=flags)
        self.filterBox.Add(self.lastHistLimBypass, 0, flag = flags)
        self.vbox.Add(self.filterBox, 0, flag = wx.ALIGN_LEFT | wx.TOP)

        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)
        
        '''bind events for buttons'''
        self.canvasLeft.mpl_connect('pick_event', self.pickScatterPointDataset)
        self.canvasRight.mpl_connect('pick_event', self.pickScatterPointProtein)

        self.exportButton.Bind(wx.EVT_BUTTON, self.on_export_button)
        self.openButton.Bind(wx.EVT_BUTTON, self.on_open_button)
        self.calcButton.Bind(wx.EVT_BUTTON, self.on_recalc)
        self.cb_grid.Bind(wx.EVT_CHECKBOX, self.on_redraw)
        self.hideCheck.Bind(wx.EVT_CHECKBOX, self.on_redraw)
        self.zoomCheck.Bind(wx.EVT_CHECKBOX, self.on_redraw)
        self.proteinZoomOn.Bind(wx.EVT_CHECKBOX, self.on_redraw)
        self.priorFilters.Bind(wx.EVT_CHECKBOX, self.on_recalc)
        
        self.ppmDiff_range_button.Bind(wx.EVT_BUTTON, self.on_filt_button)
        self.N14_range_button.Bind(wx.EVT_BUTTON, self.on_filt_button)
        self.N15_range_button.Bind(wx.EVT_BUTTON, self.on_filt_button)
        self.resid_range_button.Bind(wx.EVT_BUTTON, self.on_filt_button)
        self.ratio_range_button.Bind(wx.EVT_BUTTON, self.on_filt_button)
        self.minInt_range_button.Bind(wx.EVT_BUTTON, self.on_filt_button)
        self.lastHist_range_button.Bind(wx.EVT_BUTTON, self.on_filt_button)
        self.handSaveOn.Bind(wx.EVT_CHECKBOX, self.on_filt_button)
        self.handDeleteOn.Bind(wx.EVT_CHECKBOX, self.on_filt_button)
        self.proteinZoomOn.Bind(wx.EVT_CHECKBOX, self.on_pzoom_check)
        
        self.datasetsList.Bind(wx.EVT_LISTBOX, self.on_datasetsBoxClick)        
        self.savedList.Bind(wx.EVT_LISTBOX, self.on_savedBoxClick)
        self.filteredList.Bind(wx.EVT_LISTBOX, self.on_filteredBoxClick)
        
        self.panel.Bind(wx.EVT_KEY_UP, self.on_key_press)

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
            self.recalcAll()
            self.redrawAll()
            self.dataFrame['priorFilter'] = self.dataFrame['allFPass']
            self.dataFrame.to_csv(path, index=False)
            summaryCSVPath = path.split('.')[0] + '_median_[' + ''.join(self.calcNum) + ']_[' + ''.join(self.calcDen) + '].csv'
            self.writeSummaryCSV(summaryCSVPath)
            
        dlg.Destroy()

    def writeSummaryCSV(self, path):
        csvDict = {}
        proteins = qMS.sort_nicely(sorted(set(self.dataFrame[PROTEINHEADER].values)))
        fractions = qMS.sort_nicely(sorted(set(self.dataFrame[FILENAMEHEADER].values)))
        for frac in fractions:
            csvDict[frac] = pandas.Series([self.dataFrame[(self.dataFrame[FILENAMEHEADER] == frac) &\
                                            (self.dataFrame[PROTEINHEADER] == i) &\
                                            (self.dataFrame['allClear'] == True)]['currentCalc'].median() for i in proteins], index=proteins)
        summaryDF = pandas.DataFrame(csvDict, columns = fractions)
        summaryDF.to_csv(path)
    def on_exit(self, event):
        self.Destroy()
    
    def on_redraw(self, event):
        self.redrawAll()
    
    def on_pzoom_check(self, event):
        self.recalcAll()
        self.redrawAll()
    
    def on_filt_button(self, event):
        self.calc_filters()
        self.redrawAll()
    
    def on_recalc(self, event):
        self.recalcAll()
        self.redrawAll()

    def redrawAll(self, setZero=False):
        self.calc_figureLeft()
        self.calc_hists()
        self.calc_figureRightFit()
        self.calc_figureRightProtein()
        self.calc_lists()
        self.draw_all()
        if setZero:
            self.figLeft.tight_layout()
            self.figRight.tight_layout()
            for i in self.dataFrame[FILENAMEHEADER].unique():
                self.filterRangeHash[i] = self.findRanges()
            
        
    def recalcAll(self):
        self.calc_data()
        self.calc_filters()

    def draw_all(self):
        self.canvasLeft.draw()
        self.canvasRight.draw()
    
    def resetFilterRanges(self):
        ha = self.filterRangeHash[self.currentRow[FILENAMEHEADER]]
        self.ppmDiffRangeBypass.SetValue(str(ha['ppmDiff_low']) + ' ' + str(ha['ppmDiff_high']))
        self.N14RangeBypass.SetValue(str(ha['ppm_n14_low']) + ' ' + str(ha['ppm_n14_high']))
        self.N15RangeBypass.SetValue(str(ha['ppm_n15_low']) + ' ' + str(ha['ppm_n15_high']))
        self.residRangeBypass.SetValue(str(ha['resid']))
        self.ratioLimBypass.SetValue(str(ha['ratio_low']) + ' ' + str(ha['ratio_high']))
        self.minIntLimBypass.SetValue(str(ha['minInt']))
        self.lastHistLimBypass.SetValue(str(ha[self.lastHistField+'_low']) + ' ' + str(ha[self.lastHistField+'_high']))
        
    def on_datasetsBoxClick(self, event):
        if not (self.datasetsList.GetStringSelection() is u''):
            pickedDataset = self.datasetsList.GetStringSelection()
            pickedUID = qMS.sort_nicely(list(self.dataFrame[(self.dataFrame[FILENAMEHEADER] == pickedDataset) & \
                                                       (self.dataFrame['allClear'] == True)][UIDHEADER].values))[0]            
            self.currentRow = self.dataFrame[self.dataFrame[UIDHEADER] == pickedUID].iloc[0]
            self.resetFilterRanges()
            self.calc_filters()
            self.redrawAll()

    def on_savedBoxClick(self, event):
        if not (self.savedList.GetStringSelection() is u''):
            pickedUID = self.savedList.GetStringSelection()
            self.currentRow = self.dataFrame[self.dataFrame[UIDHEADER] == pickedUID].iloc[0]
            self.newPepSelection()
        
    def on_filteredBoxClick(self, event):
        if not (self.filteredList.GetStringSelection() is u''):
            pickedUID = self.filteredList.GetStringSelection()
            self.currentRow = self.dataFrame[self.dataFrame[UIDHEADER] == pickedUID].iloc[0]
            self.newPepSelection()
    
    def newPepSelection(self):
        self.calc_figureRightFit()
        self.calc_figureRightProtein()
        self.canvasRight.draw()
        self.selectedPointDataset.set_data(self.currentRow['currentPosDataset'], self.currentRow['currentCalc'])
        self.canvasLeft.draw()

    def on_key_press(self, event):
        event.Skip()
        c = event.GetKeyCode()
        if c is 68: #got a d keystroke
            self.dataFrame.loc[self.dataFrame[UIDHEADER] == self.currentRow[UIDHEADER],'handDelete'] = True
            self.dataFrame.loc[self.dataFrame[UIDHEADER] == self.currentRow[UIDHEADER],'handSave'] = False
            try:            
                myIndex = list(self.savedListItems).index(self.currentRow[UIDHEADER])
            except ValueError:
                myIndex = 0
            try:
                nextItem = self.savedListItems[myIndex+1]
            except IndexError:
                nextItem = self.savedListItems[myIndex-1]
            self.currentRow = self.dataFrame[self.dataFrame[UIDHEADER] == nextItem].iloc[0]
            self.calc_hand_filters()
            self.redrawAll()
            
        elif c is 83: # got a s keystroke
            self.dataFrame.loc[self.dataFrame[UIDHEADER] == self.currentRow[UIDHEADER],'handDelete'] = False
            self.dataFrame.loc[self.dataFrame[UIDHEADER] == self.currentRow[UIDHEADER],'handSave'] = True
            self.calc_hand_filters()
            self.redrawAll()
            
        elif c is 69: # got a e keystroke (delete all)
            self.dataFrame.loc[self.dataFrame['TID'] == self.currentRow['TID'], 'handDelete'] = True
            self.dataFrame.loc[self.dataFrame['TID'] == self.currentRow['TID'], 'handSave'] = False
            try:
                myIndex = list(self.savedListItems).index(self.currentRow[UIDHEADER])
            except ValueError:
                myIndex = 0
            try:
                nextItem = self.savedListItems[myIndex+1]
            except IndexError:
                nextItem = self.savedListItems[myIndex-1]
            self.currentRow = self.dataFrame[self.dataFrame[UIDHEADER] == nextItem].iloc[0]
            self.calc_hand_filters()
            self.redrawAll()
            
        elif c is 87: # got a w keystroke (save all)
            self.dataFrame.loc[self.dataFrame['TID'] == self.currentRow['TID'], 'handDelete'] = False
            self.dataFrame.loc[self.dataFrame['TID'] == self.currentRow['TID'], 'handSave'] = True
            self.calc_hand_filters()
            self.redrawAll()
    
    def pickScatterPointDataset(self, event):
        ind = event.ind
        thisline = event.artist
        ydata = thisline.get_ydata()
        xdata = thisline.get_xdata()
        self.currentRow = self.dataFrame[(self.dataFrame[FILENAMEHEADER] == self.currentRow[FILENAMEHEADER]) &\
                                            (self.dataFrame['currentCalc'] == ydata[ind][0]) & \
                                            ((self.dataFrame['currentPosDataset'] == xdata[ind][0]))].iloc[0]
        self.selectedPointDataset.set_data(self.currentRow['currentPosDataset'], self.currentRow['currentCalc'])
        self.canvasLeft.draw()
        self.calc_figureRightFit()
        self.calc_figureRightProtein()
        self.canvasRight.draw()
        self.savedList.SetStringSelection(self.currentRow[UIDHEADER])
        self.filteredList.SetStringSelection(self.currentRow[UIDHEADER])
        
    def pickScatterPointProtein(self, event):
        ind = event.ind
        thisline = event.artist
        ydata = thisline.get_ydata()
        xdata = thisline.get_xdata()
        self.currentRow = self.dataFrame[(self.dataFrame[PROTEINHEADER] == self.currentRow[PROTEINHEADER]) & \
                                            (self.dataFrame['currentCalc'] == ydata[ind][0]) & \
                                            ((self.dataFrame['currentPosProtein'] == xdata[ind][0]))].iloc[0]
        self.selectedPointProtein.set_data(self.currentRow['currentPosProtein'], self.currentRow['currentCalc'])
        theRows = self.dataFrame[(self.dataFrame[PROTEINHEADER] == self.currentRow[PROTEINHEADER]) & \
                                    (self.dataFrame['TID'] == self.currentRow['TID'])]
        self.selectedPairedPointsProtein.set_data(theRows['currentPosProtein'], theRows['currentCalc'])
        self.resetFilterRanges()
        self.calc_filters()
        self.redrawAll()

    def calc_lists(self):
        self.savedListItems = numpy.array(qMS.sort_nicely(list(self.dataFrame[(self.dataFrame[FILENAMEHEADER] == self.currentRow[FILENAMEHEADER]) & \
                                                                (self.dataFrame['allClear'] == True)][UIDHEADER].values)))
        self.filteredListItems = numpy.array(qMS.sort_nicely(list(self.dataFrame[(self.dataFrame[FILENAMEHEADER] == self.currentRow[FILENAMEHEADER]) & \
                                                                    (self.dataFrame['allClear'] == False)][UIDHEADER].values)))
        self.savedList.Set(self.savedListItems)
        self.filteredList.Set(self.filteredListItems)
        self.savedList.SetStringSelection(self.currentRow[UIDHEADER])
        self.filteredList.SetStringSelection(self.currentRow[UIDHEADER])
        self.datasetsListItems = qMS.sort_nicely(list(self.dataFrame[FILENAMEHEADER].unique()))
        self.datasetsList.Set(self.datasetsListItems)
        self.datasetsList.SetStringSelection(self.currentRow[FILENAMEHEADER])
    
    def calc_data(self):
        self.calcNum = self.getChecksNum()
        self.calcDen = self.getChecksDen()
        self.dataFrame['currentCalc'] = qMS.calcValue(self.dataFrame, self.calcNum, self.calcDen)
    
    def calc_filters(self):
        rangeHash = self.findRanges()
        curFile = self.currentRow[FILENAMEHEADER]
        if not self.priorFilters.IsChecked():
            self.dataFrame.loc[self.dataFrame[FILENAMEHEADER] == curFile, 'ppmDiff_pass'] = (self.dataFrame['ppmDiff'] >= rangeHash['ppmDiff_low']) & (self.dataFrame['ppmDiff'] <= rangeHash['ppmDiff_high'])
            self.dataFrame.loc[self.dataFrame[FILENAMEHEADER] == curFile, 'ppm_n14_pass'] = (self.dataFrame['ppm_n14'] >= rangeHash['ppm_n14_low']) & (self.dataFrame['ppm_n14'] <= rangeHash['ppm_n14_high'])
            self.dataFrame.loc[self.dataFrame[FILENAMEHEADER] == curFile, 'ppm_n15_pass'] = (self.dataFrame['ppm_n15'] >= rangeHash['ppm_n15_low']) & (self.dataFrame['ppm_n15'] <= rangeHash['ppm_n15_high'])
            self.dataFrame.loc[self.dataFrame[FILENAMEHEADER] == curFile, 'resid_pass'] = self.dataFrame['resid'] <= rangeHash['resid']
            self.dataFrame.loc[self.dataFrame[FILENAMEHEADER] == curFile, 'ratio_pass'] = (self.dataFrame['currentCalc'] >= rangeHash['ratio_low']) & (self.dataFrame['currentCalc'] <= rangeHash['ratio_high'])
            self.dataFrame.loc[self.dataFrame[FILENAMEHEADER] == curFile, 'minInt_pass'] = self.dataFrame['minIntensity'] >= rangeHash['minInt']
            self.dataFrame.loc[self.dataFrame[FILENAMEHEADER] == curFile, self.lastHistField+'_pass'] = (self.dataFrame[self.lastHistField] >= rangeHash[self.lastHistField+'_low']) & \
                                                                                                        (self.dataFrame[self.lastHistField] <= rangeHash[self.lastHistField+'_high'])

            self.dataFrame.loc[self.dataFrame[FILENAMEHEADER] == self.currentRow[FILENAMEHEADER], 'allFPass'] = self.dataFrame['ppmDiff_pass'] & \
                                                                                                                self.dataFrame['ppm_n14_pass'] & \
                                                                                                                self.dataFrame['ppm_n15_pass'] & \
                                                                                                                self.dataFrame['resid_pass'] & \
                                                                                                                self.dataFrame['ratio_pass'] & \
                                                                                                                self.dataFrame['minInt_pass'] & \
                                                                                                                self.dataFrame[self.lastHistField+'_pass']
        else:
            self.dataFrame.loc[self.dataFrame[FILENAMEHEADER] == self.currentRow[FILENAMEHEADER], 'allFPass'] = self.dataFrame['priorFilter']
    
        self.dataFrame.loc[self.dataFrame[FILENAMEHEADER] == self.currentRow[FILENAMEHEADER], 'allClear'] = self.dataFrame['allFPass']
        self.calc_hand_filters()
    def calc_hand_filters(self):
        if self.handDeleteOn.IsChecked():
            self.dataFrame['allClear'] = (self.dataFrame['allClear']) & (~self.dataFrame['handDelete'])
        if self.handSaveOn.IsChecked():
            self.dataFrame['allClear'] = (self.dataFrame['allClear']) | (self.dataFrame['handSave'])
    def findRanges(self):
        currHash = {}
        (currHash['ppmDiff_low'], currHash['ppmDiff_high']) = map(float, self.ppmDiffRangeBypass.GetValue().split(' '))
        (currHash['ppm_n14_low'], currHash['ppm_n14_high']) = map(float, self.N14RangeBypass.GetValue().split(' '))
        (currHash['ppm_n15_low'], currHash['ppm_n15_high']) = map(float, self.N15RangeBypass.GetValue().split(' '))
        (currHash['ratio_low'], currHash['ratio_high']) = map(float, self.ratioLimBypass.GetValue().split(' '))
        currHash['resid'] = float(self.residRangeBypass.GetValue())
        currHash['minInt'] = float(self.minIntLimBypass.GetValue())
        (currHash[self.lastHistField+'_low'], currHash[self.lastHistField+'_high']) = map(float, self.lastHistLimBypass.GetValue().split(' '))
        self.filterRangeHash[self.currentRow[FILENAMEHEADER]] = currHash.copy()
        return currHash

    def determineMedians(self, view, field):
        xs = list(set(view[field]))
        ys = [view[view[field]==k]['currentCalc'].median() for k in xs]
        xs = [float(i) for i in xs]
        return [xs,ys]
       
    def calc_figureLeft(self):
        self.PLPlotDataset.clear()
        self.datasetView = self.dataFrame[self.dataFrame[FILENAMEHEADER] == self.currentRow[FILENAMEHEADER]]
        self.passDatasetView = self.datasetView[self.datasetView['allClear'] == True]
        self.failDatasetView = self.datasetView[self.datasetView['allClear'] == False]

        self.selectedPointDataset, = self.PLPlotDataset.plot(self.currentRow['currentPosDataset'], self.currentRow['currentCalc'], 
                                              'o', ms=20, alpha=0.5, color='yellow', visible=True)
        self.PLPlotDataset.grid(self.cb_grid.IsChecked())
        self.PLPlotDataset.plot(self.passDatasetView['currentPosDataset'], self.passDatasetView['currentCalc'], 'ro', picker=5, label="Saved : " + str(len(self.passDatasetView['currentCalc'].values)))
        meds = self.determineMedians(self.passDatasetView, 'currentPosDataset')
        self.PLPlotDataset.plot(meds[0], meds[1], 'g-', lw=2, label="Median : " + str(round(numpy.median(meds[1]),1)))
        if not self.hideCheck.IsChecked():
            self.PLPlotDataset.plot(self.failDatasetView['currentPosDataset'], self.failDatasetView['currentCalc'], 'x', mec='grey', picker=5, label="Filtered : " + str(len(self.failDatasetView['currentCalc'].values)))
        self.PLPlotDataset.set_xticks(range(1,int(self.dataFrame['currentPosDataset'].max())+1))
        self.PLPlotDataset.set_xticklabels(self.positionLabelsDataset, rotation=90, size='small')
        self.PLPlotDataset.set_title(self.currentRow[FILENAMEHEADER] + " : " + setCurrentFrac(self.calcNum, self.calcDen))
        if self.proteinZoomOn.IsChecked():
            (x_low, x_high) = map(float, self.proteinZoomRangeBypass.GetValue().split(' '))
            self.PLPlotDataset.set_xlim([x_low,x_high])
        else:
            self.PLPlotDataset.set_xlim([0.5,self.datasetView['currentPosDataset'].max()+1.5])        
        
        if not self.zoomCheck.IsChecked():
            self.PLPlotDataset.set_ylim([0,max(self.passDatasetView['currentCalc'].max(),1)])
        else:
            self.PLPlotDataset.set_ylim([0,3])
        self.PLPlotDataset.legend()
        vizLib.cleanAxis(self.PLPlotDataset, ticksOnly=True)
    
    def calc_hists(self):
        hists = {'ppmDiff':self.hist1, 'ppm_n14':self.hist2, 'ppm_n15':self.hist3, 'resid':self.hist4, 'minIntensity':self.hist5, self.lastHistField:self.hist6}
        relSet = self.dataFrame[(self.dataFrame[FILENAMEHEADER] == self.currentRow[FILENAMEHEADER]) & (self.dataFrame['allClear'] == True)]      
        for curHist in hists.keys():
            hists[curHist].clear()
            data = relSet[curHist].values
            bin_num = min(20, len(list(set(data))))
            if curHist == 'minIntensity' or curHist == 'resid':
                bin_num = numpy.logspace(numpy.log10(data.min()), numpy.log10(data.max()), num=20)
            hists[curHist].hist(data, bins = bin_num)
            hists[curHist].text(0.05,0.75,curHist, transform=hists[curHist].transAxes)
            numTicks = 2
            spacing = (data.max()-data.min())/(numTicks+2)
            if curHist == 'minIntensity' or curHist == 'resid':
                hists[curHist].set_xscale('log')
            else:
                hists[curHist].set_xticks([data.min()+spacing*i for i in range(numTicks+3)])
                hists[curHist].set_xticklabels([str(round(data.min()+spacing*i,1)) for i in range(numTicks+3)], size='small')
                hists[curHist].set_xlim(data.min(), data.max())
            hists[curHist].set_yticks([])
            vizLib.cleanAxis(hists[curHist])
    
    def calc_figureRightFit(self):
        plotsfile = self.datapath+self.currentRow[ISOFILEHEADER]+'.plots'
        txtfile = self.datapath+self.currentRow[ISOFILEHEADER]+'.txt'
        df = pandas.read_csv(plotsfile)
        df2 = pandas.read_csv(txtfile, header=None, sep=' ')
        del df['residAdj']
        self.PNGPlot.clear()
        self.PNGPlot.plot(df2[0][0:len(df['dat'])].values, df['dat'].values, 'o', markersize=6, markerfacecolor='None', markeredgecolor='red')
        self.PNGPlot.plot(df2[0][0:len(df['dat'])].values, df['dat'].values, 'r-', linewidth=2, label='data')
        self.PNGPlot.plot(df2[0][0:len(df['fit'])].values, df['fit'].values, 'b-', linewidth=2, label='fit')
        self.PNGPlot.plot(df2[0][0:len(df['resid'])].values, df['resid'].values, 'g-', linewidth=2, label='residual')
        self.PNGPlot.set_xlim(df2[0][0:len(df['dat'])].values.min(), df2[0][0:len(df['dat'])].values.max())
        self.PNGPlot.legend()
        stringColor = 'black'
        row = self.currentRow
        if row['handDelete']:
            stringColor = 'red'
        elif row['handSave']:
            stringColor = 'green'
        seqText = str(self.currentRow['seqmod']) + " : z=" + str(self.currentRow['charge'])

        dataString1 =    "ppmDiff: " + str(round(row['ppmDiff'],1)) + " : " + str(row['ppmDiff_pass']) +\
                        "\nN14: " + str(round(row['ppm_n14'],1)) + " : " + str(row['ppm_n14_pass']) +\
                        "\nN15: " + str(round(row['ppm_n15'],1)) + " : " + str(row['ppm_n15_pass']) +\
                        "\nresid: " + str(round(row['resid'],1)) + " : " + str(row['resid_pass'])
        dataString2 =  "ratio: " + str(round(row['currentCalc'],1)) + " : " + str(row['ratio_pass']) +\
                        "\nintensity: " + str(round(row['minIntensity'],1)) + " : " + str(row['minInt_pass']) +\
                        "\n"+self.lastHistField+" "+ str(round(row[self.lastHistField],1)) + " : " + str(row[self.lastHistField+'_pass']) +\
                        "\ncurrentCalc: " + str(round(row['currentCalc'],1))
        self.PNGPlot.text(0.02, 0.98,seqText,
                          horizontalalignment='left',
                          verticalalignment='top',
                          transform = self.PNGPlot.transAxes,
                          color = stringColor,
                          size='large')        
        
        self.PNGPlot.text(0.98, 0.02,dataString1,
                          horizontalalignment='right',
                          verticalalignment='bottom',
                          transform = self.PNGPlot.transAxes,
                          color = stringColor)
        self.PNGPlot.text(0.02, 0.02,dataString2,
                          horizontalalignment='left',
                          verticalalignment='bottom',
                          transform = self.PNGPlot.transAxes,
                          color = stringColor)
            
        self.PNGPlot.set_title(self.currentRow[UIDHEADER])
        self.PNGPlot.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.0e'))
        vizLib.cleanAxis(self.PNGPlot, ticksOnly=True)
    
    def calc_figureRightProtein(self):
        self.PLPlotProtein.clear()
        self.proteinView = self.dataFrame[self.dataFrame[PROTEINHEADER] == self.currentRow[PROTEINHEADER]]
        self.passProteinView = self.proteinView[self.proteinView['allClear'] == True]
        self.failProteinView = self.proteinView[self.proteinView['allClear'] == False]

        theRows = self.proteinView[self.proteinView['TID'] == self.currentRow['TID']]
        self.selectedPairedPointsProtein, = self.PLPlotProtein.plot(theRows['currentPosProtein'], theRows['currentCalc'], 'o', ms=30*self.size, alpha=0.25, color='red')
        self.selectedPointProtein, = self.PLPlotProtein.plot(self.currentRow['currentPosProtein'], self.currentRow['currentCalc'], 
                                              'o', ms=20, alpha=0.5, color='yellow', visible=True)
        
        self.PLPlotProtein.plot(self.passProteinView['currentPosProtein'], self.passProteinView['currentCalc'], 'ro', picker=5, label="Saved : " + str(len(self.passProteinView['currentCalc'].values)))
        meds = self.determineMedians(self.passProteinView, 'currentPosProtein')
        self.PLPlotProtein.plot(meds[0], meds[1], 'g-', lw=2, label="Median : " + str(round(numpy.median(meds[1]),1)))
        if not self.hideCheck.IsChecked():
            self.PLPlotProtein.plot(self.failProteinView['currentPosProtein'], self.failProteinView['currentCalc'], 'x', mec='grey', picker=5, label="Filtered : " + str(len(self.failProteinView['currentCalc'].values)))
        self.PLPlotProtein.grid(self.cb_grid.IsChecked())        
        self.PLPlotProtein.set_xticks(range(1,int(self.dataFrame['currentPosProtein'].max())+1))
        self.PLPlotProtein.set_xticklabels(self.positionLabelsProtein, rotation=30, size='small')
        self.PLPlotProtein.set_title(self.currentRow[PROTEINHEADER] + " : " + setCurrentFrac(self.calcNum, self.calcDen))
        self.PLPlotProtein.set_xlim([0,self.dataFrame['currentPosProtein'].max()+1])
        if not self.zoomCheck.IsChecked():
            self.PLPlotProtein.set_ylim([0,max(self.passProteinView['currentCalc'].max(),1)])
        else:
            self.PLPlotProtein.set_ylim([0,3])
        self.PLPlotProtein.legend()
        vizLib.cleanAxis(self.PLPlotProtein, ticksOnly=True)
        
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
    vizLib.setRcs(scale=size*12, tickScale=0.5*size, axisLineWidth=1)
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
                    "Merged and preprocessed _iso_res.csv file (*_iso_res.csv)|*_iso_res.csv|"
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
    #pathToFile = '/home/jhdavis/data/compiledDeps/L17Unfilt/L17_1FractionListMergedNew_iso_res.csv'
    fsize=None
    size=None
    [dfr, dpa, fna, pul, vlab] = fileOpenStart(pathToFile)
    
    startApp(dfr, dpa, fna, pul, vlab, fsize=fsize, size=size)