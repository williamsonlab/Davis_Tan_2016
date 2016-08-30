import math
import matplotlib.pyplot as pylab
import scipy.cluster.hierarchy as sch
import qMS
import numpy
import matplotlib
import brewer2mpl
import mrmTools
#from mpl_toolkits.mplot3d.axes3d import Axes3D
from IPython.core.display import HTML
from mpl_toolkits.mplot3d import Axes3D

def plotMRMScatterPlot(df, samples=None, labelHash=None, yAxisLabel=None, proteins=None, alpha=1.0, legend=True, legendCols=5, median=True,
                       grid=True, yAxis=None, num=['light '], den=['heavy '], yTicks=4, sampleLocations=None, sampleField='File Name', ax=None,
                       medianMarkerInc=1.25, colors=None, figSize=(10,10), markersize=10, fitLine=False, fitR2=False, medianColor='black',
                       yMin=0, yMax=10, xMin=None, xMax=None, title=None, zOrder=None, medianOnly=None, proteinNameHeader = 'Protein Name'):
    df.loc[:,'currentCalc'] = mrmTools.calcValue(df, num, den)
    if samples is None:
        samples = qMS.sort_nicely(list(df[sampleField].unique()))
    if proteins is None:
        proteins = qMS.sort_nicely(list(df[proteinNameHeader].unique()))
    if colors is None:
        colors = pylab.cm.jet([float(i)/float(len(proteins)) for i in range(len(proteins))])
    if yAxis is None:
        interval = (yMax*1.0-yMin)/(yTicks+1.0)
        yAxis = [yMin+i*interval for i in range(0,yTicks+2)]
    if labelHash is None:
        labelHash = {i:i for i in proteins}
    if sampleLocations is None:
        sampleLocations = {samples[i]:i for i in range(len(samples))}
    if xMin is None:
        xMin = min(sampleLocations.values())
    if xMax is None:
        xMax = max(sampleLocations.values())
    if zOrder is None:
        zOrder = [i*2+1 for i in range(len(proteins))]
    if medianOnly is None:
        medianOnly = [False]*len(proteins)

    if ax is None:
        f = pylab.figure(figsize=figSize)
        ax = f.add_subplot(111)
    for i, p in enumerate(proteins):
        protDF = df[df[proteinNameHeader] == p]
        ax.plot(numpy.NaN, numpy.NaN, marker = 'o', color=colors[i], label=labelHash[p], markersize=markersize*1.5)
        ysFit = []
        xsFit = []
        for j, s in enumerate(samples):
            if not medianOnly[i]:
                ax.plot([float(sampleLocations[s])]*(len(protDF[protDF[sampleField] == s])),
                         protDF[protDF[sampleField]==s]['currentCalc'], 
                        'o', alpha=alpha, color=colors[i], markersize=markersize, zorder=zOrder[i])
            if median:
                marker='_'
                mew=markersize/3
                zorder=100
                if medianColor is None:
                    medianColor = colors[i]
                if medianOnly[i]:
                    marker='o'
                    medianColor=colors[i]
                    mew = 0
                    zorder=zOrder[i]
                ax.plot(sampleLocations[s], protDF[protDF[sampleField]==s]['currentCalc'].median(), marker=marker, 
                        mew=mew, color=medianColor, markersize=markersize*medianMarkerInc, zorder=zorder)
            ysFit.append(protDF[protDF[sampleField]==s]['currentCalc'].median())
            xsFit.append(sampleLocations[s])
        if fitLine:
            pFit = numpy.polyfit(xsFit, ysFit, 1)
            xsFit.append(xMin)
            xsFit.append(xMax)
            ax.plot(numpy.array(xsFit), pFit[0]*numpy.array(xsFit)+pFit[1], '-', color=colors[i], zorder=zOrder[i]-1)
    if legend:
        ax.legend(ncol=legendCols, loc='upper center', bbox_to_anchor=(0.5, 1.0),
                fancybox=True, shadow=False)
    ax.set_ylim(yMin, yMax)
    ax.set_yticks(yAxis)
    pylab.xticks(sampleLocations.values())
    ax.set_xlim(xMin, xMax)
    if grid:
        pylab.grid()
        cleanAxis(ax, ticksOnly=True)
    if not title is None:
        ax.set_title(title)
    if not yAxisLabel is None:
        ax.set_ylabel(yAxisLabel)
    return ax

def plotMRMCsv(df, samples = None, sampleField='File Name', colors=None, labelHash=None, yAxisLabel=None, legendLoc = 'upper center', areaField = 'Area',
               num = ['light '], den = ['heavy '], proteins = None, alpha=1.0, markersize=10, title=None, legendBBox = (0.5, 1.0), proteinField = 'Protein Name',
               median = True, medianMarker = '-', legend=True, legendCols=5, normProtein=None, medianMarkerInc = 2.5, normSample=None,
               yAxis=None, yMin = 0.0, yMax = 10.0, yTicks=2, figSize=(15,5), grid=False, scaleProt=None, markeredgewidth=1.0):
    df.loc[:,'currentCalc'] = mrmTools.calcValue(df, num, den, field=areaField)
    #print 'running code'
    if scaleProt is None:
        scaleProt = {}
    if samples is None:
        samples = qMS.sort_nicely(list(df[sampleField].unique()))
    if colors is None:
        colors = pylab.cm.jet([float(i)/float(len(samples)) for i in range(len(samples))])
    if proteins is None:
        proteins = qMS.sort_nicely(list(df[proteinField].unique()))
    if yAxis is None:
        interval = (yMax*1.0-yMin)/(yTicks+1.0)
        yAxis = [yMin+i*interval for i in range(0,yTicks+2)]
    if labelHash is None:
        labelHash = {i:i for i in samples}
    for i in proteins:
        if i not in scaleProt.keys():
            scaleProt[i] = 1.0
    if not normSample is None:
        sampDF = df[df[sampleField]==normSample]
        normValue = 1.0
        if not normProtein is None:
            normValue = sampDF[sampDF[proteinField]==normProtein]['currentCalc'].median()        
        for p in proteins:
            scaleProt[p] = (sampDF[sampDF[proteinField]==p]['currentCalc']/normValue).median()
    xOffset = 1.0/(len(samples)+1)
    f = pylab.figure(figsize=figSize)
    ax = f.add_subplot(111)
    for i, s in enumerate(samples):
        sampDF = df[df[sampleField]==s]
        ax.plot(numpy.NaN, numpy.NaN, marker = 'o', color=colors[i], label=labelHash[s], markersize=markersize*1.5, markeredgewidth=markeredgewidth)
        if not normProtein is None:
            normValue = sampDF[sampDF[proteinField]==normProtein]['currentCalc'].median()
        else:
            normValue = 1.0
        for j, p in enumerate(proteins):
            ax.plot([float(j+(i+1)*xOffset)]*len(sampDF[sampDF[proteinField]==p]),
                    sampDF[sampDF[proteinField]==p]['currentCalc']/(normValue*scaleProt[p]), 
                    'o', alpha=alpha, color=colors[i], markersize=markersize, markeredgewidth=markeredgewidth)
            if median:
                try:
                    ax.plot(float(j+(i+1)*xOffset), sampDF[sampDF[proteinField]==p]['currentCalc'].median()/(normValue*scaleProt[p]), marker='_', 
                            mew=markersize/3, color='black', markersize=markersize*medianMarkerInc)
                    #print str(s) + '\t' + str(p) + '\t' + str(sampDF[sampDF[proteinField]==p]['currentCalc'].median()/(normValue*scaleProt[p]))
                except:
                    print 'error'
    if legend:
        ax.legend(ncol=legendCols, loc=legendLoc, bbox_to_anchor=legendBBox,
                fancybox=True, shadow=False)
    ax.set_ylim(yMin, yMax)
    ax.set_yticks(yAxis)
    pylab.xticks(range(len(proteins)), proteins, rotation=45)
    pylab.xlim(0, len(proteins))
    if grid:
        pylab.grid()
        cleanAxis(ax, ticksOnly=True)
    if not title is None:
        ax.set_title(title)
    if not yAxisLabel is None:
        ax.set_ylabel(yAxisLabel)
    return ax

def css_styling(path = "/home/jhdavis/scripts/iPyNBs/custom.css", half=False):
    if half:
        path = "/home/jhdavis/scripts/iPyNBs/custom2.css"
    styles = open(path, "r").read()
    return HTML(styles)

def getBCs(a,n, name=None):
    if a is 'q':
        t='Qualitative'
        if name is None:
            name='Set1'
    elif a is 's':
        t='Sequential'
        if name is None:
            name='Blues'
    elif a is 'd':
        t='Diverging'
        if name is None:
            name='RdBu'
    elif a is 'p':
        t='Qualitative'
        if name is None:
            name='Paired'
    elif a is 'pD':
        t = 'Qualitative'
        if name is None:
            name='Paired'
    elif a is 'pL':
        t = 'Qualitative'
        if name is None:
            name='Paired'    
    else:
        print 'incorrect color type given, expected "q", "s", "p", or "d"'
        return 'ERROR'
    bmap = brewer2mpl.get_map(name, t, n).mpl_colors
    
    if a is 'pL':
        toReturn = [bmap[i] for i in range(0,len(bmap),2)]
    elif a is 'pD':
        toReturn = [bmap[i] for i in range(1,len(bmap),2)]
    else:
        toReturn = bmap
    return toReturn

def smoothListGaussian(data,degree=5):
    window=degree*2
    weight=numpy.array([1.0]*window)
    weightGauss=[]
    dataBuff = []
    for i in range(degree):
        dataBuff.append(data[0])
    for i in data:
        dataBuff.append(i)
    for i in range(degree):
        dataBuff.append(data[-1])
    data = dataBuff
    for i in range(window):
        i=i-degree+1
        frac=i/float(window)
        gauss=1/(numpy.exp((4*(frac))**2))
        weightGauss.append(gauss)
    weight=numpy.array(weightGauss)*weight
    smoothed=[0.0]*(len(data)-window)
    for i in range(len(smoothed)):
        smoothed[i]=sum(numpy.array(data[i:i+window])*weight)/sum(weight)
    
    return smoothed

def rstyle(ax): 
    """Styles an axes to appear like ggplot2
    Must be called after all plot and axis manipulation operations have been carried out (needs to know final tick spacing)
    """
    #set the style of the major and minor grid lines, filled blocks
    ax.grid(True, 'major', color='w', linestyle='-', linewidth=1.4)
    ax.grid(True, 'minor', color='0.92', linestyle='-', linewidth=0.7)
    ax.patch.set_facecolor('0.85')
    ax.set_axisbelow(True)
    
    #set minor tick spacing to 1/2 of the major ticks
    ax.xaxis.set_minor_locator(pylab.MultipleLocator( (pylab.xticks()[0][1]-pylab.xticks()[0][0]) / 2.0 ))
    ax.yaxis.set_minor_locator(pylab.MultipleLocator( (pylab.yticks()[0][1]-pylab.yticks()[0][0]) / 2.0 ))
    
    #remove axis border
    for child in ax.get_children():
        if isinstance(child, matplotlib.spines.Spine):
            child.set_alpha(0)
       
    #restyle the tick lines
    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(5)
        line.set_color("gray")
        line.set_markeredgewidth(1.4)
    
    #remove the minor tick lines    
    for line in ax.xaxis.get_ticklines(minor=True) + ax.yaxis.get_ticklines(minor=True):
        line.set_markersize(0)
    
    #only show bottom left ticks, pointing out of axis
    pylab.rcParams['xtick.direction'] = 'out'
    pylab.rcParams['ytick.direction'] = 'out'
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    
    
    if ax.legend_ <> None:
        lg = ax.legend_
        lg.get_frame().set_linewidth(0)
        lg.get_frame().set_alpha(0.5)
        
        
def setRcs(scale=None, legendScale=None, tickScale=1, axisLineWidth=3):
    """setRcs sets a series of rc params for matplotlib tomake decent looking plots. Also some useful def
    in there to cutomize things.

    :param scale: the default font scale
    :type scale: float
    
    :returns:  nothing, edits the rcparam file.
            
    """ 
    if scale is None:
        scale = 12
    if legendScale is None:
        legendScale = scale*0.85
    
    defaultFont = {'family' : 'sans-serif',
                   'variant' : 'normal',
                   'weight' : 400,
                   'size' : scale*1}

    axisFont = {'titlesize' : scale*1.5,
                'labelsize' : scale*1.25}

    xAxisTicks = {'major.size' : 8.0*tickScale,
                  'minor.size' : 4.0*tickScale,
                  'major.width' : 1.0*tickScale,
                  'minor.width' : 1.0*tickScale,
                  'labelsize' : scale*1,
                  'minor.pad' : 3,
                  'major.pad' : 3,
                  'direction' : 'out'}

    yAxisTicks = {'major.size' : 8.0*tickScale,
                  'minor.size' : 4.0*tickScale,
                  'major.width' : 1.0*tickScale,
                  'minor.width' : 1.0*tickScale,
                  'labelsize' : scale*1,
                  'minor.pad' : 3,
                  'major.pad' : 3,
                  'direction' : 'out'}

    legend = {'fancybox' : True,
              'numpoints' : 1,
              'scatterpoints' : 1,
              'fontsize' : legendScale,
              'borderaxespad' : 1}

    matplotlib.rc('font', **defaultFont)
    matplotlib.rc('axes', **axisFont)
    matplotlib.rc('xtick', **xAxisTicks)
    matplotlib.rc('ytick', **yAxisTicks)
    matplotlib.rc('legend', **legend)

    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('axes', linewidth=axisLineWidth)

def cleanAxis(axisName, ticksOnly=False, complete=False):
    """cleanAxis clears the extra ticks and borders on the axis
    
    :param axisName: a pylab axis to clean
    :type axisName: pylab axis
    :param ticksOnly: boolean if you want to just remove the ticks from the top/right (but leave the border)
    :type ticksOnly: boolean (defaults to False)
    :param complete: boolean if you want to get rid of all fo the borders
    :type complete: boolean (defaults to False)
    
    :returns:  the modified pylab axis.
    """
    axisName.xaxis.set_ticks_position('bottom')
    axisName.yaxis.set_ticks_position('left')
    axisName.xaxis.labelpad = 2
    if not ticksOnly:
        axisName.spines['top'].set_visible(False)
        axisName.spines['right'].set_visible(False)
    if complete:
        axisName.spines['top'].set_visible(False)
        axisName.spines['right'].set_visible(False)
        axisName.spines['bottom'].set_visible(False)
        axisName.spines['left'].set_visible(False)
    return axisName
    
def tickNum(a, xAxis=3, yAxis=3):
    """tickNum sets the number of ticks on a given axis
    
    :param a: a pylab axis to change the tick numbers on
    :type a: pylab axis
    :param xAxis: the number of ticks on the x axis (defaults to 3)
    :type xAxis: int
    :param yAxis: the number of ticks on the y axis (defaults to 3)
    :type yAxis: int
    
    :returns:  the modified pylab axis.
    """
    xlims = [x for x in a.axis()[0:2]]
    xRange = xlims[1] - xlims[0]
    ti = xRange/(xAxis-1)
    a.set_xticks([round(x, 1) for x in numpy.arange(xlims[0], xlims[1]+1, ti)])
    a.set_xlim([round(x, 1) for x in xlims])
    
    ylims = [i for i in a.axis()[2:]]
    yRange = ylims[1] - ylims[0]
    ti = yRange/(yAxis-1)
    a.set_yticks([round(i, 1) for i in numpy.arange(ylims[0], ylims[1]+1, ti)])
    a.set_ylim([round(i, 1) for i in ylims])


def plotStatsDict(statsDict, name='', proteins=None, offset=0.0, markerSize=12, color='#e31a1c', yMax=1.5, alpha=1.0,
                  median=False, figSize = (22,5), noFill=False, mew=1, yMin=-0.05, highlightMed=False, hms=2, hmFilled=True,
                  xTickLabels=None):
    """plotStatsDataStruct plots the contents of a stats dictionary. proteins to be plotted are 
        listed in the non-redundent list, proteins. The data is in statsDict, the name is in in name.
        Decent colors are red (['#ae2221', '#d72c2b', '#e78180']) and blue (['#25557d', '#3170a4', '#5696cc'])

    :param statsDict: a dictionary (easily created by qMS.calcStatsDict)
    :type statsDict: dictionary with keys of proteins names and values of numpy arrays of calcValues
    :param name: the name of the dataset
    :type name: string
    :param proteins: a non-redudent list of the protein names to use (should be the IDs in dataByProtein).
        If none is given, all of the keys will be plotted.
    :type proteins: list
    :param offset: where to center the point (x axis), scales 0 (left edge; default) to 1.0 (right edge)
    :type offset: float
    :param markerSize: size of the dots (default = 12)
    :type markerSize: int
    :param color: color for the dataset (default #e31a1c)
    :type color: color
    :param yMax: the max value for the y axis
    :type yMax: float
    :param median: bool to plot only the median values
    :type median: bool

    :returns:  a pyplot axis with the data plotted
    
        
    """
    if proteins is None:
        proteins = qMS.sort_nicely(statsDict.keys())

    if xTickLabels is None:
        xTickLabels = [item for item in proteins]

    if hmFilled:
        medColor=color
    else:
        medColor='none'
  
    if noFill:
        edgeColor=color
        color='none'
    else:
        edgeColor='black'
        
    xAxis = range(1,len(proteins)+1)
    fig = pylab.figure(figsize=figSize)
    ax = fig.add_subplot(111)

    xs = []
    ys = []
    for x in xAxis:
        p = proteins[x-1]
        if p in statsDict.keys():
            if median:
                xs.append(x+offset)
                ys.append(numpy.median(statsDict[p]))
            else:
                for v in statsDict[p]:
                    xs.append(x+offset)
                    ys.append(v)

    pylab.grid(b=True, which='major', color='grey', linestyle='--', axis='y', linewidth=1.5, alpha=0.5)
    pylab.grid(b=True, which='major', color='grey', linestyle='-', axis='x', linewidth=1.5, alpha=0.75)
    ax.plot(xs, ys, 'o', mfc=color, markeredgecolor=edgeColor, mew=mew, markersize=markerSize, label=name, alpha=alpha)

    if highlightMed:
        mx = []
        my = []
        for x in xAxis:
            p = proteins[x-1]
            if p in statsDict.keys():
                    mx.append(x+offset)
                    my.append(numpy.median(statsDict[p]))
        ax.plot(mx, my, '_', color='black', markeredgecolor='black', mew=2, markersize=markerSize*hms)


    pylab.xticks(xAxis, xTickLabels, rotation=45)
    pylab.xlim(1, len(proteins)+1)
    ####################################
    ####################################
    if yMin == -0.05:
        sub = 0.0
    else:
        sub = yMin
    pylab.yticks([0, (yMax-sub)/5.0, 2*(yMax-sub)/5.0, 3*(yMax-sub)/5.0, 4*(yMax-sub)/5.0, (yMax-sub)])
    pylab.ylim(yMin, yMax)
    ax.set_axisbelow(True)
    return ax

def addStatsDictToPlot(statsDict, ax, name='', offset=0.0, markerSize=12, color='#377db8', median=False, noFill=False, mew=1, 
                       highlightMed=False, hms=2, hmFilled=True, alpha=1.0, proteins=None):
    """addStatsDataStructToPlot adds the contents of a stats dictionary to an existing plot. ONLY PROTEINS
        PRESENT IN THE ORIGINAL PLOT WILL BE PLOTTED IN THE NEW PLOT
        The data is in statsDict, the axis to add to is in ax, the name is in in name.
        Decent colors are red (['#ae2221', '#d72c2b', '#e78180']) and blue (['#25557d', '#3170a4', '#5696cc'])

    :param statsDict: a dictionary (easily created by calcStatsDict)
    :type statsDict: dictionary
    :param ax: the pyplot axis to modify
    :type ax: pyplot axis
    :param name: the name of the dataset
    :type name: string
    :param offset: where to center the point (x axis), scales 0 (left edge; default) to 1.0 (right edge)
    :type offset: float
    :param markerSize: size of the dots (default = 12)
    :type markerSize: int
    :param color: color for the dataset (default #e31a1c)
    :type color: color
    :param median: bool to plot only the median values
    :type median: bool
    :param noFill: a bool if you want open circles
    :type noFill: bool
    :param mew: a float of the marker edge width
    :type mew: float
    
    :returns:  a pyplot axis with the data plotted
        
    """

    if proteins is None:
        a = ax.get_xmajorticklabels()
        proteins = [t.get_text() for t in a]

    if hmFilled:
        medColor=color
    else:
        medColor='none'
    
    if noFill:
        edgeColor=color
        color='none'
    else:
        edgeColor='black'

    xAxis = range(1,len(proteins)+1)

    xs = []
    ys = []
    for x in xAxis:
        p = proteins[x-1]
        if p in statsDict.keys():
            if median:
                xs.append(x+offset)
                ys.append(numpy.median(statsDict[p]))
            else:
                for v in statsDict[p]:
                    xs.append(x+offset)
                    ys.append(v)

    ax.plot(xs, ys, 'o', mfc=color, markeredgecolor=edgeColor, mew=mew, markersize=markerSize, label=name, alpha=alpha)
    
    if highlightMed:
        mx=[]
        my=[]
        for x in xAxis:
            p = proteins[x-1]
            if p in statsDict.keys():
                    mx.append(x+offset)
                    my.append(numpy.median(statsDict[p]))
        ax.plot(mx, my, '_', color='black', markeredgecolor='black', mew=2, markersize=markerSize*hms)

    ax.set_axisbelow(True)
    return ax

def makePlotWithFileList(isoFileList, numerator, denominator, subunits=None, normProtein=None, yMax=1.5, title=None, legendCols=4,
                         median=False, names=None, colors=None, figSize=(22,5), markerSize=None, noFill=False, legend=False,
                        mew=1, yMin=-0.05, highlightMed=False, hms=2, hmFilled=True, yAxis=None, alpha=1.0, adjProt=None):
    """makePlotWithFileList is a  helper function that plots massage-style data from a list of files

    :param isoFileList: a list of the files to be ploted (shoudl be full path to _iso.csv files)
    :type isoFileList: list of strings
    :param numerator: strings of the keys in the numerator (ampu, ampl, amps)
    :type numerator: list of strings
    :param denominator: strings of the keys in the denominator (ampu, ampl, amps)
    :type denominator: list of strings
    :param subunits: a list of the proteins to be plotted - strings must match the keys of the 'proteins' field in _iso.csv
    :type subunits: list
    :param normProtein: string of the protein to normalize to for all datasets (defaults to None)
    :type normProtein: string
    :param yMax: float of maximum y value
    :type yMax: float
    :param median: bool to plot only the median values
    :type median: bool
    :param names: a list of the names to be listed in the legend; must be same length as isoFileList
    :type names: list of strings
    :param colors: a list of the colors to be used in plotting, again must be same length as isoFileList
    :type colors: list of strings
    :param figSize: the figure size (list - (10,10))
    :type figSize: list of floats
    :param noFill: a bool if you want open circles
    :type noFill: bool
    :param mew: a float of the marker edge width
    :type mew: float

    :returns: the plotted axis
    
    """
    if names is None:
        names = isoFileList
        
    namesList = [(isoFileList[i], names[i]) for i in range(len(isoFileList))]
    allStats = qMS.multiStatsDict(isoFileList, numerator, denominator, normalization=1.0, offset=0.0, normProtein=normProtein, noProcess=True, adjProt=adjProt)
    
    return makePlotWithStatsDictDict(allStats, subunits=subunits, yMax=yMax, title=title, legend=legend, legendCols=legendCols, yAxis=yAxis,
                                     median=median, namesList=namesList, colors=colors, figSize=figSize, markerSize=markerSize,
                                     noFill=noFill, mew=mew, yMin=yMin, highlightMed=highlightMed, hms=hms, hmFilled=hmFilled, alpha=alpha)

def makePlotWithStatsDictDict(allStats, subunits=None, yMax=1.5, title=None, legend=False, legendCols=4, yAxis=None,
                         median=False, namesList=None, colors=None, figSize=(22,5), markerSize=None, 
                        noFill=False, mew=1, yMin=-0.05, highlightMed=False, hms=2, hmFilled=True, alpha=1.0,
                        xTickLabels=None):
    """makePlotWithStatsDictDict is a  helper function that plots massage-style data from a dict of statsDicts

    :param allStats: dictionary of stats dictionaries (returned by multiStatsDict)
    :type allStats: dict of statsDicts
    :param subunits: a list of the proteins to be plotted - strings must match the keys of the 'proteins' field in _iso.csv
    :type subunits: list
    :param yMax: float of maximum y value
    :type yMax: float
    :param median: bool to plot only the median values
    :type median: bool
    :param namesList: a list of lists where each sublist is a pair with the key for the 
        allStats in [0] and the name for the legend in [1]. typically, the [0] would be a full
        path.
    :type namesList: list of pairs of strings
    :param colors: a list of the colors to be used in plotting, again must be same length as allStats
    :type colors: list of strings
    :param figSize: the figure size (list - (10,10))
    :type figSize: list of floats
    :param markerSize: the size of the dots
    :type markerSize: float
    :param noFill: a bool if you want open circles
    :type noFill: bool
    :param mew: a float of the marker edge width
    :type mew: float

    :returns: the plotted axis
    
    """

    if namesList is None:
        namesList = [(k, k) for k in allStats.keys()]
    
    if colors is None:
        colors = [pylab.cm.jet(float(i)/float(len(allStats))) for i in range(len(namesList))]
    
    offsets = float(len(namesList)+1)
    if markerSize is None:
        markerSize = (20.0/offsets)+4

    ax = plotStatsDict( allStats[namesList[0][0]], name=namesList[0][1], proteins=subunits, \
                        offset=1.0/offsets, markerSize=markerSize, yMax=yMax, median=median, 
                        color=colors[0], figSize=figSize, noFill=noFill, mew=mew, yMin=yMin, 
                        highlightMed=highlightMed, hms=hms, hmFilled=hmFilled, alpha=alpha,
                        xTickLabels=xTickLabels)
        
    for i in range(1,len(namesList)):
        ax = addStatsDictToPlot(allStats[namesList[i][0]], ax, name=namesList[i][1], \
                                offset=(1.0/offsets)*(i+1.25), markerSize=markerSize, median=median, 
                                color=colors[i], noFill=noFill, mew=mew, highlightMed=highlightMed, 
                                hms=hms, hmFilled=hmFilled, alpha=alpha,
                                proteins=subunits)
        ax.set_axisbelow(True)
    if legend:
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0),
          ncol=legendCols, fancybox=True, shadow=False)
    if not title is None:
        ax.set_title(title)
    if not yAxis is None:
        ax.set_ylabel(yAxis)
    return ax

def plotCompData(xdat, ydat, proteins, title=None, xlabel='dat1', ylabel='dat2', xMax=1.5, yMax=1.5, figSize=(10,10), saveFile=None):
    """plotCompData is a  makes a scatter plot out of data to be compared. Useful for the pool size stuff. Only the median values are plotted.

    :param xdat: a dictionary of values to be plotted on the x axis (keys are protein names, values are numpy arrays of values)
    :type xdat: a dictionary of numpy arrays (keys=protein names)
    :param ydat: a dictionary of values to be plotted on the y axis (keys are protein names, values are numpy arrays of values)
    :type ydat: a dictionary of numpy arrays (keys=protein names)
    :param proteins: a list of the proteins to be plotted - strings must match the keys of the 'proteins' field in xdat, ydat
    :type proteins: list
    :param title: the plot title
    :type title: string
    :param xlabel: the plot xlabel
    :type xlabel: string
    :param ylabel: the plot ylabel
    :type ylabel: string
    :param xMax: float of maximum x value
    :type xMax: float
    :param yMax: float of maximum y value
    :type yMax: float
    :param figSize: the figure size (list - (10,10))
    :type figSize: list of floats
    :param saveFile: the full path of the file to be saved to
    :type saveFile: string

    :returns: the plotted axis
    
    """
    x = [numpy.median(xdat[i]) for i in proteins]
    y = [numpy.median(ydat[i]) for i in proteins]    
    scat = pylab.figure(figsize=figSize)
    scatAx = scat.add_subplot(111)    
    scatAx.scatter(x,y, c='b', s=150)
    scatAx.set_title(title)
    scatAx.set_xlabel(xlabel)
    scatAx.set_ylabel(ylabel)
    scatAx.set_xlim([-0.1,xMax])
    scatAx.set_ylim([-0.1,yMax])
    scatAx.set_xticks([0,xMax/5,xMax/5*2,xMax/5*3,xMax/5*4,xMax])
    scatAx.set_yticks([0,yMax/5,yMax/5*2,yMax/5*3,yMax/5*4,yMax])
    scatAx.yaxis.tick_left()
    scatAx.xaxis.tick_bottom()
    for prot, xl, yl in zip(proteins, x, y):
        scatAx.annotate(str(prot[4:]), xy = (float(xl), float(yl)), xytext = (15,15), textcoords = 'offset points', arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    pylab.tight_layout()
    scatAx.plot(numpy.linspace(0, 10), numpy.linspace(0,10))
    return scatAx

def proteinScatterPlot(yDataDict, xData, xMin=0, xMax=None, yMin=-0.1, yMax=10,
                       title=None, xLabel=None, yLabel=None, colors=None, 
                       figSize=(10,10), markerSize=10, legend=False, alpha=1.0, marker='o',
                        linestyle=None, xTicks=None, legendLoc='upper left', legendCols=2, axes=None):
    """proteinScatterPlot is a  makes a scatter plot out of ydata with a fixed x. Useful for the standard curve stuff or gradients.

    :param yDataDict: a dictionary of values to be plotted on the y axis (keys are protein names, values are numpy arrays of values)
    :type yDataDict: a dictionary of numpy arrays (keys=protein names)
    :param xData: a list of the xvalues (same for all proteins)
    :type xData: a list of floats
    :param xMin: float of minimum x value
    :type xMin: float
    :param yMin: float of minimum y value
    :type yMin: float
    :param xMax: float of maximum x value
    :type xMax: float
    :param yMax: float of maximum y value
    :type yMax: float
    :param title: the plot title
    :type title: string
    :param xlabel: the plot xlabel
    :type xlabel: string
    :param ylabel: the plot ylabel
    :type ylabel: string
    :param colors: a list of colors to be used
    :type colors: list of colors
    :param figSize: the figure size (list - (10,10))
    :type figSize: list of floats
    :param markerSize: the size of the dots
    :type markerSize: float
    :param legend: a bool of whether to have a legend
    :type legend: bool
    :param linestyle: the connector line style
    :type linestyle: the connector line style
    :param xTicks: a list of the xTicks
    :type xTicks: a list of floats
    :param legendloc: a string of where to locate the legend
    :type legendloc: string
    :param legendCols: the number of columns in the legend
    :type legendCols: int

    :returns: the plotted axis
    
    """    
    if xMax is None:
        xMax = max(xData)
    if colors is None:
        colors = [pylab.cm.jet(float(i)/float(len(yDataDict))) for i in range(len(yDataDict))]
    if axes is None:
        scat = pylab.figure(figsize=figSize)
        scatAx = scat.add_subplot(111)
    else:
        scatAx=axes
    for i,p in enumerate(qMS.sort_nicely(yDataDict.keys())):
        
        if not (linestyle is None):
            scatAx.plot(xData, yDataDict[p], c=colors[i], linestyle=linestyle, label=p, marker=marker, markersize=markerSize, alpha=alpha)
        else:
            scatAx.plot(xData, yDataDict[p], c=colors[i], markersize=markerSize, marker=marker, label=p, alpha=alpha)
    scatAx.set_title(title, multialignment='center')
    scatAx.set_xlabel(xLabel)
    scatAx.set_ylabel(yLabel)
    scatAx.set_xlim([xMin,xMax])
    scatAx.set_ylim([yMin,yMax])
    if xTicks is None:
        scatAx.set_xticks([0,xMax/4,xMax/4*2,xMax/4*3,xMax])
    else:
        scatAx.set_xticks(xTicks)
    scatAx.set_yticks([0,yMax/4,yMax/4*2,yMax/4*3,yMax])
    if legend:
        pylab.legend(loc=legendLoc, ncol=legendCols, scatterpoints=1)
    scatAx.yaxis.tick_left()
    scatAx.xaxis.tick_bottom()
    pylab.tight_layout()
    
    return scatAx

    

def plotMSSpectra3D(listOfFilesToPlot, listOfNames=None, listOfColors=None, gridLines=False, yMin=0.5, yMax=2.5, yScale = 1.0,
                    legend=True, normalizeToN15=False, subtractRef=None, legendLoc=4, lw=1.5, xMin=0, xMax=2000, scaleP=False, scaleI=0, scaleVal=1.0,
                    figsize=(10,10), tLeft=0, tRight=-1, fixedOffset=False, noTicks=False, xlabel='mass', zlabel='intensity', a14=1.0):
    """plotMSSpectra3D is a  makes a 3d plot of MS spectra

    :param listOfFilesToPlot: a list of the spectra to be plotted (full paths)
    :type listOfFilesToPlot: list of strings
    :param listOfNames: a list of the names for each dataset
    :type listOfNames: list of strings
    :param listOfColors: a list of colors to be used
    :type listOfColors: list of colors
    :param gridLines: a bool of whether to draw gridlines
    :type gridLines: bool
    :param yMin: float of minimum y value
    :type yMin: float
    :param yMax: float of maximum y value
    :type yMax: float
    :param legend: a bool of whether to have a legend
    :type legend: bool
    :param normalizeToN15: a bool of whether to nomrmalize each plot the N15 maximum
    :type normalizeToN15: bool
    :param subtractRef: a int pointing to which is the refernece spectra that should be subtracted from each 
    :type subtractRef: int

    :returns: the plotted axis
    
    """  
    if listOfNames==None:
        listOfNames = listOfFilesToPlot
    if listOfColors==None:
        listOfColors = [pylab.cm.jet(float(i)/float(len(listOfFilesToPlot))) for i in range(len(listOfFilesToPlot))]
    
    fig = pylab.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')

    yTotal = len(listOfFilesToPlot)
    top = 0.0

    if not (subtractRef is None):
        [bhah, zsRef, blah] = qMS.readMSSpectraFile(listOfFilesToPlot[subtractRef])
        #zsRef = list(numpy.array(zsRef)-fixedOffset
        zNorm = max(zsRef[len(zsRef)/2:])
        zsRef = numpy.array(zsRef)/zNorm
    
    for i,f in enumerate(listOfFilesToPlot):
        [xs, zs, name] = qMS.readMSSpectraFile(f)
        if fixedOffset:
            off = zs[len(zs)/2]
            print off
            #zs = list(numpy.array(zs)-zs[len(zs)/2])
            zs = list(numpy.array(zs)-off)
        ys = [yTotal-i]*len(xs)
        ys = numpy.array(ys)*yScale
        if normalizeToN15:
            zNorm = max(zs[len(zs)/2:])
            zs = numpy.array(zs)/zNorm
        if not (subtractRef is None):
            zNorm = max(zs[len(zs)/2:])
            zs = numpy.array(zs)/zNorm
            zs[:len(zsRef)/2] = zs[:len(zs)/2]-zsRef[:len(zsRef)/2]
            zs = zs*zNorm
            #xs = xs[:len(xs)/2]
            #ys = ys[:len(ys)/2]
            #zs = zs[:len(zs)/2]
        zs[:len(zs)/2] = numpy.array(zs[:len(zs)/2])*a14
        if (scaleP is True) and (i==scaleI):
            zs = numpy.array(zs)*scaleVal
        ax.plot(numpy.array(xs[tLeft:tRight]),numpy.array(ys[tLeft:tRight]),numpy.array(zs[tLeft:tRight]), color=listOfColors[i], lw=lw, label=listOfNames[i])
        top = max([top, float(max(zs))])


    ax.w_xaxis.pane.set_visible(False)
    ax.w_yaxis.pane.set_visible(False)
    ax.w_zaxis.pane.set_visible(False)

    if gridLines:        
        ax.w_xaxis.gridlines.set_linewidth(1)
        ax.w_yaxis.gridlines.set_linewidth(1)
        ax.w_zaxis.gridlines.set_linewidth(1)
    
    else:
        ax.w_xaxis.gridlines.set_visible(False)
        ax.w_yaxis.gridlines.set_visible(False)
        ax.w_zaxis.gridlines.set_visible(False)

    [i.set_linewidth(1) for i in ax.w_xaxis.get_ticklines()]
    [i.set_linewidth(1) for i in ax.w_yaxis.get_ticklines()]
    [i.set_linewidth(1) for i in ax.w_zaxis.get_ticklines()]

    ax.w_xaxis.line.set_linewidth(1)
    ax.w_yaxis.line.set_linewidth(1)
    ax.w_zaxis.line.set_linewidth(1)
               
    ax.set_zticks([round(i,1) for i in [0, top/3, 2*top/3, top]])
    ax.set_zlim3d([0, top])
    ax.set_ylim3d(yMin, yMax)
    ax.set_yticks(range(1,yTotal+1))
    pylab.yticks(range(1,yTotal+1), ['']*yTotal)
    ax.set_xlim3d([xMin, xMax])
    
    if noTicks:
        ax.set_zticks([])
        ax.set_xticks([])
        ax.set_yticks([])

    ax.set_xlabel(xlabel)
    ax.set_zlabel(zlabel)

    ax.view_init(15, -60)
    if legend:
        pylab.legend(loc=legendLoc)
    
    pylab.tight_layout()
    return ax

def drawHeatMap(xdat, name=None, colors=pylab.cm.Reds, dendro=False, protColors=None, cIndex=None, km=None, 
                nameDict={}, scale=None, saveName=None, colorBar=False, figSize=(6,6), topDendro=False, fig=None, axData=None):
    """drawHeatMap produces a colored heatmap in a new figure window

    :param xdat: a data object (must contain fields 'data', 'fractions', 'proteins')
    :type xdat: dict
    :param colors: a color scale (a cmap)
    :type colors: cmap
    :param name: figure name and title
    :type name: str.
    :param dendro: a boolean to draw the dendrogram on the left
    :type dendro: bool.
    :param protColors: a color map used to label the protein names with group colors
    :type protColors: cmap
    :param cIndex: a list of groupIds for the proteins
    :type cIndex: list
    :param km: if present, will draw the kmeans cluster profiles at the top of the figure- input is a 2d-matrix - rowVectors for each centroid, each column is a fraction
    :type km: matrix
    :returns:  int -- the return code.
    :raises: AttributeError, KeyError
    :returns: a figure object

    """

    data = xdat['data']
    if nameDict is None:
        nameList = [i for i in xdat['fractions']]
    else:
        nameList = [nameDict[i] for i in xdat['fractions']]

    proteins = [i for i in xdat['proteins']]
    if fig is None:
        fig = pylab.figure(figsize=figSize)
    if not (name is None):
        fig.suptitle(name)
    ##Draw heatmap
    xOffset = 0.05
    if colorBar:
        xLess = 0.10
    else:
        xLess = 0.00
    if dendro:
        xStart = 0.375
        xLength = 0.55-xLess
    else:
        xStart = 0.125
        xLength = 0.85-xLess
    if (km is None) and (topDendro is False):
        yStart = 0.05
        yLength = 0.9
    else:
        yStart = 0.05
        yLength = 0.8
    figAxes = heatMapAxes(data, dims = [xStart, yStart, xLength, yLength], columns=nameList, rows=proteins, protColors=protColors, cIndex=cIndex, fig=fig, colors=colors, axData=axData)
    ##Draw colorbar
    if colorBar:
        fig.colorbar(figAxes)

    if dendro:
        ax2Data = fig.add_axes([xOffset, yStart, xLength-0.3, yLength])
        sch.dendrogram(xdat['rightDendro'], orientation='right', color_threshold=0.0)
        ax2Data.set_xticks([])
        ax2Data.set_yticks([])
    
    if topDendro:
        ax4Data = fig.add_axes([xStart, yStart+yLength, xLength, 0.1])
        sch.dendrogram(xdat['topDendro'], orientation='down', color_threshold=0.0)
        ax4Data.set_xticks([])
        ax4Data.set_yticks([])
    
    if not km is None:
        small = data.min()
        big = data.max()
        if math.fabs(small) > math.fabs(big):
            big = 0-small
        else:
            small = 0-big
        offset=0.0
        ax3Data = fig.add_axes([xStart, yLength+offset, xLength-0.1, 0.1])
        ax3Data.matshow(km, aspect='auto', origin='lower', cmap=colors, vmin=small, vmax=big)
        for i in range(len(km)):
            ax3Data.text(-0.75, i, 'clus'+str(i), verticalalignment="center", horizontalalignment="right", fontsize=10, color=cIndex(float(i)/(protColors.max()+1)))
        ax3Data.set_xticks([])
        ax3Data.set_yticks([])
    #fig.tight_layout()
    if not (saveName is None):
        pylab.savefig(saveName)
        
    return fig

def heatMapAxes(data, dims=[0.1, 0.1, 0.7, 0.7], colors=pylab.cm.autumn, columns=None, rows=None, protColors=None, cIndex=None, fig=None, colorBar=False, axData=None):
    """heatMapAxes draws a heatmap

    :param data: a datamatrix to draw
    :type xdat: a 2D Matrix
    :param dims: the size of the plot to draw - defaults to full window
    :type dims: list (4 elements)
    :param colors: color index to use - defaults to redblue
    :type colors: cmap
    :param fractions: fraction names
    :type fractions: list
    :param proteins: protein names
    :type proteins: list
    :param protColors: a color map used to label the protein names with group colors
    :type protColors: cmap
    :param cIndex: a list of groupIds for the proteins
    :type cIndex: list
    :param fig: where to plot the axes (which figure); defaults to new figure
    :type fig: matplotlib figure
    :returns:  an axes

    """
    if fig is None:
        fig = pylab.figure()
    if axData is None:
        axData = fig.add_axes(dims)
    for i in range(len(columns)):
        axData.text(i, -0.5 , ' '+str(columns[i]), rotation=270, verticalalignment="top", horizontalalignment="center", fontsize=12)
    if protColors == None:
        for i in range(len(rows)):
            axData.text(-0.75, i, '  '+str(rows[i]), verticalalignment="center", horizontalalignment="right", fontsize=12)
    else:
        for i in range(len(rows)):
            axData.text(-0.75, i, '  '+str(rows[i]), verticalalignment="center", horizontalalignment="right", fontsize=12, color=cIndex(float(protColors[i])/(protColors.max()+1)))
    small = data.min()
    big = data.max()
    if math.fabs(small) > math.fabs(big):
        big = 0-small
    else:
        small = 0-big
    masked_array = numpy.ma.array (data, mask=numpy.isnan(data))
    colors.set_bad('grey',1.)
    figData = axData.imshow(masked_array, interpolation='nearest', cmap=colors, aspect='auto', origin='lower')
    if colorBar:
        fig.colorbar(figData, ax=axData, ticks=[0, 0.25, 0.50, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0], pad=0.01, extend='neither')
    axData.set_xticks([])
    axData.set_yticks([])
    return figData

def findIndices(g):
    """findIndices is a  helper function that likely should be deleted
    
    """
    change = [0]
    seen = [g[0]]
    for i in range(1, len(g)):
        if not g[i] in seen:
            change.append(i)
            seen.append(g[i])
    return change
    
def mapGroups(groupList, letters):
    """mapGroups is a  helper function that maps groups numbers to letters - likely should be deleted

    :param groupList: a list to be mapped
    :type d: list
    :param letters: a list to be mapped onto
    :type letters: list
    :returns: a list with elements of groupList mapped onto the letters
    
    """
    changeList = findIndices(groupList)
    i = 0
    for index in changeList:
        toReplace = groupList[index]
        groupList = qMS.listReplace(groupList, toReplace, letters[i])
        i = i+1
    return list(groupList)

def plotPeptideCount(files, labels=None, colors=None, proteins=None, num=['AMP_U'], den=['AMP_U'], legend=True, figSize=(20,10)):
    
    if labels is None:
        labels=files
    
    if colors is None:
        colors = [pylab.cm.jet(float(i)/float(len(files))) for i in range(len(files))]
    
    numFiles=len(files)
    offset = 1.0/(numFiles+1)
    
    msd = qMS.multiStatsDict(files, num, den)
    
    f = pylab.figure(figsize=figSize)
    a = f.add_subplot(111)
    for j in range(0,numFiles):
        a.bar([i+offset*j for i in range(len(proteins))], [len(msd[files[j]][i]) for i in proteins], align='edge', width=1.0/(numFiles+2), color=colors[j], label=labels[j])
    a.set_xticks([i+0.4 for i in range(0,len(proteins))])
    a.set_xticklabels(proteins, size='medium', rotation='90')
    a.set_ylabel('number of peptides')
    cleanAxis(a, ticksOnly=True)
    if legend:
        a.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=5, mode="expand", borderaxespad=0.)
    return a

def mm2inches(mm):
    return mm*0.0393701
