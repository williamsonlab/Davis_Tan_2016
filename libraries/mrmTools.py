# -*- coding: utf-8 -*-
"""
Created on Mon May  5 18:23:21 2014

@author: jhdavis
"""

import numpy
import qMS
import vizLib
import matplotlib.pyplot as plt
import pylab
from matplotlib import gridspec
import pandas
import qMSDefs

def calcErrorMRM(dataFrame):
    dataFrame['error'] = abs((dataFrame['light Area']/dataFrame['heavy Area']) - (dataFrame['light TotalArea'])/(dataFrame['heavy TotalArea']))
    dataFrame['errorFrac'] = dataFrame['error']/(dataFrame['light TotalArea']/dataFrame['heavy TotalArea'])
    dataFrame['RTOffset'] = dataFrame['light RetentionTime'] - dataFrame['heavy RetentionTime']
    dataFrame['light FracTotalArea'] = dataFrame['light Area']/(dataFrame['light TotalArea'])
    dataFrame['heavy FracTotalArea'] = dataFrame['heavy Area']/(dataFrame['heavy TotalArea'])
    return dataFrame

def scoreDatasetsMRM(df, lppml, lppmh, hppml, hppmh,
                   la, ha,
                   dp, fal, fah,
                    e, ef, rtoL, rtoH):
    pandas.options.mode.chained_assignment = None
    df['score'] = 0
    df['score'] = df['score'] + ((df['light MassErrorPPM'] > lppml) & (df['light MassErrorPPM'] < lppmh))*1 + \
                                ((df['heavy MassErrorPPM'] > hppml) & (df['heavy MassErrorPPM'] < hppmh))*1 + \
                                (df['light Area'] > la)*1 + \
                                (df['heavy Area'] > ha)*1 + \
                                (df['DotProductLightToHeavy'] > dp)*1 + \
                                (df['light FracTotalArea'] > fal)*1 + \
                                (df['heavy FracTotalArea'] > fah)*1 + \
                                (df['error'] < e)*1 + \
                                (df['errorFrac'] < ef)*1 + \
                                ((df['RTOffset'] > rtoL) & (df['RTOffset'] < rtoH))*1
    return df    

def correctOccupancyMRM(ref, toCorr, maxZero=True):
    corr = {}
    for p in ref.keys():
        med = numpy.median(ref[p])
        if maxZero:
            corr[p] = numpy.array([max(i-med, 0.0) for i in toCorr[p]])
        else:
            corr[p] = numpy.array([i-med for i in toCorr[p]])
    return corr

def correctStatsDictDictMRM(sdd, ref):
    corrSdd = {}
    for k in sdd.keys():
        corrSdd[k] = correctOccupancyMRM(ref, sdd[k])
    return corrSdd

def normStatsDictMRM(sd, normValue=1.0, normProtein=None):
    normSD = {}
    if not normProtein==None:
        normValue = 1.0/numpy.median(sd[normProtein])
    for p in sd:
        normSD[p] = numpy.array([i*normValue for i in sd[p]])
    return normSD

def normStatsDictDictMRM(sdd, normValue=1.0, normProtein=None):
    normSDD = {}
    for k in sdd.keys():
        normSDD[k] = normStatsDictMRM(sdd[k], normValue, normProtein)
    return normSDD

def makeStatsDictMRM(dataFrame, fileName, field = 'currentCalc', proteinList = None, filterField = 'allClear'):
    if proteinList is None:
        proteinList = list(dataFrame['Protein Name'].unique())
    pDict = {}
    fileFrame = dataFrame[dataFrame['File Name'] == fileName]
    for p in proteinList:
        pDict[p] = fileFrame[(fileFrame['Protein Name'] == p) & (fileFrame[filterField])][field].values
    return pDict

def makeFileStatsDictMRM(dataFrame, listOfFiles=None, field = 'currentCalc', proteinList = None, filterField = 'allClear'):
    if listOfFiles is None:
        listOfFiles = list(dataFrame['File Name'].unique())
    filePDict = {}
    for i in listOfFiles:
        filePDict[i] = makeStatsDictMRM(dataFrame, i, field=field, proteinList = proteinList, filterField = filterField)
    return filePDict
    
def getAllOccupancyMRM(dataFrame, fileName, listOfProteins=None,
                       num=['light'], den=['heavy'], total=False,
                        normProtein=None, normValue=1.0, offset=0.0, selectField='File Name'):
    if listOfProteins is None:
        listOfProteins = list(dataFrame['Protein Name'].unique())
    pDict = {}
    for p in listOfProteins:
        pDict[p] = getOccupancyMRM(dataFrame, fileName, p, num=num, den=den, 
                                    total=total, normValue=normValue, 
                                    offset=offset, allData=False, selectField=selectField).values
    if not (normProtein is None):
        normValue = 1/numpy.median(pDict[normProtein])
        for p in listOfProteins:
            pDict[p] = getOccupancyMRM(dataFrame, fileName, p, num=num, den=den, 
                                        total=total, normValue=normValue, 
                                        offset=offset, allData=False, selectField=selectField).values
    return pDict
   
def getAllOccupancyFileListMRM(dataFrame, fileList, listOfProteins=None, 
                               num=['light'], den=['heavy'], total=False, 
                                normProtein=None, normValue=1.0, offset=0.0, selectField='File Name'):
    filePDict = {}
    for i in fileList:
        filePDict[i] = getAllOccupancyMRM(dataFrame, i, listOfProteins=listOfProteins, 
                                            num=num, den=den, total=total, normProtein=normProtein, 
                                            normValue=normValue, offset=offset, selectField=selectField)
    return filePDict

def getOccupancyMRM(dataFrame, fileName, proteinName, num=['light'], den=['heavy'], normValue=1.0, offset=0.0, total=False, allData=False, selectField='File Name'):
    """getOccupancy takes a pandas dataframe generated from reading an MRM CSV file
        as well as a fileName (eg.wiff) and a a proteinName. It calculates a 
        protein occupancy using the fields specified in the numerator and 
        denominator lists.

    :param dataFrame: the dataFrame to work on. Should be generated from a 
        Skyline MRM export results CSV (e.g. dataFrame = pandas.read_csv(csvFile, na_values='#N/A'))
        Must bear the columns "FileName, ProteinName, Area, TotalArea
    :type dataFrame: a pandas dataframe
    :param fileName: the file to consider (the wiff file)
    :type fileName: a string
    :param proteinName: the protein to consider
    :type proteinName: a string
    :param proteinName: the protein to consider
    :type proteinName: a string
    :param num: species to calc in the numerator
    :type num: a list of strings (light or heavy)
    :param den: species to calc in the denominator
    :type den: a list of strings (light or heavy)
    :param total: boolean to calculate the total or on a product by product basis
    :type total: boolean (defaults to false)
    :param allData: a boolean to return the full dataframe or just the calculated value
    :type allData: boolean (defaults to false)
    
    :returns:  a pandas dataframe with the calculated value (either appended or on its own)

    """
    allProducts = dataFrame[(dataFrame[selectField]==fileName) & (dataFrame['Protein Name']==proteinName)]
    
    if total:
        stringAppend = ' Total Area'
    else:
        stringAppend = ' Area'
    
    allProducts['num'] = 0.0
    allProducts['den'] = 0.0
    
    for i in num:
        allProducts['num'] = allProducts['num'] + allProducts[i + stringAppend]
    for i in den:
        allProducts['den'] = allProducts['den'] + allProducts[i + stringAppend]
    
    allProducts['calcValue'] = (allProducts['num']/allProducts['den'])*normValue+offset

    if allData:
        return allProducts
    else:
        if total:
            return qMS.dropDuplicatesPandas(allProducts)['calcValue'].dropna()
        else:
            return allProducts['calcValue'].dropna()
            
def getInfoMRM(dataFrame, index):
    """getInfoMRM takes a pandas dataframe (read from the tsv output from skyline),
        and an index into that dataframe. Returns info related to that index 
        in the order [filename, mod. pep. seq, precharge, fragIon, prodCharge, isotope]

    :param dataFrame: pandas dataFrame - generated using tsv from skyline and the command
        dataFrame = pandas.read_csv(file, na_values='#N/A', sep='\t')
    :type dataFrame: pandas dataframe
    :param index: a int of what index to inspect
    :type index: int
    :returns:  a list of data, [filename, mod. pep. seq, precharge, fragIon, prodCharge, isotope]

    """
    fn = dataFrame.ix[index]['FileName']
    pepSeq = dataFrame.ix[index]['Peptide Modified Sequence']
    precursorCharge = dataFrame.ix[index]['Precursor Charge']
    fragIon = dataFrame.ix[index]['Fragment Ion']
    prodCharge = dataFrame.ix[index]['Product Charge']
    isotope = dataFrame.ix[index]['IsotopeLabel Type']
    return [fn, pepSeq, precursorCharge, fragIon, prodCharge, isotope]

def getTotalChromatographMRM(dataFrame, index, fragIon=None):
    [fn, pms, preC, fi, proC, isotope] = getInfoMRM(dataFrame, index)
    subFrame = dataFrame[(dataFrame['File Name'] == fn) &
                         (dataFrame['Peptide Modified Sequence'] == pms) &
                         (dataFrame['Precursor Charge'] == preC) &
                         (dataFrame['Isotope Label Type'] == isotope)]
    if fragIon is None:
        first = subFrame['Intensities'].values[0].split(',')
        totArray = numpy.array([float(i) for i in first])
        for a in subFrame['Intensities'].values[1:]:
            b = numpy.array([float(x) for x in a.split(',')])
            totArray = numpy.add(totArray,b)
    else:
        print "NEED TO WRITE THIS CODE TO DEAL WITH SUBSETS OF FRAGMENT IONS"
    
    return [numpy.array([float(x) for x in subFrame['Times'].values[0].split(',')]), totArray]

def getPairedIonMRM(dataFrame, index):
    [fn, pms, preC, fi, proC, isotope] = getInfoMRM(dataFrame, index)
    if isotope == 'light':
        isotope = 'heavy'
    else:
        isotope = 'light'
    return dataFrame[(dataFrame['File Name'] == fn) &
                         (dataFrame['Peptide Modified Sequence'] == pms) &
                         (dataFrame['Precursor Charge'] == preC) &
                         (dataFrame['Fragment Ion'] == fi) &
                         (dataFrame['Product Charge'] == proC) &
                         (dataFrame['Isotope Label Type'] == isotope)].index.values[0]

def getIsotopPairTotalsMRM(dataFrame, index):
    [fn, pms, preC, fi, proC, isotope] = getInfoMRM(dataFrame, index)
    indexOpp = getPairedIonMRM(dataFrame, index)
    if isotope == 'light':
        [lX, lY] = getTotalChromatographMRM(dataFrame, index)
        [hX, hY] = getTotalChromatographMRM(dataFrame, indexOpp)
    
    else:
        [lX, lY] = getTotalChromatographMRM(dataFrame, indexOpp)
        [hX, hY] = getTotalChromatographMRM(dataFrame, index)

    return [lX, lY, hX, hY]

def getRelatedIndeciesMRM(dataFrame, index):
    [fn, ps, preC, fi, proC, i] = getInfoMRM(dataFrame, index)
    return dataFrame[(dataFrame['Isotope Label Type'] == i) &
                        (dataFrame['File Name'] == fn) &
                        (dataFrame['Peptide Modified Sequence'] == ps) &
                        (dataFrame['Precursor Charge'] == preC)].index.values


def plotTotalChromPairsMRM(dataFrame, toPlotIndex, axis, 
                            colors=['blue', 'red'], smooth=0, zoom=False):
    [lightX, lightY, heavyX, heavyY] = getIsotopPairTotalsMRM(dataFrame, toPlotIndex)
    if smooth > 0:
        lightY = vizLib.smoothListGaussian([float(i) for i in lightY], degree=smooth)
        heavyY = vizLib.smoothListGaussian([float(i) for i in heavyY], degree=smooth)
    axis.plot(lightX, lightY, colors[0], label='light')
    axis.plot(heavyX, heavyY, colors[1], label='heavy')
    if zoom:
        maxIndex = numpy.argmax(numpy.array(lightY))
        axis.set_xlim(lightX[maxIndex]-2, lightX[maxIndex]+2)
    return axis

                        
def plotAllTransitionsMRM(dataFrame, index, a, colors=None, smooth=0, zoom=False):
    allIndicies = getRelatedIndeciesMRM(dataFrame, index)
    colors = vizLib.getBCs('q', min(len(allIndicies), 9))
    colors.append('grey')
    for i in range(len(allIndicies)):
        [fn, ps, preC, fi, proC, iso] = getInfoMRM(dataFrame, allIndicies[i])
        a = plotMRM_index(dataFrame, allIndicies[i], a, color=colors[i], 
                          smooth=smooth, zoom=zoom)
    return a

def plotMRM(dataFrame, fileName, pepSeq, precursorCharge, fragIon, prodCharge, 
            isotopeLabel, ax, color='grey', smooth=0, zoom=False):
    subDF = dataFrame[(dataFrame['Peptide Modified Sequence'] == pepSeq) & 
                        (dataFrame['Fragment Ion'] == fragIon) & 
                        (dataFrame['Precursor Charge'] == precursorCharge) & 
                        (dataFrame['Isotope Label Type'] == isotopeLabel) & 
                        (dataFrame['Product Charge'] == prodCharge) & 
                        (dataFrame['File Name'] == fileName)]
    X = subDF['Times'].values[0].split(',')
    Y = subDF['Intensities'].values[0].split(',')
    
    if smooth > 0:
        Y = vizLib.smoothListGaussian([float(i) for i in Y], degree=smooth)
    ax.plot(X, Y, color=color)
    if zoom:
        if numpy.max(Y) > 100:
            maxIndex = numpy.argmax(numpy.array(Y))
            ax.set_xlim(float(X[maxIndex])-2, float(X[maxIndex])+2)

    return ax
    
def plotMRM_index(dataFrame, index, ax, color='grey', smooth=0, zoom=False):
    [fn, pms, preC, fi, proC, isotope] = getInfoMRM(dataFrame, index)
    plotMRM(dataFrame, fn, pms, preC, fi, proC, isotope, ax, color=color, smooth=smooth, zoom=zoom)
    return ax

def plotLHMRM_index(dataFrame, index, ax, smooth=0):
    [fn, pms, preC, fi, proC, isotope] = getInfoMRM(dataFrame, index)
    indexOpp = getPairedIonMRM(dataFrame, index)
    if isotope == 'heavy':
        [index, indexOpp] = [indexOpp, index]
    ax = plotMRM_index(dataFrame, index, ax, color='blue', smooth=smooth)
    ax = plotMRM_index(dataFrame, indexOpp, ax, color='red', smooth=smooth)
    return ax

def plotLHMRM(dataFrame, fileName, pepSeq, precursorCharge, prodCharge, fragIon, ax, smooth=0):
    ax = plotMRM(dataFrame, fileName, pepSeq, precursorCharge, prodCharge, fragIon, 'light', ax, color='blue', smooth=smooth)
    ax = plotMRM(dataFrame, fileName, pepSeq, precursorCharge, prodCharge, fragIon, 'heavy', ax, color='red', smooth=smooth)
    return ax

def prettyPlot3TransMRM(dataFrame, toPlotIndex, figsize=(33,8.5)):
    axisArray = []
    pylab.figure(figsize=figsize)
    for i in range(len(toPlotIndex)):
        a = plt.subplot2grid((2,9), (0, i*3), colspan=1, rowspan=2)
        a = plotTotalChromPairsMRM(dataFrame, toPlotIndex[i], a, smooth=5, zoom=True)
        info = getInfoMRM(dataFrame, toPlotIndex[i])
        a.set_title(info[1])
        vizLib.tickNum(a, xAxis=4, yAxis=4)
        vizLib.cleanAxis(a, ticksOnly=True)
        axisArray.append(a)
    
        b = plt.subplot2grid((2,9), (0,i*3+1), colspan=2, rowspan=1)
        b = plotAllTransitionsMRM(dataFrame, toPlotIndex[i], b, smooth=5, zoom=True)
        info = getInfoMRM(dataFrame, toPlotIndex[i])
        b.set_title('light')
        vizLib.tickNum(b, xAxis=4, yAxis=4)
        vizLib.cleanAxis(b, ticksOnly=True)
        axisArray.append(b)
    
        c = plt.subplot2grid((2,9), (1,i*3+1), colspan=2, rowspan=1)
        opp = getPairedIonMRM(dataFrame, toPlotIndex[i])
        c = plotAllTransitionsMRM(dataFrame, opp, c, smooth=5, zoom=True)
        info = getInfoMRM(dataFrame, toPlotIndex[i])
        c.set_title('heavy')
        vizLib.tickNum(c, xAxis=4, yAxis=4)
        vizLib.cleanAxis(c, ticksOnly=True)
        axisArray.append(c)
    pylab.tight_layout()
    return axisArray

def readMRMCSV(path, l = 'light ', h = 'heavy ', fileNameHeader = 'File Name', proteinNameHeader = 'Protein Name'):
    fileName = path.split('/')[-1].split('.')
    dataFrame = pandas.read_csv(path)
    s='_'
    dataFrame['shortName'] = dataFrame[fileNameHeader].str.split('.').str[0]
    dataFrame['UID'] =  dataFrame['shortName'] +s+ dataFrame[proteinNameHeader] +s+ dataFrame['Begin Pos'].map(str) +s+\
                        dataFrame['End Pos'].map(str) +s+ dataFrame[l + 'Precursor Mz'].map(str).str.split('.').str[0] +s+ \
                        dataFrame['Product Charge'].map(str) +s+ dataFrame['Fragment Ion'].str[-3:]
    dataFrame['TID'] =  dataFrame[proteinNameHeader] +s+ dataFrame['Begin Pos'].map(str) +s+\
                        dataFrame['End Pos'].map(str) +s+ dataFrame[l + 'Precursor Mz'].map(str).str.split('.').str[0] +s+ \
                        dataFrame['Product Charge'].map(str) +s+ dataFrame['Fragment Ion'].str[-3:]
    dataFrame['PID'] =  dataFrame['shortName'] +s+ dataFrame[proteinNameHeader] +s+ dataFrame['Begin Pos'].map(str) +s+\
                        dataFrame['End Pos'].map(str) +s+ dataFrame[l + 'Precursor Mz'].map(str).str.split('.').str[0]
    dataFrame[l+'AdjArea'] = dataFrame[l+'Area']
    dataFrame[h+'AdjArea'] = dataFrame[h+'Area']
    dataFrame['currentCalc'] = calcValue(dataFrame, [l], [h])
    dataFrame['ratio'] = calcValue(dataFrame, [l],[h])
    dataFrame = dataFrame[pandas.notnull(dataFrame[fileNameHeader])]
    
    positionOtherDict = {key:int(value)+1 for value, key in enumerate(qMS.sort_nicely(sorted(set(dataFrame[proteinNameHeader].values))))}
    positionLookupOther = pandas.Series(positionOtherDict)
    dataFrame['otherpos']=positionLookupOther[dataFrame[proteinNameHeader]].values
    dataFrame['handDelete'] = False
    dataFrame['handSave'] = False
    dataFrame['PPMtranLH'] = abs(dataFrame[l+'Mass Error PPM'] - dataFrame[h+'Mass Error PPM'])
    dataFrame['PPMtranTRANALL_light'] = abs(dataFrame[l+'Mass Error PPM'] - dataFrame[l+'Average Mass Error PPM'])
    dataFrame['PPMtranTRANALL_heavy'] = abs(dataFrame[h+'Mass Error PPM'] - dataFrame[h+'Average Mass Error PPM'])
    dataFrame['PPMtranTRANALL'] = dataFrame[['PPMtranTRANALL_light', 'PPMtranTRANALL_heavy']].max(axis=1)
    
    dataFrame['RTdsLH'] = abs(dataFrame[l+'Retention Time'] - dataFrame[h+'Retention Time'])
    dataFrame['RTdsTRANALL_light'] = abs(dataFrame[l+'Retention Time'] - dataFrame[l+'Best Retention Time'])
    dataFrame['RTdsTRANALL_heavy'] = abs(dataFrame[h+'Retention Time'] - dataFrame[h+'Best Retention Time'])
    dataFrame['RTdsTRANALL'] = dataFrame[['RTdsTRANALL_light', 'RTdsTRANALL_heavy']].max(axis=1)
    
    dataFrame['RTpepTRANALL_light'] = abs(dataFrame[l+'Retention Time'] - dataFrame['Average Measured Retention Time'])
    dataFrame['RTpepTRANALL_heavy'] = abs(dataFrame[h+'Retention Time'] - dataFrame['Average Measured Retention Time'])
    dataFrame['RTpepTRANALL'] = dataFrame[['RTpepTRANALL_light', 'RTpepTRANALL_heavy']].max(axis=1)
    if not 'priorFilter' in dataFrame.columns:
        dataFrame['priorFilter'] = True
    if not 'allClear' in dataFrame.columns:
        dataFrame['allClear'] = dataFrame['priorFilter']
    
    for p in list(dataFrame[proteinNameHeader].unique()):
        peptides = list(dataFrame.loc[dataFrame[proteinNameHeader] == p, 'Peptide Modified Sequence'].unique())
        for n, pep in enumerate(peptides):
            l = float(len(peptides))
            ind = n/l
            dataFrame.loc[(dataFrame['Peptide Modified Sequence'] == pep) & (dataFrame[proteinNameHeader] == p),'colOff'] = ind
    
    dataFrame['currentPosDataset']=dataFrame['otherpos']
    positionFileDict = {key:int(value)+1 for value,key in enumerate(qMS.sort_nicely(sorted(dataFrame[fileNameHeader].unique())))}
    positionFileLookup = pandas.Series(positionFileDict)
    dataFrame['currentPosProtein']=positionFileLookup[dataFrame[fileNameHeader]].values
    dataFrame.loc[dataFrame['currentCalc'] == numpy.inf, 'allClear'] = False
    return dataFrame

def calcValue(df, num, den, field = 'AdjArea', offset=0.0, func=qMS.unity):
    nsDF = df[num[0]+field]
    for x in num[1:]:
        nsDF = nsDF + df[x+field]
    if den == 'label-free':
        value = nsDF + offset    
    else:
        dsDF = df[den[0]+field]
        for x in den[1:]:
            dsDF = dsDF + df[x+field]
        try:
            value = nsDF/dsDF + offset
        except TypeError:
            print "Error in calculating values - some entry must contain strings." 
            print "This can be fixed by deleting this row in vi (you'll see a bunch of NaN values there)." 
            print "Until this is fixed, all values set to 0.0" 
            value = -10.0
    return func(value)

def calcAdjValue(df, isotopeLabels=['light ', 'heavy '], field='Area'):
    for i in isotopeLabels:
        df.loc[:, i+'AdjArea'] = df[i+field]
    return df

def calcMRMTotalProtOcc(dfTotal, proteins, files, num=['light '], den=['heavy ']):
    dfTotal.loc[:,'currentCalc'] = calcValue(dfTotal, num, den)
    d = {}
    for i in proteins:
        d[i]={}
        for j in files:
            peps = list(set(dfTotal[(dfTotal['Protein Name']==i) & (dfTotal['File Name']==j)]['Peptide Modified Sequence']))
            d[i][j]=[]
            for k in peps:
                d[i][j].append(dfTotal[(dfTotal['Protein Name']==i) & 
                                       (dfTotal['File Name']==j) & 
                                       (dfTotal['Peptide Modified Sequence']==k)]['currentCalc'].values[0])
    return d

def calcForUniquePeps(mdf, num = ['light'], den = ['heavy']):
    assert len(mdf['File Name'].unique()) == 1, 'Please pass a datafram with only one file name to calcForUniquePeps'
    
    ups = mdf['Peptide Modified Sequence With Charge'].unique()
    upsDict = {}
    
    for u in ups:
        upsDict[u] = {'precursorCalc':numpy.nan, 'productCalc':numpy.nan, 'totalCalc':numpy.nan}
    for pep in ups:
        numTotal = [0.0,0.0,0.0]
        denTotal = [0.0,0.0,0.0]
        for n in num:
            numTotal[0] = float(numTotal[0] + mdf.loc[(mdf['Peptide Modified Sequence With Charge'] == pep) & 
                                                      (mdf['Fragment Ion Type'] == 'precursor'), n+' Area'].sum())
            numTotal[1] = float(numTotal[1] + mdf.loc[(mdf['Peptide Modified Sequence With Charge'] == pep) & 
                                                      (mdf['Fragment Ion Type'] != 'product'), n+' Area'].sum())
        numTotal[2] = numTotal[0]+numTotal[1]
        for d in den:
            denTotal[0] = float(denTotal[0] + mdf.loc[(mdf['Peptide Modified Sequence With Charge'] == pep) & 
                                      (mdf['Fragment Ion Type'] == 'precursor'), d+' Area'].sum())
            denTotal[1] = float(denTotal[1] + mdf.loc[(mdf['Peptide Modified Sequence With Charge'] == pep) & 
                                      (mdf['Fragment Ion Type'] != 'product'), d+' Area'].sum())
        denTotal[2] = denTotal[0]+denTotal[1]
        try:
            upsDict[pep]['precursorCalc'] = numTotal[0]/denTotal[0]
        except ZeroDivisionError:
            upsDict[pep]['precursorCalc'] = numpy.nan
        try:
            upsDict[pep]['productCalc'] = numTotal[1]/denTotal[1]
        except ZeroDivisionError:
            upsDict[pep]['productCalc'] = numpy.nan
        try:
            upsDict[pep]['totalCalc'] = numTotal[2]/denTotal[2]
        except ZeroDivisionError:
            upsDict[pep]['totalCalc'] = numpy.nan
    return upsDict
    
def getValByProt(mdf, proteinName, upsDict, averageMethod = 'median', field = 'productCalc', allValues = False):
    assert len(mdf['File Name'].unique()) == 1, 'Please pass a dataframe with only a single file name to getValByProt'
    assert (averageMethod == 'median' or averageMethod == 'mean')
    
    pepsToAvg = mdf[mdf['Protein Name'] == proteinName]['Peptide Modified Sequence With Charge'].unique()
    valList = []
    for p in pepsToAvg:
        valList.append(upsDict[p][field])
    
    if averageMethod == 'median':
        toReturn = numpy.nanmedian(valList)
    else:
        toReturn = numpy.nanmean(valList)
        
    if allValues:
        toReturn = [toReturn]
        toReturn.append(valList)
    else:
        return toReturn
    
def calcFullOccupancyHeatMap(df, allFiles = None, allProteins = None, num = ['light'], den = ['heavy'], averageMethod = 'median', field = 'productCalc', allValues = False, nameField='File Name'):
    if allProteins == None:    
        allProteins = df['Protein Name'].unique()
    if allFiles == None:
        allFiles = df[nameField].unique()
    prot_file_dict = {}
    for p in allProteins:
        prot_file_dict[p] = {}
        for file_dict in allFiles:
            prot_file_dict[p][file_dict] = [numpy.nan]
    for f in allFiles:
        subDF = df[df[nameField] == f]
        upCalc = calcForUniquePeps(subDF)
        for p in allProteins:
            prot_file_dict[p][f] = getValByProt(subDF, p, upCalc)
    prot_file_dataFrame = pandas.DataFrame(prot_file_dict)
    return prot_file_dataFrame

def filterDF(df, sigmas=2, ppmDiff=15, RTDiff=0.5, ldp=0.2, lhDP=0.6, nameField='File Name'):
    df['Peptide Modified Sequence With Charge'] = df['Peptide Modified Sequence'] + '[+' + df['Precursor Charge'].map(str) + ']'
    fileNames = qMS.sort_nicely(list(df[nameField].unique()))
    for fn in fileNames:
        uH = df[df[nameField] == fn]['heavy Mass Error PPM'].mean()
        uL = df[df[nameField] == fn]['light Mass Error PPM'].mean()
        sH = df[df[nameField] == fn]['heavy Mass Error PPM'].std()
        sL = df[df[nameField] == fn]['light Mass Error PPM'].std()

        df.loc[(df[nameField] == fn), 'hPPM'] = abs(df.loc[(df[nameField] == fn), 'heavy Mass Error PPM'] - uH) < sigmas * sH
        df.loc[(df[nameField] == fn), 'lPPM'] = abs(df.loc[(df[nameField] == fn), 'light Mass Error PPM'] - uL) < sigmas * sL

        df.loc[(df[nameField] == fn), 'ppmDiff'] = abs(df.loc[(df[nameField] == fn), 'heavy Mass Error PPM'] - 
                                                      df.loc[(df[nameField] == fn), 'light Mass Error PPM']) < ppmDiff

        df.loc[(df[nameField] == fn), 'RTDiff'] = abs(df.loc[(df[nameField] == fn), 'light Retention Time'] - 
                                                      df.loc[(df[nameField] == fn), 'heavy Retention Time']) < RTDiff

        df.loc[(df[nameField] == fn), 'hLDP'] = df.loc[(df[nameField] == fn), 'heavy Library Dot Product'] > ldp
        df.loc[(df[nameField] == fn), 'lLDP'] = df.loc[(df[nameField] == fn), 'light Library Dot Product'] > ldp
    
        df.loc[(df[nameField] == fn), 'dpLTH'] = df.loc[(df[nameField] == fn), 'DotProductLightToHeavy'] > lhDP
    
        df.loc[(df[nameField] == fn), 'allClear'] = ((df.loc[(df[nameField] == fn), 'hPPM']) &
                                                (df.loc[(df[nameField] == fn), 'lPPM']) &
                                                (df.loc[(df[nameField] == fn), 'ppmDiff']) &
                                                (df.loc[(df[nameField] == fn), 'RTDiff']) &
                                                (df.loc[(df[nameField] == fn), 'hLDP']) & 
                                                (df.loc[(df[nameField] == fn), 'lLDP']) & 
                                                (df.loc[(df[nameField] == fn), 'dpLTH']))
    return df