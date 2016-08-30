"""
.. module:: qMS 
   :platform: Any
   :synopsis: A collection of functions for qMS data processing and analysis

.. moduleauthor:: Joey Davis <joeydavis@gmail.com>; Josh Silverman <josh.silverman@gmail.com>

"""

import os
import csv
import sys
import urllib2
import numpy
import qMSDefs
import pandas as pd
import re

#qMS utilities

path = os.getcwd()
sys.setrecursionlimit(10000000)
[startindex, pep_seq_index, exp_mr_index, csvopen] = [0, 0, 0, ''];
referencePath = '/home/jhdavis/scripts/python/modules/qMSmodule/'

#######Nice sorting and printing utilities######
def tryint(s):
    try:
        return int(s)
    except:
        return s
     
def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)
    return l

def calculateMedian(statsDictDict, orderedListOfDictKeys, subunits, defaultValue=0.0):
    proteinList = {}    
    for p in subunits:
            proteinList[p] = []
    
    for n in orderedListOfDictKeys:
        for p in subunits:
            try:
                proteinList[p].append(getProtMedValue(statsDictDict[n], p))
            except KeyError:
                proteinList[p].append(defaultValue)
    return proteinList

def readDataFile(filename, scale=1.0, delimiter='\t'):
    """readDataFile reads a datafile. The datafile should be tab separated and have both column and row headers listing the proteins/fractions.
    Empty values should be empty, they will be treated as np.NAN - see lambda below
    readDataFile takes an optional scale that will multiply the data by the specified factor - useful for "scaling" error

    :param filename: a string of the filename to be read
    :type filename: string
    :param scale: a float indicating how to scale the data, default is unscaled (1)
    :type scale: float
    :returns:  a data ditionary. 'data', 'fractions', 'proteins' are filled in, others are None

    """
    with open(filename, 'rb') as inputFile:
        csvFile = list(csv.reader(inputFile, delimiter = delimiter))
        header = csvFile[0]
        proteins = [row[0] for row in csvFile[1:]]
        #insertValue2 = lambda x, y: numpy.NAN if x=='' else float(x)*y
        insertValue2 = lambda x, y: 0.0 if x=='' else float(x)*y
        data = numpy.array([[insertValue2(col, scale) for col in row[1:]] for row in csvFile[1:]])
        inputFile.close()
        return {'fractions':header[1:], 'proteins':proteins, 'fi':None, 'pi':None, 'data':data}

def readIsoCSV(filename, columns=None, noProcess=False):
    """readIsoCSV takes a filename pointing to a _iso.csv file. It returns the calculated
        pandas dataFrame. Optional argument columns can be used to specify specific column headers

    :param filename: full path to the _iso.csv file
    :type filename: string
    :param columns: optional list of strings with the columns to incorporate into the dataFrame
    :type columns: list of strings
    :returns:  a pandas dataFrame with the relevant contents of the _iso.csv. Function
        automatically determines if the dataset is a pulse (bears AMP_S) or variable labeling (bears FRC_NX)

    """

    if columns is None:
        r = csv.reader(open(filename))
        header = r.next()
        pulse = 'AMP_S' in header
        varLab = 'FRC_NX' in header
        al = 'AMP_L' in header
        allClear = 'allClear' in header
        origin = 'originFile' in header
        TID = 'TID' in header
        UID = 'UID' in header
        priorFilt = 'priorFilter' in header
        allFPass = 'allFPass' in header
        shortName = 'shortName' in header
        
        if noProcess:
            columns=['isofile', 'isoprotein', 'isopep', 'mw', 'isoz_charge', 'tim', 
                     'chisq', 'symb', 'mz', 'B', 'OFF', 'GW', 'AMP_U', 
                     'rt_n14', 'rt_n15', 'mz_n14', 'mz_n15',
                     'ppm_n14', 'ppm_n15', 'n14mass', 'n15mass', 'protein', 'startres',
                     'endres', 'charge', 'missed', 'seq', 'mod', 'seqmod', 'file']
                     
        else:
            columns=['isofile', 'isoprotein', 'isopep', 'mw', 'isoz_charge', 'tim', 
                     'chisq', 'symb', 'mz', 'B', 'OFF', 'GW', 'AMP_U', 
                     'rt_n14', 'rt_n15', 'mz_n14', 'mz_n15',
                     'ppm_n14', 'ppm_n15', 'n14mass', 'n15mass', 'protein', 'startres',
                     'endres', 'charge', 'missed', 'seq', 'mod', 'seqmod', 'file',
                     'resid', 'minIntensity', 'ratio', 'currentCalc',
                     '70Spos', '50Spos', '30Spos', 'otherpos', 'currentPos',
                     'ppmDiff', 'rtDiff', 'handDelete', 'handSave']
        
        
        if pulse:
            columns.append('AMP_S')
        if varLab:
            columns.append('FRC_NX')
        if al:
            columns.append('AMP_L')
        if origin:
            columns.append('originFile')
        if allClear:
            columns.append('allClear')
        if TID:
            columns.append('TID')
        if UID:
            columns.append('UID')
        if priorFilt:
            columns.append('priorFilter')
        if allFPass:
            columns.append('allFPass')
        if shortName:
            columns.append('shortName')
    
    ttof = 'matchmode' in header
    if ttof:
        columns.append('pks_orig_rt_n14')
        columns.append('pks_orig_rt_n15')
    data = pd.read_csv(filename, usecols=columns)
    if not pulse:
        data = data.rename(columns={'AMP_L': 'AMP_S'})

    positionOtherDict = {key:int(value)+1 for value, key in enumerate(sort_nicely(sorted(set(data['protein'].values))))}
    positionLookupOther = pd.Series(positionOtherDict)
    data['70Spos']=qMSDefs.positionLookup70S[data['protein']].values
    data['50Spos']=qMSDefs.positionLookup50S[data['protein']].values
    data['30Spos']=qMSDefs.positionLookup30S[data['protein']].values
    data['otherpos']=positionLookupOther[data['protein']].values
    
    data['currentPos']=data['otherpos']
    data['ppmDiff']=data['ppm_n14'] - data['ppm_n15']
    data['rtDiff']=data['rt_n14'] - data['rt_n15']
    if not allClear:
        data['allClear'] = True
        data['handDelete'] = False
        data['handSave'] = False
    if not priorFilt:
        data['priorFilter'] = True
    if not allFPass:
        data['allFPass'] = True
    if not origin:
        data['originFile'] = filename
    if ttof:
        data['missed'] = data.apply(lambda row: calcMissedTTOF(row), axis=1)
        data['ppm_n14'] = data.apply(lambda row: calcPPMTTOF(row), axis=1)
        data['rtDiff'] = data.apply(lambda row: calcRTDiffTTOF(row), axis=1)
        data['ppm_n15'] = data['ppm_n14']
        data['ppmDiff'] = data['ppm_n14'] - data['ppm_n14'].median()
    return data

def calcRTDiffTTOF(row, origField = 'pks_orig_rt_n14', finalField='rt_n14'):
    return row[origField] - row[finalField]

def calcPPMTTOF(row, field='mz_n14', useCharge=True):
    try:
        return 1e6*(row['OFF']/row[field])
    except TypeError:
        return -999

def calcMissedTTOF(row, cutters = ['R', 'K'], field='seq'):
    """calcMissed looks at the sequence and calculates the number of missed cleavages.
        If there is a non-correct C terminus, it returns -1 (non tryptic, for example)
       
    :param row: a pandas series that must contain a the sequence in the correct field
    :type row: pandas series
    :param cutters: an array with the letters to look for in the sequence as missed cleavages.
    :type cutters: an array of strings
    :param field: a string with the index name for the field to search.
    :type field: string
    
    """
    totalMissed = 0    
    for c in cutters:
        totalMissed+=row[field].count(c)
    return totalMissed-1
    

def subtractDoubleSpike(refDF, dataDF, num='AMP_U', den='AMP_S'):
    """subtractDoubleSpike takes two pandas dataFrames, it first divides AMP_U by AMP_S.
        It then finds the median value for the AMP_U/S field on a protein by protein basis from the reference set.
        This is subtracted from the experimetnal fields, and they are mupliplied back by AMP_S to give the proper value.
        A corrected DoubleSpikeDF is returned

    :param refDF: a pandas dataFrame as calculated from an _iso.csv file (qMS.readIsoCSV).
    :type refDF: pandas dataFrame.
    :param dataDF: a list of strings identifying the numerator (must be AMP_U, AMP_L, and/or AMP_S).
    :type dataDF: list of strings.
    :returns:  a pandas DF that has been corrected for a double spike. Should look the same as the input
        dataDF with the AMP_U fields now corrected.

    """

    dataDF[num] = dataDF[num]/dataDF[den]
    refDF[num] = refDF[num]/refDF[den]
    
    for i in dataDF.index:
        p = dataDF.ix[i]['protein']
        med = refDF[refDF['protein']==p][num].median()
        dataDF.ix[i,num] = max(dataDF.ix[i][num] - med, 0.0)
    dataDF[num] = dataDF[num]*dataDF[den]
    refDF[num] = refDF[num]*refDF[den]
    
    return dataDF

def correctFileForDoubleSpike(expPath, refDF=None, refPath=None, num='AMP_U', den='AMP_S'):
    """correctFileForDoubleSpike takes a path to a dataset to be corrected and either a reference dataframe
        or a string pointing to the iso_csv file for the reference set (the double spike alone dataset).
        It makes a pandas DF out of the experimental data, corrects them for the double spike, and returns 
        the corrected pandas dataframe

    :param expPath: a path to the experimental dataset
    :type expPath: strng
    :param refDF: a pandas dataframe of the reference dataset
    :type refDF: pandas dataframe
    :param refPath: a string pointing to the reference datapath
    :type refPath: a string
    :returns:  a pandas DF that has been corrected for a double spike.

    """
    if refDF is None:
        refDF = readIsoCSV(refPath)
    currentDF = readIsoCSV(expPath)
    currentDF = subtractDoubleSpike(refDF, currentDF, num=num, den=den)
    return currentDF

def correctListOfFiles(refPath, listOfFiles, extension=None, savePath=None, num='AMP_U', den='AMP_S'):
    """correctListOfFiles takes a path to a reference set (double spike alone) and a list of paths to the files to be corrected.
        It makes pandas DFs out of the list of files, corrects them for the double spike, and returns a dictionary
        with the keys as the path to the file and the value as the corrected pandas dataframe

    :param refPath: a path to the reference dataset
    :type refPath: strng
    :param listOfFiles: a list of paths to the files to be corrected
    :type listOfFiles: list of strings.
    :returns:  a dictionary of a pandas DF that has been corrected for a double spike. Each key is the dict
        provided in the listOfFiles

    """
    refDF = readIsoCSV(refPath)
    DFDict = {}
    for n,i in enumerate(listOfFiles):
        currentDF = correctFileForDoubleSpike(i, refDF=refDF, refPath=None, num=num, den=den)
        DFDict[i] = currentDF.copy()
        if not (extension is None):
            if (savePath is None):
                currentDF.to_csv(i+extension, index=False)
            else:
                fileName = i.split('/')[-1]
                currentDF.to_csv(savePath+fileName+extension, index=False)
    return DFDict

def generateDFDict(listOfFiles):
    """generateDFDict takes list of paths to the files that will generate the dictionary of dataframes.
        It makes pandas DFs out of the list of files, and returns a dictionary
        with the keys as the path to the file and the value as the pandas dataframe

    :param listOfFiles: a list of paths to the files to be corrected
    :type listOfFiles: list of strings.
    :returns:  a dictionary of a pandas DF. Each key is the dict
        provided in the listOfFiles

    """
    DFDict = {i:readIsoCSV(i) for i in listOfFiles}
    return DFDict

def openListOfFiles(listOfFiles):
    """openListOfFiles takes a list of files and resturns a dictionary of pandas dataframes with the contents

    :param listOfFiles: a list of paths to the files to be corrected
    :type listOfFiles: list of strings.
    :returns:  a dictionary of a pandas DF that has been corrected for a double spike. Each key is the dict
        provided in the listOfFiles

    """
    DFDict = {}
    for n,i in enumerate(listOfFiles):
        currentDF = readIsoCSV(i)
        DFDict[i] = currentDF.copy()
    return DFDict

def unity(x):
    return x      

def calcStatsDict(dataFrame, numerator, denominator, normalization=1.0, offset=0.0, func=unity, adjProt=None):
    """calcStatsDict takes a dataFrame and , a numerator, a denominator, an offset (applied to value first)
        and a normalization factor (scaling factor applied last). It returns a dictionary with keys of 
        protein names and values as a numpy array of calculated values based on numerator and denominator keys.

    :param dataFrame: a pandas dataFrame as calculated from an _iso.csv file (qMS.readIsoCSV)
    :type dataFrame: pandas dataFrame
    :param numerator: a list of strings identifying the numerator (must be AMP_U, AMP_L, and/or AMP_S)
    :type numerator: list of strings
    :param denominator: a list of strings identifying the numerator (must be AMP_U, AMP_L, and/or AMP_S)
    :type denominator: list of strings
    :param normalization: a float normalization factor if you want to scale all of the values uniformly (applied last)
    :type normalization: float
    :param offset: a float offset factor if you want to offset the values (applied first)
    :type offset: float
    :returns:  a dictionary of numpy arrays. First key is the protein name (one of those given in protein_set).
        This leads to a numpy array with the list of values

    """
    ps = list(set(dataFrame['protein'].values))
    ps.sort()
    toReturn = {p:calcValue(dataFrame[dataFrame['protein']==p], numerator, denominator, offset=offset, func=func).values*normalization for p in ps}
    if not (adjProt is None):
        toReturn[adjProt[0]] = toReturn[adjProt[0]]*adjProt[1]
    return toReturn

def multiStatsDictFromDF(dFDict, num, den, namesList=None, normalization=1.0, offset=0.0, normProtein=None, adjProt=None):
    """multiStatsDictFromDF takes a list of dataframesand a list of nums and dens.
        It returns a dict of dict (first key is the file name) - this leads to a statsDict
        that is keyed by protein names. All of the statsDicts contain a full compliment of keys 
        (based on the file list), with empty numpy arrays if there were no values in the original
        dataset

    :param dFDict: a dict of pandas dataframes (keys are the names of the dataframes - output of qMS.correctListOfFiles)
    :type dFDict: a dict of pandas dataframes
    :param num: a list of strings identifying the numerator (must be AMP_U, AMP_L, and/or AMP_S)
    :type num: list of strings
    :param den: a list of strings identifying the numerator (must be AMP_U, AMP_L, and/or AMP_S)
    :type den: list of strings
    :param namesList: a list of strings identifying the keys in the dfDict to use
    :type namesList: list of strings
    :param normalization: a float normalization factor if you want to scale all of the values uniformly
    :type normalization: float
    :param offset: a float offset factor if you want to alter the values uniformly
    :type offset: float
    :param normProtein: string of the protein to normalize to (will be to the median)
    :type normProtein: string
    :returns:  a dictionary of of dicationaries of numpy arrays. First key is the file name, this leads
        to a dictionary where the first key is the protein name (one of those given in protein_set).
        This leads to a numpy array with the list of values

    """
    if namesList is None:
        namesList = dFDict.keys()
        
    allPs = calcStatsDict(dFDict[namesList[0]], num, den)
    dFStatsDict = dict()
    for name in namesList:
        df = dFDict[name]
        dFStatsDict[name] = calcStatsDict(df, num, den, normalization=normalization, offset=offset)
        if not (normProtein is None):
            normValue = 1/numpy.median(dFStatsDict[name][normProtein])
            dFStatsDict[name] = calcStatsDict(df, num, den, normalization=normValue, offset=offset)
            
    for name in namesList[1:]:
        allPs = appendKeys(allPs, dFStatsDict[name])
    
    for name in namesList:
        dFStatsDict[name] = appendKeys(dFStatsDict[name], allPs) 

    if not (adjProt is None):
        for name in namesList:
            dFStatsDict[name][adjProt[0]] = dFStatsDict[name][adjProt[0]]*adjProt[1]
    
    return dFStatsDict

def multiStatsDict(isoFileList, num, den, normalization=1.0, offset=0.0, normProtein=None, noProcess=False, adjProt=None):
    """multiStatsDict takes a list of _iso.csv files and a list of nums and dens.
        It returns a dict of dict (first key is the file name) - this leads to a statsDict
        that is keyed by protein names. All of the statsDicts contain a full compliment of keys 
        (based on the file list), with empty numpy arrays if there were no values in the original
        dataset

    :param isoFileList: a list of file paths (full paths with entire file name)
    :type isoFileList: a list of strings
    :param num: a list of strings identifying the numerator (must be AMP_U, AMP_L, and/or AMP_S)
    :type num: list of strings
    :param den: a list of strings identifying the numerator (must be AMP_U, AMP_L, and/or AMP_S)
    :type den: list of strings
    :param normalization: a float normalization factor if you want to scale all of the values uniformly
    :type normalization: float
    :param offset: a float offset factor if you want to alter the values uniformly
    :type offset: float
    :param normProtein: string of the protein to normalize to (will be to the median)
    :type normProtein: string
    :returns:  a dictionary of of dicationaries of numpy arrays. First key is the file name, this leads
        to a dictionary where the first key is the protein name (one of those given in protein_set).
        This leads to a numpy array with the list of values

    """
    dFDict = {}
    for i in isoFileList:
        dFDict[i] = readIsoCSV(i, noProcess=noProcess)
    return multiStatsDictFromDF(dFDict, num, den, namesList=isoFileList, normalization=normalization, offset=offset, normProtein=normProtein, adjProt=adjProt)

def mergeFiles(fileList, numerator, denominator, normProtein=None):
    """mergeFiles takes a list of _iso.csv files and returns a merged statsFile data dictionary structure

    :param fileList: a list of strings with the full path for reach .csv file
    :type fileList: list (of strings)
    :param numerator: a list of strings for what elements in the numerator (ampu, ampl, amps)
    :type numerator: list of strings
    :param denominator: a list of strings for what elements in the numerator (ampu, ampl, amps)
    :type denominator: list of strings
    :param proteinToNormalizeTo: string of what protein to normalize to (defaults to None)
    :type proteinToNormalizeTo: strings
    :returns:  the median value in the 'vals' field for that protein
    
    """

    splitDict = multiStatsDict(fileList, numerator,denominator, normProtein=normProtein)
    mergedDict = splitDict[sorted(splitDict.keys())[0]]
    for k in sorted(splitDict.keys())[1:]:
        cd = splitDict[k]
        for p in sorted(cd.keys()):
            mergedDict[p] = numpy.concatenate((mergedDict[p],cd[p]))
    return mergedDict
    


def getProtMedValue(statsFileDict, proteinName):
    """getProtMedValue gets the median value of proteins values from a StatsFile-style dictionary

    :param statsFileDict: a statsFile-style dictionary (can be generated by calcStatsDict 
        (which takes a pandas dataFrame from readIsoCSV))
    :type statsFileDict: dictionary (statsFilte type (needs keys of protein names (strings) and values of numpy arrays))
    :param proteinName: name of the protein to get the med from
    :type proteinName: string
    :returns:  the median value
    
    """
    return numpy.median(statsFileDict[proteinName])

def calcPercent(f, sigfig=2):
    """calcPercent takes a floating point number, rounds it to the number of sigfigs (default 2) and
        returns a string of that number multiplied by 100 and with a ' %'

    :param f: a floating point number to be converted to percentage
    :type f: float
    :param sigfig: the number of sigfigs to round to (default 2)
    :type sigfig: int
    :returns:  a string of the float rounded + ' %'

    """
    s = str(round(f, sigfig)*100)
    return s+" %"

def maxLabFunc(k,t):
    """maxLabFunc is a function to calculate the max labeling based on the equation from Stephen Chen's paper

    :param k: the growth rate (calculated as ln(2)/doubling time)
    :type k: float
    :param t: the time (if a single number, you get back teh max lab at that point), 
        can accept an array and will give back the curve
    :type t: float/int/array
    :returns:  the max lab (array or single value)
    
    """
    return 1 - numpy.exp(-k*t)

def poolFunc(k,t,P):
    """poolFunc is a function to calculate the labeling kinetics for a protein with a given pool size
        Equation from Stephen Chen's paper. USE THIS TO FIT THE LABELING OF TERMINAL (70S) pool.

    :param k: the growth rate (calculated as ln(2)/doubling time)
    :type k: float
    :param t: the time (if a single number, you get back teh max lab at that point), 
        can accept an array and will give back the curve    
    :type t: float/int/array
    :param P: The pool size (expressed as precursorPool/completedRibosomePool)
    :type t: float
    :returns:  the expected labeling (array or single value)
    
    """
    
    return 1.0 + (P*numpy.exp((0.0-k)*(1.0+(1.0/P))*t)) - ((1.0 + P)*numpy.exp((0.0-k)*t))

def poolInterFunc(k,t,P):
    """poolInterFunc is a function to calculate the labeling kinetics for a protein using the 
        overlabeling of an intermediate. Derived from the differential equation in Stephen Chen's paper

    :param k: the growth rate (calculated as ln(2)/doubling time)
    :type k: float
    :param t: the time (if a single number, you get back teh max lab at that point), 
        can accept an array and will give back the curve    
    :type t: float/int/array
    :param d: the turnover rate
    :type d: float
    :returns:  the expected labeling lab (array or single value)
    
    """
    return 1.0 - (numpy.exp((0.0-k)*(1.0+(1.0/P))*t))
    
def poolInterFracXFunc(k,t,P,X=0.65):
    """poolInterFunc is a function to calculate the labeling kinetics for a protein using the 
        overlabeling of an intermediate. Derived from the differential equation in Stephen Chen's paper

    :param k: the growth rate (calculated as ln(2)/doubling time)
    :type k: float
    :param t: the time (if a single number, you get back teh max lab at that point), 
        can accept an array and will give back the curve    
    :type t: float/int/array
    :param d: the turnover rate
    :type d: float
    :returns:  the expected labeling lab (array or single value)
    
    """
    return X*poolInterFunc(k,t,P) + (1.0-X)*poolFunc(k,t,P)

def overLabelingFunc(k,t,d):
    """overLabelingFunc is a function to calculate the labeling kinetics for a protein with a given turnover rate
        Equation from Stephen Chen's paper

    :param k: the growth rate (calculated as ln(2)/doubling time)
    :type k: float
    :param t: the time (if a single number, you get back teh max lab at that point), 
        can accept an array and will give back the curve    
    :type t: float/int/array
    :param d: the turnover rate
    :type d: float
    :returns:  the expected labeling lab (array or single value)
    
    """
    
    return 1.0 - (numpy.exp(-(d+k)*t))

def growthCurve(doublingTime, t, Ai):
    """growthRate calculates a growth rate k given a doubling time (units of k can be used in equations)
        listed in Stephen Chen's paper

    :param doublingTime: the doubling time of the cell in units of time (mins, secs)
    :type doublingTime: floatz
    :returns:  the growth rate k (units of inverse time) - calculated at ln(2)/k
    
    """
    
    return Ai*numpy.exp2(t/float(doublingTime))

def growthRate(doublingTime):
    """growthRate calculates a growth rate k given a doubling time (units of k can be used in equations)
        listed in Stephen Chen's paper

    :param doublingTime: the doubling time of the cell in units of time (mins, secs)
    :type doublingTime: floatz
    :returns:  the growth rate k (units of inverse time) - calculated at ln(2)/k
    
    """
    
    return numpy.log(2)/float(doublingTime)

def addBlankKey(d, k):
    """addBlankKey takes a stats dict and a key. It checks if the key in in the stats file,
        if it is not, it make a new dict entry with the key k and the value as an empty numpy array

    :param d: statsDictionary
    :type d: a dict of numpy arrays - keyed by protein names
    :param k: a string name for the key to check
    :type k: string
    :returns:  statsDictionary (a dictionary of numpy arrays). First key is the protein name (one of those given in protein_set).
        This leads to a numpy array with the list of values

    """
    if not d.has_key(k):
        d[k]=numpy.array([])
    return d
    
def appendKeys(d1, d2):
    """appendKeys a helper function that adds all the keys in d2 to d1 (making empty entries into d1)

    :param d1: statsDictionary
    :type d1: a dict of numpy arrays - keyed by protein names
    :param d2: statsDictionary
    :type d2: a dict of numpy arrays - keyed by protein names
    :returns:  statsDictionary (a dictionary of numpy arrays). First key is the protein name (one of those given in protein_set).
        This leads to a numpy array with the list of values

    """

    for k in d2:
        d1 = addBlankKey(d1, k)
    return d1

def calcValue(df, num, den, offset=0.0, func=unity):
    """calcVale takes a pandas dataFrame bearing the keys to be used and calculates
        the ratio of the num/den (specified as lists of the keys - AMP_U, AMP_L, AMP_S)

    :param df: pandas dataFrame with the information from an _iso.csv
    :type df: pandas dataFrame
    :param num: a list of the strings in the numerator (AMP_U, AMP_L, AMP_S)
    :type den: list of strings
    :param den: a list of strings identifying the denominator (AMP_U, AMP_L, AMP_S)
    :type den: list of strings
    :param offset: an optional offset float to move the point up or down
    :type offset: float
    :returns:  a float of the [(sum of the nums) divided by the (sum of the dens)] + the offset
    
    """
    nsDF = df[num[0]]
    dsDF = df[den[0]]
    for x in num[1:]:
        nsDF = nsDF + df[x]
    for x in den[1:]:
        dsDF = dsDF + df[x]
    try:
        value = nsDF/dsDF + offset
    except TypeError:
        print "Error in calculating values - some entry must contain strings." 
        print "This can be fixed by deleting this row in vi (you'll see a bunch of NaN values there)." 
        print "Until this is fixed, all values set to 0.0" 
        value = 0.0
    return func(value)

def boolParse(s):
    """boolParse takes a string and returns a bool. Any capitilization of "true" results in True
        all other strings result in False

    :param s: string to process
    :type s: string
    :returns:  bool (see description)
    
    """
    return s.upper()=='TRUE'

def preProcessIsoCSV(isoPath, genPlots=True):
    """preProcessIsoCSV processes an _iso.csv and _plots directory to append the 
        residual column as well as a "currentCalc" column showing protein levels.
        Function is a helper that calls calcResidual

    :param isoPath: full path to the _iso.csv file
    :type isoPath: string
    :param genPlots: bool indicating if .plots files should be generated
    :type genPlots: list of strings
    :returns: a complteted dataFrame with resids and currentCalc columns. Has
        externality of generating .plots files in the _peaks directory if genPlots is true
    
    """
    print "reading : " + isoPath + "..."
    dataFrame = readIsoCSV(isoPath, noProcess=True)
    dataFrame['currentCalc'] = calcValue(dataFrame, ['AMP_U'], ['AMP_U', 'AMP_S'])
    dataFrame['ratio'] = calcValue(dataFrame, ['AMP_U'], ['AMP_S'])
    rootPath = '/'.join(isoPath.split('/')[:-1])+'/'
    dataFrame = calcResidual(rootPath, dataFrame, genPlots=genPlots)

    fileName = isoPath.split('/')[-1:][0].replace('.', '_res.')
    print "writing : " + fileName + "..."
    dataFrame.to_csv(rootPath+fileName, index=False)
    return dataFrame

def cleanPlotsDir(path, extensions=['.newcon', '.fit', '.png']):
    """cleanPlotsDir executes a system command to remove files from a _peaks directory

    :param path: full path to the _peaks directory (no trailing /)
    :type path: path
    :param extensions: list of strings with extension types
    :type extensions: list of strings
    :returns: no return, outputs the command issued and any output from the command.
        Has externalities of deleting all *.ext files in the path directory
    
    """
    for ext in extensions:
        command = 'rm -rf ' + path + '/*' + ext
        print command
        output = os.popen(command).read()
        print output

def calcResidual(datapath, dataFrame, genPlots=False):
    """calcResidual takes a path to a _peaks directory and a pandas dataFrame containing the contents
        of the _iso.csv file. It appends a column to the dataFrame with the calculated residual
        (the difference between the fit and data in a .dat file) for each peptide. It also appends a
        column to the dataFrame with the max fit intensity
        Optional paramter genPlots will also generate .plots files that can be used to plot the datasets

    :param datapath: A string with the full path to the _plots directory
    :type datapath: string
    :param dataFrame: A pandas dataframe with ['isofile'] at the minimum (must point to the .dat files in 
        the _plots directory)
    :type dataFrame: pandas dataframe
    :param genPlots: Optional boolean telling function if it should write .plots files
    :type genPlots: boolean

    :returns:  the dataFrame modified to include the 'resid' and 'minIntensity' columns. .Dats that cause any errors
        are given fits with constant 666, resids with constant 666, and minIntensity with constant -666.
    
    """
    dataFrame['resid']=0
    dataFrame['minIntensity']=0
    for iso in dataFrame['isofile'].values:
        datFileName = datapath+iso+'.dat'
        try:
            datPd = pd.read_csv(datFileName, names=['offset', 'dat', 'resid'], header=None)
            del datPd['offset']
            datPd['residAdj'] = datPd['resid']+(datPd['dat'].median()-datPd['resid'].median())
            datPd['fit'] = datPd['residAdj']+datPd['dat']
            calcResid = datPd['residAdj'].abs().sum()/min([datPd['fit'].max(), datPd['dat'].max()])
            calcMinIntensity = min([datPd['fit'].max(), datPd['dat'].max()])
        except (IOError, TypeError, Exception) as e:
            print "Error " + e.message + " in " + datFileName
            sys.stdout.flush()
            datPd = pd.DataFrame({'fit' : 666}, index=[0])
            datPd['residAdj'] = 666
            calcResid = 666
            calcMinIntensity = -666
        #row = dataFrame[dataFrame['isofile']==iso]
        #row['resid'] = calcResid
        #row['minIntensity'] = calcMinIntensity
        #dataFrame.update(row)
        rowIX = dataFrame[dataFrame['isofile']==iso].index.values[0]
        dataFrame.loc[rowIX, 'resid'] = calcResid
        dataFrame.loc[rowIX, 'minIntensity'] = calcMinIntensity
        if genPlots:
            datPd.to_csv(datapath+iso+'.plots', index=False)
    return dataFrame

def getRPSeqData(ID):
    """getRPSeqData is a helper function to get the AA and cDNA seqs for ribosomal proteins

    :param ID: string with the geneID
    :type ID: string
    :returns:  a list of strings, first is by cDNA, second by AA;
        
    """
    address = 'http://ribosome.med.miyazaki-u.ac.jp/rpg.cgi?id='+ID+'&mode=seq'
    website = urllib2.urlopen(address)
    website_html = website.read()
    cDNAi = website_html.find('>cDNA Sequence')
    cDNAj = website_html.find('</textarea></td></tr>\n<tr><td align="center" bgcolor="#FFFF80" width="150">Amino Acids Sequence</td><td>')
    cDNA = website_html[cDNAi+74:cDNAj-1]

    AAj = website_html.find('</textarea></td></tr>\n</table>\n<div class="footer">')
    AA = website_html[cDNAj+155:AAj]
    
    return [cDNA, AA] 

def getRPInfo(numberGenes, start=10, baseURL='http://ribosome.med.miyazaki-u.ac.jp/rpg.cgi?mode=gene&id=ECO100'):
    """getRPInfo generates two dictionaries with relevant ribosomal protein information from a database.
        It prefers the url listed above, but can be used with other organisms so long as the
        find commands still work

    :param numberGenes: the number of genes to look for
    :type numberGenes: int
    :param baseURL: a string pointing to the base url, defaults to the japanese database
    :type baseURL: string
    :returns:  a list of dictionaries, first is by geneNames, second by geneProduct;
        each dictionary is a dict of dicts with subkeys size, cDNA, AA and either GP or GN (opposite
        of the base key)
        
    """
    #gs = ['01', '02', '03', '04', '05', '06', '07', '08', '09']
    genes = range(start, start+numberGenes)
    #genes = gs + genes
    base = baseURL
    addys = [base + str(i) for i in genes]
    print addys
    rpdictGN = {}
    rpdictGP = {}
    for address in addys:
        website = urllib2.urlopen(address)
        website_html = website.read()
        gni = website_html.find('>Gene Name</td><td>')
        gn = website_html[gni+19:gni+26]
        gpi = website_html.find('ibosomal protein ')
        gp = website_html[gpi+17:gpi+20]
        gp = gp.upper()
        if gp[-1] == '<':
            gp = gp[:-1]
        if gp == 'ITL':
            gp = 'L9'
        if gp == 'L7/':
            gp = 'L7/L12'
        si = website_html.find('Gene Size [bp]</td><td>')
        s = website_html[si+23:si+27]
        if s[-1] == '<':
            s = s[:-1]
        cDNA, AA = getRPSeqData(address[-8:])
        rpdictGN[gn] = {'GP':gp, 'size':s, 'cDNA':cDNA, 'AA':AA}
        rpdictGP[gp] = {'GN':gn, 'size':s, 'cDNA':cDNA, 'AA':AA}
    return [rpdictGN, rpdictGP]
    

#fetches protein sequence from uniprot ID
def getsequence(uniprot):
	try:
		urlbase = "http://www.uniprot.org/uniprot/"
		req = urllib2.Request(urlbase + uniprot + ".fasta")
		response = urllib2.urlopen(req)
		lines = [x for i, x in enumerate(response) if i > 0]
		return "".join(map(lambda x: x.rstrip("\n"), lines))
	except urllib2.URLError:
		print "No internet connection or bad Uniprot id"
		
#Returns the MAD of a list
def MAD(list):
	temp_median = numpy.median(list)
	temp_absdevs = map(lambda item: abs(item - temp_median), list)
	return numpy.median(temp_absdevs)

def listReplace(l, to, rv):
    """listReplace is a helper function replaces all occurances of to with rv in a list l

    :param l: list to be replaced
    :type l: list
    :param to: item to be replaced
    :type to: string
    :param rv: item to replace with
    :type rv: string
    :returns: a list with all occurances of to replaced with rv
    
    """
    tr = []
    for i in l:
        if i == to:
            tr.append(rv)
        else:
            tr.append(i)
    return tr

def printSortedDict(d):
    """printSortedDict is a helper function to print a dictionary

    :param d: a dictionary to be printed
    :type d: dict
    :returns: a string of the dictionary
    
    """
    k = d.keys()
    k.sort()
    tp = ''
    for i in k:
        tp = tp + str(i) + ":" + str(d[str(i)]) + ", "
    return tp

def calcMW(seq):
    """calcMW is a helper function that calculates the mass of an AA sequence

    :param seq: the sequence to be calculated
    :type seq: string
    :returns: a float with the mass of the sequence
    
    """
    mw = 0.0
    for aa in seq:
        mw = mw+qMSDefs.aaweights[aa.upper()]
    return mw

def dropDuplicatesPandas(df, cols_to_consider=None):
    """dropDuplicatesPandas is a helper to remove duplicated rows based on a set
        of columns that must all contain identical information.

    :param df: the dataframe to be inspected
    :type df: a pandas dataframe
    :param cols_to_consider: the list of columns to be inspected. this defaults
        to an extensive list that is used for the MRM material, although
        technically an optional parameter, this should almost always be passed.
    :type cols_to_consider: a list of strings
    :returns: the unique dataframe
    
    """
    if cols_to_consider is None:
        cols_to_consider=['Protein Name', 'Begin Pos', 'End Pos', 'Missed Cleavages', 
        'Precursor Charge', 'File Name', 'Peptide Modified Sequence', 'Average Measured Retention Time', 
        'light Total Area', 'heavy Total Area', 'light Precursor Mz', 'heavy Precursor Mz']
    grouped = df.groupby(cols_to_consider)
    index = [gp_keys[0] for gp_keys in grouped.groups.values()]
    unique_df = df.reindex(index)
    return unique_df

def print_fullPandas(x):
    """print_fullPandas is a helper function for printing a full dataframe

    :param x: the pandas dataframe to be printed
    :type x: a pandas dataframe
    :returns: no return, simply prints the full dataframe and then resets the
        default pandas print values.
    
    """
    pd.set_option('display.max_rows', len(x))
    pd.set_option('display.max_columns', 200)
    print(x)
    pd.reset_option('display.max_rows')
    pd.reset_option('display.max_columns')

def makeNameConvDict(set1, set2):
    """calcMW is a helper function that generates a name conversion dictionary.
        Two lists of strings are input and a list with the two dictionaries is output.

    :param set1: a list of strings in the fist set
    :type set1: list of strings
    :param set2: a list of strings in the second set
    :type set2: list of strings
    :returns: a list of dictionaries - [0] contains from set1 to set2; 
        [1] contains from set2 to set1
    
    """
    set1toset2 = {}
    set2toset1 = {}
    for i in range(len(set1)):
        set1toset2[set1[i]] = set2[i]
        set2toset1[set2[i]] = set1[i]
    return [set1toset2, set2toset1]

def convertNames(fullpath, oldNames, newNames):
    """convertNames is a helper function that replaces the protein fields in an iso csv files
        with a converted set of names. The old names are given as a list in oldNames.
        The new names are given as a list in newNames. Returns the pandas dataframe (with new names).
        **HAS EXTERNALITY - CREATES A NEW _ISO.CSV FILE IN THE PATH DIRECTORY THAT HAS THE NEW NAMES.

    :param fullpath: a path string to the file to convert
    :type fullpath: string
    :param oldNames: a list of strings of the old names
    :type oldNames: list of strings
    :param newNames: a list of strings of the old names
    :type newNames: list of strings
    :returns: a pandas DF with the corrected names. 
        **HAS EXTERNALITY - CREATES A NEW _ISO.CSV FILE IN THE PATH DIRECTORY THAT HAS THE NEW NAMES.
    
    """
    [oldToNew, newToOld] = makeNameConvDict(oldNames, newNames)
    dataFrame = readIsoCSV(fullpath)
    dataFrame.update(dataFrame['protein'].replace(oldToNew))
    dataFrame.to_csv(fullpath[:-4]+'_newNames.csv', index=False)
    return dataFrame

def concatonateIsoCSVFiles(fileList, outFile='mergedOutTemp.csv'):
    """concatonateIsoCSVFiles is a helper function that merges a series of isocsv files.
        It does outputs these merged files in the working directory as "mergedOutTmp.csv"

    :param fileList: a list of full path strings to the files to be merged
    :type fileList: a list of strings
    :returns:  no return. 
        **HAS EXTERNALITY - CREATES A NEW _ISO.CSV FILE IN THE WORKING DIRECTORY CALLED mergedOutTmp.csv.
    
    """
    initialData = open(fileList[0], 'r').read()
    for i in fileList[1:]:
        with open(i,'r') as f:
            next(f)
            for line in f:
                initialData = initialData + line
        
    fout = open(outFile, "w")
    fout.write(initialData)
    fout.close()

def addProts(addFrom, addTo, output, listOfProts):
    """Adds specific proteins from one iso csv to another

    :param fileList: a list of full path strings to the files to be merged
    :type fileList: a list of strings
    :returns:  no return. 
        **HAS EXTERNALITY - CREATES A NEW _ISO.CSV FILE IN THE WORKING DIRECTORY CALLED mergedOutTmp.csv.
    
    """
    startingData = open(addTo, 'r').read()
    with open(addFrom,'r') as f:
        next(f)
        for line in f:
            for p in listOfProts:
                if p in line:
                    startingData = startingData+line
    fout = open(output, "w")
    fout.write(startingData)
    fout.close()

def readMSSpectraFile(datapath):
    """readMSSpectraFile is a helper function that reads a spectra .txt file (saved in the _plots directory by massacre)

    :param datapath: a full path string to the file to be read
    :type datapath: a sting
    :returns:  a list with arrays for the x values [0], y values [1], and the datapath [3]
    
    """
    data = list(csv.reader( open(datapath, 'rU'), delimiter=' '))
    xs = []
    ys = []
    for line in data:
        xs.append(float(line[0]))
        ys.append(float(line[1]))
    return [xs, ys, datapath]

def offsetLog(x, offset=1.0):
    return numpy.log(x+offset)

def readMascotCSV(path, header=66):
    """readMascotCSV reads a mascot csv output into a pandas dataframe. All columns are read in.
        Also appends two new columns - isotope (determined by findIsotope), and sc (spectral counts for that
        prot_acc, isotope pair)

    :param path: full path to the mascot csv file
    :type path: string
    :param header: number of initial lines to skip
    :type header: int
    
    :returns:  A string telling whether the isotope is "Light", "Heavy", or "Unknown"
    
    """ 
    df = pd.read_csv(path, header=66)
    df['isotope']=df.apply(findIsotope, axis=1)
    g = df.groupby(['prot_acc','isotope'])
    g = df.groupby(['prot_acc','isotope'])
    df1 = df.set_index(['prot_acc','isotope'])
    df1['sc'] = g.size()
    df = df1.reset_index()
    return df

def findIsotope(row):
    """findIsotope determines the isotope labeling of a particular mascot hit using the mascot calculated mr and 
        my calculated mr (uses calcMonoistopicMass). Assumes 0 or 100% 15N, doesn't do other isotopes.

    :param row: the mascot pandas dataframe for that row (must have 'pep_var_mod', 'pep_seq', 'pep_calc_mr')
    :type row: a row of a pandas dataframe generated by qMS.readMascotCSV
    :returns:  A string telling whether the isotope is "Light", "Heavy", or "Unknown"
    
    """ 
    metOx=row['pep_var_mod']
    if metOx is numpy.NAN:
        metOx = 0
    else:
        try: 
            metOx = int(metOx[0])
        except ValueError:
            metOx = 1
            
            
    if abs(calcMonoisotopicMass(row['pep_seq'], N15=False, cysIAA=True, metOx=metOx) - row['pep_exp_mr']) < .1:
        return 'Light'
    elif abs(calcMonoisotopicMass(row['pep_seq'], N15=True, cysIAA=True, metOx=metOx) - row['pep_exp_mr']) < .1:
        return 'Heavy'
    else:
        return 'Unknown'

def outputSpectralCounts(path, pathOut):
    """outputSpectralCounts generates a csv file with the light and heavy spectral counts found in a mascot csv file
        It also returns a dictionary with the same information (light in the 0 position, heavy in the 1 position)

    :param path: full path to the mascot csv file
    :type path: string
    :param pathOut: full path to the desired output csv file
    :type pathOut: string
    :returns:  A string telling whether the isotope is "Light", "Heavy", or "Unknown"
    
    """ 
    outfile = open(pathOut, 'w')
    df = readMascotCSV(path)
    protsDict = {}
    prots = list(set(df['prot_acc'].values))
    for i in prots:
        protsDict[i] = []
        try:
            protsDict[i].append(df[(df['prot_acc']==i) & (df['isotope']=='Light')]['sc'].values[0])
        except:
            protsDict[i].append(0)
        try:
            protsDict[i].append(df[(df['prot_acc']==i) & (df['isotope']=='Heavy')]['sc'].values[0])
        except:
            protsDict[i].append(0)
    outfile.write('prot_acc,Light_sc,Heavy_sc,'+path+'\n')
    for i in sort_nicely(protsDict.keys()):
        strOut = i+','+str(protsDict[i][0])+','+str(protsDict[i][1])+'\n'
        outfile.write(strOut)
    return protsDict
     

def outputDataMatrixFile(filePath, dfStatsDict, keyList, num, den, subunits=None, normalization=1.0, offset=0.0, normProtein=None, 
                         rowNorms=None, delimiter=',', func=unity, adjProt=None):
    """outputDataMatrixFile outputs a datamatrix that can bea easily plotted as a heat map or clustered or etc.
        It uses a dfStatsDict, keyList (the list of the dataset names), proteins to plot, the num, den, nommalization, offset, and
        normProtein if desired

    :param filePath: a full path string to the file to be written
    :type filePath: a sting
    :param dfStatsDict: a dict of pandas dataframes (keys are the names of the dataframes - output of qMS.correctListOfFiles)
    :type dfStatsDict: a dict of pandas dataFrames
    :param keyList: a list of strings with the keys for the dict
    :type keyList: list of strings
    :param subunits: a list of strings of the proteins to be plotted
    :type subunits: list of strings
    :param num: a list of strings identifying the numerator (must be AMP_U, AMP_L, and/or AMP_S)
    :type num: list of strings
    :param den: a list of strings identifying the numerator (must be AMP_U, AMP_L, and/or AMP_S)
    :type den: list of strings
    :param normalization: a float normalization factor if you want to scale all of the values uniformly
    :type normalization: float
    :param offset: a float offset factor if you want to alter the values uniformly
    :type offset: float
    :param normProtein: string of the protein to normalize to (will be to the median)
    :type normProtein: string
    :returns:  a dictionary of of dicationaries of numpy arrays. First key is the file name, this leads
        to a dictionary where the first key is the protein name (one of those given in protein_set).
        This leads to a numpy array with the list of calculated values.
        ***************HAS EXTERNALITY - CREATES A FILE WITH ALL OF THIS DATA**********************
    
    """       
    outFile = open(filePath, 'w')
    header = 'Protein'+delimiter
    statsDictDict = {}
    if rowNorms is None:
        rowNorms = {k:1.0 for k in keyList}
    for i in keyList:
        header = header + i + delimiter
        statsDictDict[i] = calcStatsDict(dfStatsDict[i], num, den, normalization=1.0/rowNorms[i], offset=offset, adjProt=adjProt)
        if not (normProtein is None):
            if normProtein == 'mean':
                normValue = 1/numpy.mean(numpy.array([numpy.median(statsDictDict[i][p]) for p in statsDictDict[i].keys()]))
            elif normProtein == 'median':
                normValue = 1/numpy.median(numpy.array([j for j in numpy.array([numpy.median(statsDictDict[i][p]) for p in statsDictDict[i].keys()]) if j > 0.1]))
            else:
                normValue = 1/numpy.median(statsDictDict[i][normProtein])
            statsDictDict[i] = calcStatsDict(dfStatsDict[i], num, den, normalization=normValue, offset=offset, func=func, adjProt=adjProt)
    outFile.write(header[:-1] + '\n')

    if subunits is None:
        pList = []
        for i in keyList:
            cList = statsDictDict[i].keys()
            for prot in cList:
                pList.append(prot)
        subunits = sort_nicely(list(set(pList)))
    for p in subunits:
        line = p + delimiter
        for k in keyList:
            try:
                line = line + str(numpy.median(statsDictDict[k][p])) + delimiter
            except KeyError:
                line = line + delimiter
        outFile.write(line[:-1] + '\n')
    return statsDictDict

def outputDataMatrixFileFromSDD(filePath, statsDictDict, keyList, subunits=None, normalization=1.0, 
                                offset=0.0, delimiter=',', func=unity):
    outFile = open(filePath, 'w')
    header = 'Protein'+delimiter
    for i in keyList:
        header = header + i + delimiter
    outFile.write(header[:-1] + '\n')

    if subunits is None:
        pList = []
        for i in keyList:
            cList = statsDictDict[i].keys()
            for prot in cList:
                pList.append(prot)
        subunits = sort_nicely(list(set(pList)))
    for p in subunits:
        line = p + delimiter
        for k in keyList:
            try:
                line = line + str(numpy.median(statsDictDict[k][p])) + delimiter
            except KeyError:
                line = line + delimiter
        outFile.write(line[:-1] + '\n')
    return statsDictDict

def calcMonoisotopicMass(seq, N15=False, waterTermini=True, cysIAA=True, metOx=0):
    total = 0.000
    seq.upper()
    if cysIAA:
        total = seq.count('C')*57.02146
    total = total+metOx*15.99491
    for a in seq:
        total = total+qMSDefs.aaweights[a]

    if N15:
        nitDiff = qMSDefs.N15Mass - qMSDefs.N14Mass
        totalNits = 0
        for a in seq:
            totalNits = totalNits+qMSDefs.aaNitrogens[a]
        total = total + totalNits*nitDiff
    if waterTermini:
        total = total+18.010565
    return total

def to_precision(x,p):
    """
    returns a string representation of x formatted with a precision of p

    Based on the webkit javascript implementation taken from here:
    https://code.google.com/p/webkit-mirror/source/browse/JavaScriptCore/kjs/number_object.cpp
    """


    import math
    x = float(x)

    if x == 0.:
        return "0." + "0"*(p-1)

    out = []

    if x < 0:
        out.append("-")
        x = -x

    e = int(math.log10(x))
    tens = math.pow(10, e - p + 1)
    n = math.floor(x/tens)

    if n < math.pow(10, p - 1):
        e = e -1
        tens = math.pow(10, e - p+1)
        n = math.floor(x / tens)

    if abs((n + 1.) * tens - x) <= abs(n * tens -x):
        n = n + 1

    if n >= math.pow(10,p):
        n = n / 10.
        e = e + 1


    m = "%.*g" % (p, n)

    if e < -2 or e >= p:
        out.append(m[0])
        if p > 1:
            out.append(".")
            out.extend(m[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p -1):
        out.append(m)
    elif e >= 0:
        out.append(m[:e+1])
        if e+1 < len(m):
            out.append(".")
            out.extend(m[e+1:])
    else:
        out.append("0.")
        out.extend(["0"]*-(e+1))
        out.append(m)

    return "".join(out)

def makeStatsDictMRM(df, pList=[], fList=[], fileNameField='File Name', proteinNameField = 'Protein Name', 
                     outputField='RatioLightToHeavy', median=True):
    statsDict = {}
    for prot in pList:
        statsDict[prot]={}
        for f in fList:
            hold = df[(df[fileNameField]==f) & (df[proteinNameField]==prot)][outputField]
            if median:
                statsDict[prot][f]=hold.median()
            else:
                statsDict[prot][f]=list(hold.values)
    return statsDict
  
def saveRPOnly(inputFile, outputFile, stringContains=['RS', 'RL'], field='protein'):
    init = readIsoCSV(inputFile)
    goodIdx = init[field].str.contains('|'.join([i for i in stringContains]))
    init[goodIdx].to_csv(outputFile, index=False)
 
def filterCSV(inputCSVPath, outputCSVPath, regExp):
    ##Example usage: regExp = 'R[SL][0-9]' will give ribosomal proteins
    f = open(inputCSVPath, 'r')
    outFile = open(outputCSVPath, 'w')
    lines = f.readlines()
    outFile.write(lines[0])
    for l in lines[1:]:
        protein = l.split(',')[0].split('/')[1].split('_')[0]
        if bool(re.search(regExp, protein)):
            outFile.write(l)
    outFile.close() 
