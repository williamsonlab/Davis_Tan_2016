# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 16:30:13 2015

@author: joeydavis <joeydavis@gmail.com>
"""

import pandas
import qMS
import numpy
import scipy
import sys


def readFileListing(fileListString):
    fileList = [i.rstrip().split(',') for i in open(fileListString).readlines()]
    return fileList

def readSingleIsoCSV(path, shortName=None):
    dataset = qMS.readIsoCSV(path)
    if shortName is None:
        shortName = path
    dataset['shortName'] = shortName

    s = '_'
    dataset['UID'] =  dataset['originFile'] +s+ dataset['protein'] +s+ dataset['startres'].map(str).str.split('.').str[0] +s+\
                        dataset['endres'].map(str).str.split('.').str[0] +s+ dataset['mw'].map(str).str.split('.').str[0] +s+ \
                        dataset['charge'].map(str).str.split('.').str[0] +s+ dataset['rt_n14'].map(str).str.split('.').str[0]
    dataset['TID'] =  dataset['protein'] +s+ dataset['startres'].map(str).str.split('.').str[0] +s+\
                        dataset['endres'].map(str).str.split('.').str[0] +s+ dataset['mw'].map(str).str.split('.').str[0] +s+ \
                        dataset['charge'].map(str).str.split('.').str[0]
    return dataset

if __name__ == '__main__':
    if len(sys.argv) > 1:
        fileListString=sys.argv[1]
    else:
        fileListString='/home/jhdavis/data/emDeps/emrSet/testFileList.txt'

    if len(sys.argv) > 2:
        outputFileString=sys.argv[2]
    else:
        outputFileString='/home/jhdavis/data/emDeps/emrSet/testOutput_mergedCSV.csv'
        
    fileList = readFileListing(fileListString)
    
    fullDataset = pandas.DataFrame()
    for f in fileList:
	print f
        sn = f[1]
        fullDataset = fullDataset.append(readSingleIsoCSV(f[0], shortName=sn), ignore_index=True)
    fullDataset.to_csv(outputFileString, index=False)
