# -*- coding: utf-8 -*-
"""
Created on Sun Jun  8 10:25:03 2014

@author: jhdavis
"""

import qMS
import glob
import sys
import pandas as pd
pd.options.mode.chained_assignment = None
processDirectory = str(sys.argv[1])
listToProcess=glob.glob(processDirectory+'/*_iso.csv')
for f in listToProcess:
    print 'processing' + f
    sys.stdout.flush()
    qMS.preProcessIsoCSV(f, True)
