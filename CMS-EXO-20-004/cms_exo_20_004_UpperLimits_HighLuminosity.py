#!/usr/bin/env python3

import numpy as np
import pandas as pd
import time,os,sys
import progressbar as P
sys.path.append('../')
from statisticalTools.simplifiedLikelihoods import Data,UpperLimitComputer



def computeULs(inputFile,outputFile,deltas=0.0):

    metBins = [250,  280,  310,  340,  370,  400,  430,  470,  510, 550,  590,  640,  690,  
                740,  790,  840,  900,  960, 1020, 1090, 1160, 1250, 1400]


    # ### Load CMS data
    cmsYieldFile = './CMS_data/HEPData-ins1894408-v2-csv/Simplifiedlikelihood:Yields(Monojet).csv'
    bgYields = np.genfromtxt(cmsYieldFile,delimiter=',',skip_header=5,
                            skip_footer=67,names=True,dtype=None,encoding=None)

    covarianceFile = './CMS_data/HEPData-ins1894408-v2-csv/Simplifiedlikelihood:covariancematrix(Monojet).csv'
    cov = np.genfromtxt(covarianceFile,delimiter=',',skip_header=5,
                            names=True,dtype=None,encoding='UTF-8')


    lumi_2018 = 59.7
    lumi = 3000.0
    # ### For simplicity only rescale the 2018 data set
    srOrder = cov['Second_Bin'][:66]
    srOrder = np.array([sr for sr in srOrder[::] if '2018' in sr])
    covMatrix = np.zeros((len(srOrder),len(srOrder)))
    for pt in cov:
        if pt['First_Bin'] not in srOrder:
            continue
        if pt['Second_Bin'] not in srOrder:
            continue
        i = np.where(srOrder == pt['First_Bin'])
        j = np.where(srOrder == pt['Second_Bin'])
        covMatrix[i,j] = pt['Covariance']*((lumi/lumi_2018)**2) # Recale by lumi**2


    # ### Get number of observed and expected (BG) events
    nbg = [float(bgYields[bgYields['Bin'] == sr]['Background_yield'][0])*(lumi/lumi_2018) for sr in srOrder] # Rescale by lumi
    nobs = [int(n) for n in nbg] # Assume nobs = nbg


    # ### Load Recast Data
    recastData = pd.read_pickle(inputFile)


    # ### Get all model points
    models = []

    mCols = ['mLLP','mLSP']
    
    for row in recastData[mCols].values:
        m = dict(zip(mCols,row.tolist()))
        if m not in models:
            models.append(m)
    
    # ### Loop over model points and compute UL on mu
    ulComp = UpperLimitComputer()

    progressbar = P.ProgressBar(widgets=["Computing upper limits: ", P.Percentage(),
                                P.Bar(marker=P.RotatingMarker()), P.ETA()])
    progressbar.maxval = len(models)
    progressbar.start()


    for im,m in enumerate(models):

        progressbar.update(im)

        dfModel = recastData.loc[(recastData[list(m)] == pd.Series(m)).all(axis=1)]
        ns = []
        for sr in srOrder:
            dataset,ibin = sr.split('_')[1:]
            dataset = eval(dataset)
            ibin = eval(ibin.replace('bin',''))
            binLabel = 'bin_%1.1f_%1.1f' %(metBins[ibin],metBins[ibin+1])
            binLumi = dfModel[dfModel['Data-takingperiod'] == dataset]['Luminosity (1/fb)'].iloc[0]
            signalYield = dfModel[dfModel['Data-takingperiod'] == dataset][binLabel].iloc[0]
            signalYield *= (lumi/lumi_2018)  # Rescale by lumi
            ns.append(signalYield)
        ns = np.array(ns)    
        data = Data(observed=nobs, backgrounds=nbg, covariance=covMatrix, 
                    nsignal=ns,deltas_rel=deltas)
        dataExp = Data(observed=[int(b) for b in nbg], backgrounds=nbg, covariance=covMatrix, 
                    nsignal=ns,deltas_rel=deltas)    
        try:
            ul = ulComp.getUpperLimitOnMu(data)
            recastData.loc[dfModel.index,'robs'] = 1./ul
        except:
            print('Error computing ul for model:\n',m,'\n')
            ul = None
            recastData.loc[dfModel.index,'robs'] = 0.0
        try:
            ulExp = ulComp.getUpperLimitOnMu(dataExp)
            recastData.loc[dfModel.index,'rexp'] = 1./ulExp 
        except:
            print('Error computing ulExp for model:\n',m,'\n')
            ulExp = None
            recastData.loc[dfModel.index,'rexp'] = 0.0
        
        
        

    progressbar.finish()
    recastData.to_pickle(outputFile)



if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser(description=
            "Compute the upper limit on mu from CMS-EXO-20-004 for the recast data stored in the input file.")
    ap.add_argument('-f', '--inputFile', required=True,
            help='path to the pickle file containing the Pandas DataFrame with the recasting results for the models')
    ap.add_argument('-e', '--signalError', required=False, default = 0.0, type=float,
            help='value for the relative error on the signal [default = 0.0]')
    ap.add_argument('-o', '--outputFile', required=False,
            help='path to output file. If not defined the upper limits will be stored in the input file.',
            default = None)

    t0 = time.time()

    # # Set output file
    args = ap.parse_args()
    inputFile = args.inputFile
    if not os.path.isfile(inputFile):
        print("File %s not found" %inputFile)
        sys.exit()
    outputFile = args.outputFile
    if outputFile is None:
        outputFile = inputFile.replace('.pcl','_HL.pcl')

    deltas = args.signalError
    computeULs(inputFile,outputFile,deltas)

    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))





