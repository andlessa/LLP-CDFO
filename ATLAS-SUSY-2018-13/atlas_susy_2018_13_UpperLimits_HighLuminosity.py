#!/usr/bin/env python3

import pandas as pd
import time,os,sys
import pandas as pd
import numpy as np
import sys
sys.path.append('../')
from statisticalTools.simplifiedLikelihoods import Data,UpperLimitComputer

def computeULs(inputFile,outputFile,deltas=0.25):

    # ### Load Recast Data
    recastData = pd.read_pickle(inputFile)

    # First compute ULs for each individual SR
    # ATLAS limits:
    # HighPT: nobs = nbg = 9.9, nbgErr = 6.5, lumi = 3000/fb
    # Trackless: nobs = nbg = 17.9, nbgErr = 10.8, lumi = 3000/fb
    atlasUL = {'HighPT' : {'S95_obs' : 11.4, 'S95_exp' : 11.4},
               'Trackless' : {'S95_obs' : 18.0, 'S95_exp' : 18.0}}
    lumi = 3000.0
    # Set r-value
    robs = []
    rexp = []
    # Set r-value
    for _,row in recastData.iterrows():
        S95obs = atlasUL[row['SR']]['S95_obs']
        S95exp = atlasUL[row['SR']]['S95_exp']
        robs.append((row['$N_s$']*(lumi/row['Luminosity (1/fb)'])/S95obs,row['$N_s$ Err']*(lumi/row['Luminosity (1/fb)'])/S95obs))
        rexp.append((row['$N_s$']*(lumi/row['Luminosity (1/fb)'])/S95exp,row['$N_s$ Err']*(lumi/row['Luminosity (1/fb)'])/S95exp))


    robs = np.array(robs)
    rexp = np.array(rexp)
    recastData['robs'] = robs[:,0]
    recastData['rexp'] = rexp[:,0]
    recastData['robsErr'] = robs[:,1]
    recastData['rexpErr'] = rexp[:,1]


    # Now compute combined UL assuming uncorrelated SRs
    # ATLAS data
    atlasBG = {'HighPT' : 9.9, 'Trackless' : 17.9}
    atlasBGErr = {'HighPT' : 6.5, 'Trackless' : 10.8}
    atlasObs =  {'HighPT' : 10, 'Trackless' : 18}
    srOrder = ['HighPT', 'Trackless']

    # Define covariance matrix
    covMatrix = np.zeros((2,2))
    for isr,sr in enumerate(srOrder):
        for jsr,_ in enumerate(srOrder):
            if isr != jsr: continue
            covMatrix[isr,jsr] = atlasBGErr[sr]**2
    
    nobs = [atlasObs[sr] for sr in srOrder]
    nbg = [atlasBG[sr] for sr in srOrder]

     # ### Get all model points
    models = []
    mCols = ['mLLP','mLSP','tau_ns']
    for row in recastData[mCols].values:
        m = dict(zip(mCols,row.tolist()))
        if m not in models:
            models.append(m)

    # ### Loop over model points and compute UL on mu
    ulComp = UpperLimitComputer()

    for m in models:

        dfModel = recastData.loc[(recastData[list(m)] == pd.Series(m)).all(axis=1)]
        ns = []
        for sr in srOrder:
            dataset = dfModel[dfModel['SR'] == sr]
            if len(dataset) != 1:
                print('Error finding single SR/model')
                return None
            ns.append(dataset['$N_s$'].tolist()[0])
        ns = np.array(ns)*(lumi/dataset.iloc[0]['Luminosity (1/fb)'])
        data = Data(observed=nobs, backgrounds=nbg, covariance=covMatrix, 
                    nsignal=ns,deltas_rel=deltas)
        dataExp = Data(observed=[int(b) for b in nbg], backgrounds=nbg, covariance=covMatrix, 
                    nsignal=ns,deltas_rel=deltas)    
        try:
            ul = ulComp.getUpperLimitOnMu(data)
        except:
            print('Error computing ul for model:\n',m,'\n')
            ul = None
        try:
            ulExp = ulComp.getUpperLimitOnMu(dataExp)    
        except:
            print('Error computing ulExp for model:\n',m,'\n')
            ulExp = None
        
        if ul is not None:
            recastData.loc[dfModel.index,'robs_comb'] = 1.0/ul
        else:
            recastData.loc[dfModel.index,'robs_comb'] = 0.0
        if ulExp is not None:
            recastData.loc[dfModel.index,'rexp_comb'] = 1.0/ulExp
        else:
            recastData.loc[dfModel.index,'rexp_comb'] = 0.0
       

    # Store all data to the output file
    recastData.to_pickle(outputFile)

   

if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser(description=
            "Compute the r-values for ATLAS-SUSY-2018-13 for the recast data stored in the input file. \
            The points without any signal will be removed from the output.")
    ap.add_argument('-f', '--inputFile', required=True,
            help='path to the pickle file containing the Pandas DataFrame with the recasting results for the models')
    ap.add_argument('-k', '--kfactor', required=False,default=1.0,type=float,
            help='constant k-factor for all points')    
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

    computeULs(inputFile,outputFile,args.kfactor)
    
    
    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))





