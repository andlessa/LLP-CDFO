#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import time
import sys
sys.path.append('../')
import progressbar as P
from helper import getLLPs,getModelDict,splitModels
from ATLAS_data.effFunctions import (getMuonRecoEff,getTriggerEff,getTrackEff,
                                     getSelectionEff,getTargetMass,getMassSelEff,
                                     massLong,massShort)

delphesDir = os.path.abspath("../DelphesLLP")
os.environ['ROOT_INCLUDE_PATH'] = os.path.join(delphesDir,"external")

import ROOT
import xml.etree.ElementTree as ET


ROOT.gSystem.Load(os.path.join(delphesDir,"libDelphes.so"))

ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')


def getHSCPCandidates(rhadrons,rhadronDaughters,llps):

    candidates = []    
    for ip in range(rhadrons.GetEntries()):
        rhadron = rhadrons.At(ip)
        if abs(rhadron.Charge) != 1: # Skip neutral or multicharged particles
            continue
        
        # Map rhadron to llp (parton)
        # (check which llp came from the R-hadron)
        hscpCandidate = None
        for jp in range(rhadronDaughters.GetEntries()):
            d = rhadronDaughters.At(jp)
            if d.M1 != ip:
                continue # Daughter does not belong to current rhadron
            for llp in llps:
                if d is llp._candidate:
                    hscpCandidate = llp
                    break
            if hscpCandidate is not None:
                break
        candidates.append(hscpCandidate)

    return candidates

def applyHSCPSelection(hscpList,pT=50.,eta=2.4,r=500.):

    selHSCPs = []
    for hscp in hscpList:
        if hscp.PT < pT: continue 
        if abs(hscp.Eta) > eta: continue
        if hscp.r_decay < r: continue
        selHSCPs.append(hscp)
    
    return selHSCPs

def applyIsolation(hscpList,pTmax=5.0):

    isoHSCPs = []
    # Apply isolation requirement for HSCP tracks
    for hscp in hscpList:
        sumPT = hscp.SumPtCharged
        if sumPT > pTmax: continue
        isoHSCPs.append(hscp)
    return isoHSCPs

def applyMuonTagging(hscpList):

    """
    Computes the probability of reconstructing the hscp as a muon.
    :param hscpList: List of GenParticle objects
    """

    muonsLLP = []
    for hscp in hscpList:
        if hscp.r_decay < 3.9e3 and hscp.z_decay < 6e3: # Skip decays before MS
            continue
        
        beta = hscp.beta
        eta = abs(hscp.Eta)
        eff = getMuonRecoEff(beta,eta,hscp.PID)
        # Randomly reconstrunct the HSCP as a muon
        if np.random.uniform() < eff:
            continue
        muonsLLP.append(hscp)
    
    return muonsLLP
def removeFromMET(particles,METobj):
    """
    Removes the contribution from the particles in the list
    to the total MET.
    """

    metx = METobj.MET*np.cos(METobj.Phi)
    mety = METobj.MET*np.sin(METobj.Phi)

    if particles:
        # Remove particles from MET:            
        pxTot = sum([p.Px for p in particles])
        pyTot = sum([p.Py for p in particles])        
        metx = (metx-pxTot)
        mety = (mety-pyTot)
    
    return [metx,mety]

def applyMTCut(hscps,METvector):
    """
    Remove tracks which have mT < 130 GeV
    """

    selHSCPs = []
    met = np.sqrt(METvector[0]**2 + METvector[1]**2)    
    for hscp in hscps:
        pThscp = [hscp.Px,hscp.Py]
        cosdphi = np.dot(pThscp,METvector)/(hscp.PT*met)
        mT = np.sqrt(2*hscp.PT*met*(1-cosdphi))
        if mT < 130.: continue
        selHSCPs.append(hscp)
    
    return selHSCPs

# ### Define dictionary to store data
def getRecastData(inputFiles,model='sbottom',modelDict=None,normalize=True):

    if len(inputFiles) > 1:
        print('Combining files:')
        for f in inputFiles:
            print(f)

    if modelDict is None:
        modelDict = getModelDict(inputFiles[0],model)
    if not modelDict:
        modelDict = {}

    # Select mass windows accoding to lifetime
    # (use long lifetime by default)
    if 'width' not in modelDict or (modelDict['width'] > 0 and (6.582e-25/modelDict['width']) < 1e-9):
        massWindows = massShort
        lifetimeRegime = 'Short Lifetime'
    else:
        massWindows = massLong
        lifetimeRegime = 'Long Lifetime'

    massWindows = massWindows[['Mass_window_Low','Mass_window_High','Target_Mass_GeV']]

    modelDict['Total MC Events'] = 0

    nevtsDict = {}
    # Get total number of events:
    for inputFile in inputFiles:
        f = ROOT.TFile(inputFile,'read')
        tree = f.Get("Delphes")
        nevts = tree.GetEntries()
        modelDict['Total MC Events'] += nevts        
        nevtsDict[inputFile] = nevts
        f.Close()


    lumi = 139.0
    yields = {}
    totalweightPB = 0.0


    progressbar = P.ProgressBar(widgets=["Reading %i Events: " %modelDict['Total MC Events'], 
                                P.Percentage(),P.Bar(marker=P.RotatingMarker()), P.ETA()])
    progressbar.maxval = modelDict['Total MC Events']
    progressbar.start()

    ntotal = 0
    for inputFile in inputFiles:
        f = ROOT.TFile(inputFile,'read')
        tree = f.Get("Delphes")
        nevts = tree.GetEntries()
        # If normalize = True: 
        # assume multiple files correspond to equivalent samplings
        # of the same distributions
        # If normalize = False: directly add events
        if normalize:
            norm =nevtsDict[inputFile]/modelDict['Total MC Events']
        else:
            norm = 1.0

        for ievt in range(nevts):    
            
            ntotal += 1
            progressbar.update(ntotal)
            tree.GetEntry(ievt)   
            weightPB = tree.Weight.At(1).Weight     
            weightPB = weightPB*norm
            totalweightPB += weightPB
            ns = weightPB*1e3*lumi # number of signal events

            metCalo = tree.MissingETCalo.At(0).MET
            llps = getLLPs(tree.bsm,tree.bsmDirectDaughters,tree.bsmFinalDaughters)
            hscpCandidates = getHSCPCandidates(tree.isoRhadrons,tree.rhadronDaughters,llps)

            muonsLLP = applyMuonTagging(hscpCandidates)
            hscps = [hscp for hscp in hscpCandidates if hscp not in muonsLLP]
            newMETv = removeFromMET(muonsLLP,tree.MissingET.At(0))
            newMET = np.sqrt(newMETv[0]**2+newMETv[1]**2)
            triggerEff = getTriggerEff(metCalo)
            
            ns = ns*triggerEff
            if not ns: continue
            # Apply event selection efficiency
            eventEff = getSelectionEff(newMET)
            ns = ns*eventEff
            if not ns: continue

            # Apply selection to HSCP candidates (following ATLAS snippet)
            hscps = applyHSCPSelection(hscps,pT=50.,eta=3.0,r=500.0)
            if not hscps: continue
            # hscps = applyIsolation(hscps,tree.Track) # Already included in trackEff
            if not hscps: continue
            hscps = applyHSCPSelection(hscps,pT=120.)
            if not hscps: continue
            if not hscps: continue
            hscps = applyHSCPSelection(hscps,pT=120.,eta=1.8)
            if not hscps: continue
            # hscps = applyMTCut(hscps,newMETv) # Already included in trackEff
            if not hscps: continue

            gbetas = [h.gbeta for h in hscps]
            trackEffHigh = getTrackEff(gbetas,sr='High')
            trackEffLow =  getTrackEff(gbetas,sr='Low')

            # Assume hscp masses can be approximated
            # by the closest target mass for computing the
            # mass window efficiency:
            masses = [getTargetMass(hscp.Mass) for hscp in hscps]
            if all(m is None for m in masses):
                continue            
            # masses = [hscp.Mass for hscp in hscps] # Use target mass or real HSCP mass?

            # Loop over masses 
            # (if HSCP have different masses, each one can contribute to distinct
            # mass windows/target masses)
            for mTarget in set(masses):
                if not mTarget: continue
                # For the HSCP masses not matching the target mass, set them to None
                # so they will have zero mass selection efficiencies and will not be
                # counted to the target mass SR
                selectedMasses = [h.Mass if masses[ih] == mTarget else None 
                                    for ih,h in enumerate(hscps)]
                # Get mass window selection efficiency (zero if mass = None)
                wmassSRHigh = getMassSelEff(selectedMasses,sr='High')
                wmassSRLow = getMassSelEff(selectedMasses,sr='Low')
                
                # Get total yield for respective target mass and 
                # the High and Low SRs:
                yieldHigh = ns*(1-np.prod(1.0-trackEffHigh*wmassSRHigh))            
                yieldLow = ns*(1-np.prod(1.0-trackEffLow*wmassSRLow))
            
                # In case there are distinct masses, use the largest one (correct?)
                # (if there are no good masses, store in 0.)
                if not mTarget in yields:
                    yields[mTarget] = {'SR-Inclusive_Low': [], 'SR-Inclusive_High' : []}
                    
                # Store event for a given target window
                yields[mTarget]['SR-Inclusive_Low'].append(yieldLow)
                yields[mTarget]['SR-Inclusive_High'].append(yieldHigh)
            
        f.Close()
    progressbar.finish()

    modelDict['Total xsec (pb)'] = totalweightPB
    print('\nCross-section (pb) = %1.3e\n' %totalweightPB)
    
    # Create a dictionary for storing data
    dataDict = {'Luminosity (1/fb)' : lumi, 'Regime' : lifetimeRegime}
    # Signal regions
    if not yields:
        dataDict['SR'] = ['SR-Inclusive_Low', 'SR-Inclusive_High']
        dataDict['Target Mass [GeV]'] = [0.0,0.0]
        dataDict['$N_s$'] = [0.0,0.0]
        dataDict['$\sigma_{Ns}$'] = [0.0,0.0]
    else:
        dataDict['SR'] = []
        dataDict['Target Mass [GeV]'] = []
        dataDict['$N_s$'] = []
        dataDict['$\sigma_{Ns}$'] = []

        # Total signal regions = mass targets * (Low, High)
        for targetMass in yields.keys():
            for sr in yields[targetMass]:
                ns = np.array(yields[targetMass][sr])
                nsTot = sum(ns)
                nsError = np.sqrt(sum(ns**2))
                if not nsTot:
                    continue # Skip empty bins
                dataDict['SR'].append(sr)
                dataDict['Target Mass [GeV]'].append(targetMass)
                dataDict['$N_s$'].append(nsTot)
                dataDict['$\sigma_{Ns}$'].append(nsError)    

    # Expand cutflow and modelDict to match number of rows in dataDict:
    for key,val in modelDict.items():
        modelDict[key] = [val]*len(dataDict['SR'])

    # Create a dictionary for storing data
    dataDict.update(modelDict)
   

    return dataDict


if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Run the recasting for ATLAS-SUSY-2018-42 using one or multiple Delphes ROOT files as input. "
            + "If multiple files are given as argument, add them (the samples weights will be normalized if -n is given)."
            + " Store the cutflow and SR bins in a pickle (Pandas DataFrame) file." )
    ap.add_argument('-f', '--inputFile', required=True,nargs='+',
            help='path to the ROOT event file(s) generated by Delphes.', default =[])
    ap.add_argument('-o', '--outputFile', required=False,
            help='path to output file storing the DataFrame with the recasting data.'
                 + 'If not defined, will use the name of the first input file', 
            default = None)
    ap.add_argument('-n', '--normalize', required=False,action='store_true',
            help='If set, the input files will be considered to refer to multiple samples of the same process and their weights will be normalized.')
    ap.add_argument('-m', '--model', required=False,type=str,default='wino',
            help='Defines which model should be considered for extracting model parameters (stau,wino,gluino).')
    ap.add_argument('-U', '--update', required=False,action='store_true',
            help='If the flag is set only the model points containing data newer than the dataframe will be read.')


    ap.add_argument('-v', '--verbose', default='info',
            help='verbose level (debug, info, warning or error). Default is info')


    # First make sure the correct env variables have been set:
    import subprocess
    import sys
    from datetime import datetime as dt
    LDPATH = subprocess.check_output('echo $LD_LIBRARY_PATH',shell=True,text=True)
    ROOTINC = subprocess.check_output('echo $ROOT_INCLUDE_PATH',shell=True,text=True)
    pythiaDir = os.path.abspath('../MG5/HEPTools/pythia8/lib')
    delphesDir = os.path.abspath('../DelphesLLP/external')
    if pythiaDir not in LDPATH or delphesDir not in ROOTINC:
        print('Enviroment variables not properly set. Run source setenv.sh first.')
        sys.exit()

    t0 = time.time()

    # # Set output file
    args = ap.parse_args()
    inputFiles = args.inputFile
    outputFile = args.outputFile

    if args.update:
        print('\n\n======= Updating files with new event data ==========\n')

    # Split input files by distinct models and get recast data for
    # the set of files from the same model:
    for fileList,mDict in splitModels(inputFiles,args.model):

        # Set output file
        if outputFile is None:
            outFile = fileList[0].replace('delphes_events.root','atlas_2018_42.pcl')
        else:
            outFile = outputFile[:]

        if os.path.splitext(outFile)[1] != '.pcl':
            outFile = os.path.splitext(outFile)[0] + '.pcl'

        skipModel = False
        if args.update and os.path.isfile(outFile):
            outFile_date = dt.fromtimestamp(os.path.getctime(outFile))
            inputFiles_date = max([dt.fromtimestamp(os.path.getctime(f)) for f in fileList])
            if inputFiles_date <= outFile_date:
                skipModel = True
        if skipModel:
            print('\nSkipping',mDict,'\n')
            # print('files=',fileList)
            # sys.exit()
            continue

        print('----------------------------------')
        print('\t Model: %s (%i files)' %(mDict,len(fileList)))

        dataDict = getRecastData(fileList,args.model,mDict,normalize=args.normalize)
        if args.verbose == 'debug':
            for k,v in dataDict.items():
                print(k,v)
        

        # #### Create pandas DataFrame
        df = pd.DataFrame.from_dict(dataDict)
        # ### Save DataFrame to pickle file
        print('Saving to',outFile)
        df.to_pickle(outFile)
        print('\n')



    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
