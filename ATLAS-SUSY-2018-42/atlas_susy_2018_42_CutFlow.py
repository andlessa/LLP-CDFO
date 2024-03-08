#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import glob
import time
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
def getCutFlow(inputFiles,model='sbottom',modelDict=None,normalize=True):

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
    else:
        massWindows = massLong

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
    totalweightPB = 0.0
    # Keep track of yields for each dataset
    keys = ["Total","Trigger","$E_{T}^{miss}>170$ GeV","$p_{T} > 50$ GeV","Track isolation","$p_{T} > 120$ GeV","$|\eta|<1.8$","$m_{T}({track},{p}_{{T}}^{{ miss}}) > 130$ GeV","(Acceptance)","(SR-Low - no mass Window)","(SR-High - no mass Window)"]
    cutFlow = {k  : np.zeros(2) for k in keys}    


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

            cutFlow["Total"] += (ns,ns**2)

            hscpsFilter = applyHSCPSelection(hscpCandidates,pT=120.,eta=1.8,r=500.0)
            if hscpsFilter:
                cutFlow['(Acceptance)'] += (ns,ns**2)

            muonsLLP = applyMuonTagging(hscpCandidates)
            hscps = [hscp for hscp in hscpCandidates if hscp not in muonsLLP]
            newMETv = removeFromMET(muonsLLP,tree.MissingET.At(0))
            newMET = np.sqrt(newMETv[0]**2+newMETv[1]**2)
            triggerEff = getTriggerEff(metCalo)
            
            ns = ns*triggerEff
            if not ns: continue
            cutFlow['Trigger'] += (ns,ns**2)
            # Apply event selection efficiency
            eventEff = getSelectionEff(newMET)
            ns = ns*eventEff
            if not ns: continue
            cutFlow["$E_{T}^{miss}>170$ GeV"] += (ns,ns**2)

            # Apply selection to HSCP candidates (following ATLAS snippet)
            hscps = applyHSCPSelection(hscps,pT=50.,eta=3.0,r=500.0)
            if not hscps: continue
            cutFlow['$p_{T} > 50$ GeV'] += (ns,ns**2)        
            # hscps = applyIsolation(hscps,tree.Track) # Already included in trackEff
            if not hscps: continue
            cutFlow['Track isolation'] += (ns,ns**2)     
            hscps = applyHSCPSelection(hscps,pT=120.)
            if not hscps: continue
            cutFlow['$p_{T} > 120$ GeV'] += (ns,ns**2)  
            if not hscps: continue
            hscps = applyHSCPSelection(hscps,pT=120.,eta=1.8)
            if not hscps: continue
            cutFlow['$|\eta|<1.8$'] += (ns,ns**2)
            # hscps = applyMTCut(hscps,newMETv) # Already included in trackEff
            if not hscps: continue
            cutFlow['$m_{T}({track},{p}_{{T}}^{{ miss}}) > 130$ GeV'] += (ns,ns**2)

            gbetas = [h.gbeta for h in hscps]
            trackEffHigh = getTrackEff(gbetas,sr='High')
            trackEffLow =  getTrackEff(gbetas,sr='Low')
            nsHigh = ns*(1-np.prod(1.0-trackEffHigh))
            nsLow = ns*(1-np.prod(1.0-trackEffLow))
            cutFlow['(SR-High - no mass Window)'] += (nsHigh,nsHigh**2)
            cutFlow['(SR-Low - no mass Window)'] += (nsLow,nsLow**2)
            
        f.Close()
    progressbar.finish()

    modelDict['Total xsec (pb)'] = totalweightPB
    print('\nCross-section (pb) = %1.3e\n' %totalweightPB)

    cutFlowErr = {k : np.sqrt(v[1]) for k,v in cutFlow.items()}
    cutFlow = {k : v[0]  for k,v in cutFlow.items()}
    
    # Compute normalized cutflow
    print('Cutflow:')
    for key,val in cutFlow.items():
        if key == 'Total':
            continue
        valNorm = float('%1.3e' %(val/cutFlow['Total']))
        errNorm = float('%1.3e' %(cutFlowErr[key]/cutFlow['Total']))
        cutFlow[key] = valNorm
        cutFlowErr[key] = errNorm
    cutFlow['Total'] = 1.0
    cutFlowErr['Total'] = 0.0

    for k,v in cutFlow.items():
        if v != 0.0:
            print('%s : %1.3e +- %1.1f%%' %(k,v,1e2*cutFlowErr[k]/v))
        else:
            print('%s : %1.3e +- ??' %(k,v))

    
    return modelDict,cutFlow,cutFlowErr


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

    # Split input files by distinct models and get recast data for
    # the set of files from the same model:
    for fileList,mDict in splitModels(inputFiles,args.model):        
        modelDict,cutFlow,cutFlowErr = getCutFlow(fileList,args.model,mDict,normalize=args.normalize)
        dataDict = {key : [val] for key,val in modelDict.items()}
        for key,val in cutFlow.items():
            dataDict[key] = [(val,cutFlowErr[key])]

        if outputFile is None:
            outFile = fileList[0].replace('delphes_events.root','atlas_2018_42_cutflow.pcl')
        else:
            outFile = outputFile[:]

        if os.path.splitext(outFile)[1] != '.pcl':
            outFile = os.path.splitext(outFile)[0] + '.pcl'

        # #### Create pandas DataFrame
        df = pd.DataFrame.from_dict(dataDict)

        # ### Save DataFrame to pickle file
        print('Saving to',outFile)
        df.to_pickle(outFile)
        print('\n\n')

    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
