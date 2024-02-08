#!/usr/bin/env python3

import os,glob
from typing import Any
import numpy as np
import pandas as pd
import glob
import pyslha
import time
import progressbar as P
import sys
from ATLAS_data.effFunctions import eventEff,vertexEff
from atlas_susy_2016_08_Recast import (getLLPs, getJets, eventAcc, 
                                       vertexAcc, getModelDict)


delphesDir = os.path.abspath("../DelphesLLP")
os.environ['ROOT_INCLUDE_PATH'] = os.path.join(delphesDir,"external")

import ROOT
import xml.etree.ElementTree as ET


ROOT.gSystem.Load(os.path.join(delphesDir,"libDelphes.so"))

ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')

# ### Define dictionary to store data
def getCutFlow(inputFiles,model='bb'):

    if len(inputFiles) > 1:
        print('Combining files:')
        for f in inputFiles:
            print(f)

    modelDict = getModelDict(inputFiles,model)
    if not modelDict:
        modelDict = {}

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


    lumi = 32.8
    totalweightPB = 0.0
    # Keep track of yields for each dataset
    cutFlow = { "Total" : 0.0,
                "Jet Selection" : 0.0,
                "MET > 200" : 0.0,
                "$R_{xy},z <$ 300 mm" : 0.0,
                "$R_{DV} > 4$ mm" : 0.0,
                "$d_0 > 2$ mm" : 0.0,
                "$nTracks >= 5$" : 0.0,
                "mDV > 10 GeV" : 0.0,
                "+Evt Eff" : 0.0,
                "+DV Eff" : 0.0
                }
    

    progressbar = P.ProgressBar(widgets=["Reading %i Events: " %modelDict['Total MC Events'], 
                                P.Percentage(),P.Bar(marker=P.RotatingMarker()), P.ETA()])
    progressbar.maxval = modelDict['Total MC Events']
    progressbar.start()

    ntotal = 0
    totalweightPB = 0.0
    for inputFile in inputFiles:
        f = ROOT.TFile(inputFile,'read')
        tree = f.Get("Delphes")
        nevts = tree.GetEntries()
        norm =nevtsDict[inputFile]/modelDict['Total MC Events']

        for ievt in range(nevts):    
            
            ntotal += 1
            progressbar.update(ntotal)
            tree.GetEntry(ievt)   
            weightPB = tree.Weight.At(1).Weight     
            weightPB = weightPB*norm
            totalweightPB += weightPB
            ns = weightPB*1e3*lumi # number of signal events

            jets = getJets(tree.GenJet,pTmin=25.,etaMax=5.0)
            met = tree.GenMissingET.At(0).MET

            cutFlow["Total"] += ns

            # Event acceptance
            evt_acc = eventAcc(jets,met,metCut=0.0,
                               maxJetChargedPT=5.0,minJetPt1=70.,
                               minJetPt2=25.,minPVdistance=4.0)
            
            if (not evt_acc):
                continue

            ns = ns*evt_acc
            cutFlow["Jet Selection"] += ns

            if met < 200.0:
                continue
            
            cutFlow["MET > 200"] += ns

            llps = getLLPs(tree.bsm,tree.bsmDirectDaughters,tree.bsmFinalDaughters)


            llpsSel = [llp for llp in llps if vertexAcc(llp,Rmax=300.0,zmax=300.0)]
            if not llpsSel: continue
            cutFlow["$R_{xy},z <$ 300 mm"] += ns

            llpsSel = [llp for llp in llps if vertexAcc(llp,Rmax=300.0,zmax=300.0,Rmin=4.0)]
            if not llpsSel: continue
            cutFlow["$R_{DV} > 4$ mm"] += ns

            llpsSel = [llp for llp in llps if vertexAcc(llp,Rmax=300.0,zmax=300.0,Rmin=4.0,d0min=2.0)]
            if not llpsSel: continue
            cutFlow["$d_0 > 2$ mm"] += ns

            llpsSel = [llp for llp in llps if vertexAcc(llp,Rmax=300.0,zmax=300.0,Rmin=4.0,d0min=2.0,nmin=5)]
            if not llpsSel: continue
            cutFlow["$nTracks >= 5$"] += ns

            llpsSel = [llp for llp in llps if vertexAcc(llp,Rmax=300.0,zmax=300.0,Rmin=4.0,d0min=2.0,nmin=5,mDVmin=10.0)]
            if not llpsSel: continue
            cutFlow["mDV > 10 GeV"] += ns
           
            # Event efficiency
            evt_eff = eventEff(met,llpsSel)
            # Vertex acceptances:
            v_acc = np.array([vertexAcc(llp,Rmax=300.0,zmax=300.0,Rmin=4.0,d0min=2.0,nmin=5,mDVmin=10.0) for llp in llps])
            # Vertex efficiencies:
            v_eff = np.array([vertexEff(llp) for llp in llps])

            ns = ns*evt_eff

            cutFlow["+Evt Eff"] += ns
            
            
            wvertex = 1.0-np.prod(1.0-v_acc*v_eff)
            
            # Add to the total weight in each SR:
            cutFlow["+DV Eff"] += ns*wvertex

        f.Close()
    progressbar.finish()

    modelDict['Total xsec (pb)'] = totalweightPB
    print('\nCross-section (pb) = %1.3e\n' %totalweightPB)

    # Compute normalized cutflow
    print('Cutflow:')
    for key,val in cutFlow.items():
        if key == 'Total':
            continue
        valNorm = float('%1.3e' %(val/cutFlow['Total']))
        cutFlow[key] = valNorm
    cutFlow['Total'] = 1.0

    for k,v in cutFlow.items():
        print('%s : %1.3e' %(k,v))


    return cutFlow


if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Run the recasting for ATLAS-SUSY-2018-13 using one or multiple Delphes ROOT files as input. "
            + "If multiple files are given as argument, add them (the samples weights will be normalized if -n is given)."
            + " Store the cutflow and SR bins in a pickle (Pandas DataFrame) file." )
    ap.add_argument('-f', '--inputFile', required=True,nargs='+',
            help='path to the ROOT event file(s) generated by Delphes.', default =[])
    ap.add_argument('-o', '--outputFile', required=False,
            help='path to output file storing the DataFrame with the recasting data. '
                 + 'If not defined, will use the name of the first input file', 
            default = None)
    ap.add_argument('-m', '--model', required=False,type=str,default='sbottom',
            help='Defines which model should be considered for extracting model parameters (strong,ewk,gluino,sbottom).')


    # First make sure the correct env variables have been set:
    import subprocess
    import sys
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
    if outputFile is None:
        outputFile = inputFiles[0].replace('delphes_events.root','atlas_2016_08_cutflow.pcl')

    if os.path.splitext(outputFile)[1] != '.pcl':
        outputFile = os.path.splitext(outputFile)[0] + '.pcl'

    cutFlow = getCutFlow(inputFiles,args.model)
    for key,val in cutFlow.items():
        cutFlow[key] = [val]

    # #### Create pandas DataFrame
    df = pd.DataFrame.from_dict(cutFlow)

    # ### Save DataFrame to pickle file
    print('Saving to',outputFile)
    df.to_pickle(outputFile)

    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))