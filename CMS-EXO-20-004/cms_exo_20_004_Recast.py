#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import time
import sys
sys.path.append('../')
from helper import getModelDict,splitModels,getDisplacedJets,getLLPs,getHSCPCandidates
import progressbar as P

delphesDir = os.path.abspath("../DelphesLLP")
os.environ['ROOT_INCLUDE_PATH'] = os.path.join(delphesDir,"external")

import ROOT
import xml.etree.ElementTree as ET


ROOT.gSystem.Load(os.path.join(delphesDir,"libDelphes.so"))

ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')


def passVetoJets2018(jets):
    """
    Calorimeter mitigation failure I for the 2018
    data set.

    :param jets: List of Jet objects

    :return: True if event should be accepted, False otherwise.
    """

    for jet in jets:
        if jet.PT < 30.0:
            continue
        if jet.Phi < -1.57 or jet.Phi > -0.87:
            continue
        if jet.Eta < -3.0 or jet.Eta > -1.3:
            continue
        return False
    return True

def passVetoPtMiss2018(met):
    """
    Calorimeter mitigation failure II for the 2018
    data set.

    :param met: MET object

    :return: True if event should be accepted, False otherwise.
    """        

    if met.MET > 470.:
        return True
    if met.Phi < -1.62 or met.Phi > -0.62:
        return True
    return False

# ### Define dictionary to store data
def getRecastData(inputFiles,llpVeto=False,
                  model='sbottom',modelDict=None,addweights=False):

    if len(inputFiles) > 1:
        print('Combining files:')
        for f in inputFiles:
            print(f)

    if modelDict is None:
        modelDict = getModelDict(inputFiles[0],model)
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

    # ## Cuts
    ## jets
    pTj1min = 100.
    pTjmin = 20.
    etamax = 2.4
    ## MET
    minMET = 250.
    ## Electrons
    pTmin_el = 10.
    etamax_el = 2.5
    nMax_el = 0
    ## Photons
    pTmin_a = 15.
    etamax_a = 2.5
    nMax_a = 0
    ## Muons
    pTmin_mu = 10.
    etamax_mu = 2.4
    nMax_mu = 0
    ## Tau jets
    nMax_tau = 0
    etatau_max = 2.3
    pTtau_min = 18.0
    ## b jets
    nMax_b = 0
    etab_max = 2.4
    pTb_min = 20.0

    luminosities = {2016 : 36.0, 2017 : 41.5, 2018: 59.7}
    lumTot = sum(luminosities.values())
    yields = {ds : [] for ds in luminosities}
    metAll = {ds : [] for ds in luminosities}
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
        # If addweights = Fakse: 
        # assume multiple files correspond to equivalent samplings
        # of the same distributions
        # If addweights = True: directly add events
        if not addweights:
            norm =nevtsDict[inputFile]/modelDict['Total MC Events']
        else:
            norm = 1.0

        for ievt in range(nevts):
            ntotal += 1
            progressbar.update(ntotal)
            tree.GetEntry(ievt)        

            jets = tree.Jet
            try:
                weightPB = tree.Event.At(0).Weight/nevts
            except:
                weightPB = tree.Weight.At(0).Weight
            weightPB = weightPB*norm
            ns = weightPB*1e3*lumTot # number of signal events
            totalweightPB += weightPB

            missingET = tree.MissingET.At(0)
            electrons = tree.Electron
            muons = tree.Muon
            photons = tree.Photon

            # Filter electrons:
            electronList = []
            for iel in range(electrons.GetEntries()):
                electron = electrons.At(iel)
                if electron.PT < pTmin_el:
                    continue
                if abs(electron.Eta) > etamax_el:
                    continue
                electronList.append(electron)

            # Filter muons:
            muonList = []
            for imu in range(muons.GetEntries()):
                muon = muons.At(imu)
                if muon.PT < pTmin_mu:
                    continue
                if abs(muon.Eta) > etamax_mu:
                    continue
                muonList.append(muon)

            # Filter photons:
            photonList = []
            for ia in range(photons.GetEntries()):
                photon = photons.At(ia)
                if photon.PT < pTmin_a:
                    continue
                if abs(photon.Eta) > etamax_a:
                    continue
                photonList.append(photon)            

            # Filter jets
            jetList = []
            bjetList = []
            taujetList = []
            for ijet in range(jets.GetEntries()):
                jet = jets.At(ijet)
                if jet.BTag and jet.PT > pTb_min and abs(jet.Eta) < etab_max:
                    bjetList.append(jet)
                elif jet.TauTag and jet.PT > pTtau_min and abs(jet.Eta) < etatau_max:
                    taujetList.append(jet)
                elif jet.PT > pTjmin and abs(jet.Eta) < etamax:
                    jetList.append(jet)  
            jetList = sorted(jetList, key = lambda j: j.PT, reverse=True)    

            if len(jetList) > 0:
                deltaPhi = np.abs(jetList[0].Phi-missingET.Phi) 
            else:
                deltaPhi = 0.0
            

            # Split event into datasets:
            lumRnd = np.random.uniform(0.,lumTot)
            if lumRnd < luminosities[2016]:
                useDataSet = 2016
            elif lumRnd < luminosities[2016]+luminosities[2017]:
                useDataSet = 2017
            else:
                useDataSet = 2018

            # Apply cuts:
            ## Apply trigger efficiency
            # ns = ns*triggerEff
            if missingET.MET < 120.0:
                continue

            ## Cut on MET
            if missingET.MET < minMET: continue              
            ## Veto electrons
            if len(electronList) > nMax_el: continue  
            ## Veto muons
            if len(muonList) > nMax_mu: continue  
            ## Veto tau jets
            if len(taujetList) > nMax_tau: continue  
            ## Veto b jets
            if len(bjetList) > nMax_b: continue  
            ## Veto photons
            if len(photonList) > nMax_a: continue  
            ## Delta Phi cut
            if deltaPhi < 0.5: continue
            ## Jet cuts
            if len(jetList) < 1 or jetList[0].PT < pTj1min: continue
            if abs(jetList[0].Eta) > etamax: continue

            if llpVeto:            
                llps = getLLPs(tree.bsm,tree.bsmDirectDaughters,
                               tree.bsmFinalDaughters,mothers=tree.bsmMothers)
                hscps = getHSCPCandidates(llps) # Select charged LLPs
                # Veto HSCPs decaying after 1m
                if any(hscp.r_decay > 1e3 for hscp in hscps):
                    continue
                # Veto displaced jets
                jetsDisp = getDisplacedJets(jetList,llps)
                maxR = max([0.0]+[j.llp.r_decay for j in jetsDisp])
                # Veto displacements larger than 2mm
                if maxR > 2.0:
                    continue
            

            if useDataSet == 2018 and not passVetoJets2018(jetList):
                continue            
            if useDataSet == 2018 and not passVetoPtMiss2018(missingET):
                continue
            
            # Store relevant data        
            yields[useDataSet].append(ns)
            metAll[useDataSet].append(missingET.MET)  
            
        f.Close()
    progressbar.finish()

    modelDict['Total xsec-pT150 (pb)'] = 0.0
    # Store total (combined xsec)
    modelDict['Total xsec (pb)'] = totalweightPB
    print('\nCross-section (pb) = %1.3e\n' %totalweightPB)

        

    metBins = [250,  280,  310,  340,  370,  400,  430,  470,  510, 550,  590,  640,  690,  
            740,  790,  840,  900,  960, 1020, 1090, 1160, 1250, 99999]

    # Create a dictionary with lists for each datasets
    dataDict = {'Data-takingperiod' : [2016,2017,2018]}
    dataDict.update({'Luminosity (1/fb)' : [luminosities[ds] 
                      for ds in dataDict['Data-takingperiod']]})
    
    for ibin,b in enumerate(metBins[:-1]):
        label = 'bin_%1.1f_%1.1f'%(b,min(1400.0,metBins[ibin+1]))
        dataDict[label] = []
        dataDict[label+'_ErrorPlus'] = []
        dataDict[label+'_ErrorMinus'] = []

    # Split results into 3 data taking periods:
    dataDict.update({key : [] for key in modelDict})

    for ds in dataDict['Data-takingperiod']:
        # Store common values to all datasets:
        for key,val in modelDict.items():
            dataDict[key].append(val)
        met = metAll[ds]
        ns = np.array(yields[ds])
        binc,binEdges = np.histogram(met,bins=metBins, weights=ns)
        binc2,_ = np.histogram(met,bins=metBins, weights=ns**2)
        for ibin,b in enumerate(binc):
            label = 'bin_%1.1f_%1.1f'%(binEdges[ibin],min(1400.0,binEdges[ibin+1]))    
            dataDict[label].append(b)
            dataDict[label+'_ErrorPlus'].append(np.sqrt(binc2[ibin]))
            dataDict[label+'_ErrorMinus'].append(np.sqrt(binc2[ibin]))

    return dataDict


if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Run the recasting for CMS-EXO-20-004 using one or multiple Delphes ROOT files as input. "
            + "If multiple files are given as argument, add them (the samples weights will be normalized if -n is given)."
            + " Store the signal yields for the SR bins in a pickle (Pandas DataFrame) file." )
    ap.add_argument('-f', '--inputFile', required=True,nargs='+',
            help='path to the ROOT event file(s) generated by Delphes.', default =[])
    ap.add_argument('-o', '--outputFile', required=False,
            help='path to output file storing the DataFrame with the recasting data.'
                 + 'If not defined, will use the name of the first input file', 
            default = None)
    ap.add_argument('-A', '--add', required=False,action='store_true',default=False,
            help='If set, the input files will be considered to refer to samples of the orthogonal processes and their weights will be added.')    
    ap.add_argument('-pt', '--pTcut', required=False,default=0.0,type=float,
            help='Gen level MET cut for computing partial cross-sections.')
    ap.add_argument('-m', '--model', required=False,type=str,default='sbottom',
            help='Defines which model should be considered for extracting model parameters (strong,ewk,gluino,sbottom).')
    ap.add_argument('-llpveto', '--llpVeto', required=False,action='store_true',default=False,
            help='If set, applies a veto on displaced jets matched to a LLP and HSPCs.')

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

    np.random.seed(15)
    t0 = time.time()

    # # Set output file
    args = ap.parse_args()
    inputFiles = args.inputFile
    outputFile = args.outputFile
    
    # Split input files by distinct models and get recast data for
    # the set of files from the same model:
    for fileList,mDict in splitModels(inputFiles,args.model):


        if outputFile is None:
            if not args.llpVeto:
                outFile = fileList[0].replace('delphes_events.root','cms_exo_20_004.pcl')
            else:
                outFile = fileList[0].replace('delphes_events.root','cms_exo_20_004_llpVeto.pcl')
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

        dataDict = getRecastData(fileList,args.llpVeto,
                                 args.model,mDict,addweights=args.add)
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
