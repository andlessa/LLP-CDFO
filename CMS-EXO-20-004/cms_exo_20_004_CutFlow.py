#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import time
import sys
sys.path.append('../')
from helper import getModelDict,splitModels,getDisplacedJets,getLLPs,getHSCPCandidates
from cms_exo_20_004_Recast import passVetoJets2018,passVetoPtMiss2018
import progressbar as P

delphesDir = os.path.abspath("../DelphesLLP")
os.environ['ROOT_INCLUDE_PATH'] = os.path.join(delphesDir,"external")

import ROOT
import xml.etree.ElementTree as ET


ROOT.gSystem.Load(os.path.join(delphesDir,"libDelphes.so"))

ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')


# ### Define dictionary to store data
def getCutFlow(inputFiles,llpVeto=False,model='sbottom',modelDict=None,addweights=False):

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

    lumi =  59.7
    totalweightPB = 0.0
    keys = ['Total','Triggeremulation','$MET > 250$ GeV', 'Electronveto','Muonveto', 'Tauveto', 'Bjetveto', 'Photonveto','$\Delta \phi (jet,p_{T}^{miss})>0.5$ rad','LeadingAK4jet$p_{T}>100$GeV', 'LeadingAK4jet$\eta<2.4$']
    if llpVeto:
        keys += ['LLP veto']
    keys += ['HCALmitigation(jets)','HCALmitigation($\phi^{miss}$)']
    cutFlow = { k : np.zeros(2) for k in keys}


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
            weightPB = tree.Event.At(0).Weight/nevts
            weightPB = weightPB*norm
            ns = weightPB*1e3*lumi # number of signal events
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
                        
            cutFlow['Total'] += (ns,ns**2)

            # Apply cuts:
            ## Apply trigger efficiency
            # ns = ns*triggerEff
            if missingET.MET < 120.0:
                continue
            cutFlow['Triggeremulation'] += (ns,ns**2)

            ## Cut on MET
            if missingET.MET < minMET: continue              
            cutFlow['$MET > 250$ GeV'] += (ns,ns**2)
            ## Veto electrons
            if len(electronList) > nMax_el: continue  
            cutFlow['Electronveto'] += (ns,ns**2)
            ## Veto muons
            if len(muonList) > nMax_mu: continue  
            cutFlow['Muonveto'] += (ns,ns**2)
            ## Veto tau jets
            if len(taujetList) > nMax_tau: continue  
            cutFlow['Tauveto'] += (ns,ns**2)
            ## Veto b jets
            if len(bjetList) > nMax_b: continue  
            cutFlow['Bjetveto'] += (ns,ns**2)
            ## Veto photons
            if len(photonList) > nMax_a: continue  
            cutFlow['Photonveto'] += (ns,ns**2)
            ## Delta Phi cut
            if deltaPhi < 0.5: continue
            cutFlow['$\Delta \phi (jet,p_{T}^{miss})>0.5$ rad'] += (ns,ns**2)
            ## Jet cuts
            if len(jetList) < 1 or jetList[0].PT < pTj1min: continue
            cutFlow['LeadingAK4jet$p_{T}>100$GeV'] += (ns,ns**2)
            if abs(jetList[0].Eta) > etamax: continue
            cutFlow['LeadingAK4jet$\eta<2.4$'] += (ns,ns**2)

            if llpVeto:            
                llps = getLLPs(tree.bsm,tree.bsmDirectDaughters,tree.bsmFinalDaughters)
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
            
                cutFlow['LLP veto'] += (ns,ns**2)


            if not passVetoJets2018(jetList):
                continue
            cutFlow['HCALmitigation(jets)'] += (ns,ns**2)
            if not passVetoPtMiss2018(missingET):
                continue
            cutFlow['HCALmitigation($\phi^{miss}$)'] += (ns,ns**2)

        f.Close()
    progressbar.finish()


    # Store total (combined xsec)
    modelDict['Total xsec (pb)'] = totalweightPB
    # print('\nCross-section (pb) = %1.3e\n' %totalweightPB)

    cutFlowErr = {k : np.sqrt(v[1]) for k,v in cutFlow.items()}
    cutFlow = {k : v[0]  for k,v in cutFlow.items()}
    
    print('-'*10)
    print('Model:')
    for key,val in modelDict.items():
        print("%s = %1.5e" %(key,val))

    # Compute normalized cutflow
    print('-'*10)
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
            "Run the recasting for CMS-EXO-20-004 using one or multiple Delphes ROOT files as input. "
            + "If multiple files are given as argument, add them (the samples weights will be normalized if -n is given)."
            + " Store the cutflow and SR bins in a pickle (Pandas DataFrame) file." )
    ap.add_argument('-f', '--inputFile', required=True,nargs='+',
            help='path to the ROOT event file(s) generated by Delphes.', default =[])
    ap.add_argument('-o', '--outputFile', required=False,
            help='path to output file storing the DataFrame with the recasting data.'
                 + 'If not defined, will use the name of the first input file', 
            default = None)
    ap.add_argument('-A', '--add', required=False,action='store_true',default=False,
            help='If set, the input files will be considered to refer to samples of the orthogonal processes and their weights will be added.')    
    ap.add_argument('-m', '--model', required=False,type=str,default='sbottom',
            help='Defines which model should be considered for extracting model parameters (strong,ewk,gluino,sbottom).')
    ap.add_argument('-llpveto', '--llpVeto', required=False,action='store_true',default=False,
            help='If set, applies a veto on displaced jets matched to a LLP and HSPCs.')



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

    np.random.seed(15)
    t0 = time.time()

    # # Set output file
    args = ap.parse_args()
    inputFiles = args.inputFile
    outputFile = args.outputFile
    
    # Split input files by distinct models and get recast data for
    # the set of files from the same model:
    for fileList,mDict in splitModels(inputFiles,args.model):
        modelDict,cutFlow,cutFlowErr = getCutFlow(fileList,args.llpVeto,
                                                  args.model,mDict,addweights=args.add)
        dataDict = {key : [val] for key,val in modelDict.items()}
        for key,val in cutFlow.items():
            dataDict[key] = [(val,cutFlowErr[key])]

        if outputFile is None:
            if args.maxJetR < 0.0:
                outFile = fileList[0].replace('delphes_events.root','cms_exo_20_004_cutflow.pcl')
            else:
                outFile = fileList[0].replace('delphes_events.root','cms_exo_20_004_cutflow_llpVeto.pcl')
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
