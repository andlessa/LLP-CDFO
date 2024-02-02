#!/usr/bin/env python3

import os,glob
import numpy as np
import pandas as pd
import glob
import pyslha
import time
import progressbar as P

delphesDir = os.path.abspath("../../MG5/Delphes")
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


def getModelDict(inputFiles):

    # cmsColumns = ['Coupling', 'Mode', '$m_{med}$', '$M_{D}$' '$m_{DM}$', '$d$' '$g_{DM}$', '$g_{q}$']
    modelDict = {}
    # for column in cmsColumns:
        # modelDict[column] = None
    
    # ## Get Model Parameters
    banner = sorted(glob.glob(os.path.dirname(inputFiles[0])+'/*banner.txt'),key=os.path.getmtime,reverse=True)
    if len(banner) == 0:
        print('Banner not found for %s' %inputFiles[0])
        return None
    elif len(banner) > 1:        
        print('\n%i banner files found in %s.' 
            %(len(banner),os.path.dirname(inputFiles[0])))
        f = os.path.basename(inputFiles[0])
        matches = []
        for b in banner:
            if os.path.splitext(f)[0] in b:
                banner = b
                break
            matches.append(set(f).intersection(set(os.path.basename(b))))
        if isinstance(banner,list):
            banner = banner[np.argmax(matches)]
        print('Using banner %s'%banner)
    else:
        banner = banner[0]

    xtree = ET.parse(banner)
    xroot = xtree.getroot()
    slha = xroot.find('header').find('slha').text
    pars = pyslha.readSLHA(slha)
    if 55 in pars.blocks['MASS']:
        model = 'spin1'
        mMed = pars.blocks['MASS'][55]
    elif 54 in pars.blocks['MASS']:
        model = 'spin0'
        mMed = pars.blocks['MASS'][54]
    elif 5000039 in pars.blocks['MASS']:
        model = 'ADD'
        mMed = pars.blocks['MASS'][5000039]
    elif 5000002 in pars.blocks['MASS']:
        model = 'stop'
        mMed = pars.blocks['MASS'][5000002]

    if model in ['spin1','spin0']:
        mDM = pars.blocks['MASS'][52]
    elif model in ['stop']:
        mDM = pars.blocks['MASS'][5000012]
    else:
        mDM = None

    if model == 'spin1':
        gVq = pars.blocks['DMINPUTS'][4] # Mediator-quark vector coupling
        gAq = pars.blocks['DMINPUTS'][10] # Mediator-quark axial coupling
        gVx = pars.blocks['DMINPUTS'][2] # Mediator-DM vector coupling
        gAx = pars.blocks['DMINPUTS'][3] # Mediator-DM axial coupling
    elif model == 'spin0':
        gVq = pars.blocks['DMINPUTS'][6] # Mediator-quark scalar coupling
        gAq = pars.blocks['DMINPUTS'][12] # Mediator-quark pseudoscalar coupling
        gVx = pars.blocks['DMINPUTS'][3] # Mediator-DM scalar coupling
        gAx = pars.blocks['DMINPUTS'][4] # Mediator-DM pseudoscalar coupling
    elif model == 'ADD':
        MD = pars.blocks['ADDINPUTS'][1] # Fundamental Planck scale in large extra dimensions
        d = pars.blocks['ADDINPUTS'][2] # Number of extra dimentions
    elif model == 'stop':
        yDM = pars.blocks['FRBLOCK'][1] # Fundamental Planck scale in large extra dimensions
    
    print('\nModel parameters:')
    if model in ['spin1','spin0']:
        print('mMed = %1.2f GeV, mDM = %1.2f GeV, gVq = %1.2f, gAq = %1.2f, gVx = %1.2f, gAx = %1.2f\n' 
            %(mMed,mDM,gVq,gAq,gVx,gAx))
    elif model == 'ADD':
        print('MD = %1.2f GeV, d = %1.2f\n'  %(MD,d))
    elif model == 'stop':
        print('Mstop = %1.2f GeV, mDM = %1.2f, yDM = %1.1f \n'  %(mMed,mDM,yDM))


    # #### Store data
    if model in ['spin1','spin0']:
        if gVx != 0:
            if model == 'spin1':
                modelDict['Coupling'] = 'Vector'
            elif model == 'spin0':
                modelDict['Coupling'] = 'Scalar'
        else:
            if model == 'spin1':
                modelDict['Coupling'] = 'Axial'
            elif model == 'spin0':
                modelDict['Coupling'] = 'Pseudoscalar'
    elif model == 'ADD':
        modelDict['Coupling'] = 'ADD'
    elif model == 'stop':
        modelDict['Coupling'] = 'Stop'


    modelDict['Mode'] = 'DM+QCDjets'

    if model in ['spin1','spin0']:
        modelDict['$m_{med}$'] = mMed
        modelDict['$m_{DM}$'] = mDM
        if (modelDict['Coupling'] == 'Vector') or (modelDict['Coupling'] == 'Scalar'):
            modelDict['$g_{DM}$'] = gVx
            modelDict['$g_{q}$'] = gVq
        else:
            modelDict['$g_{DM}$'] = gAx
            modelDict['$g_{q}$'] = gAq
    elif model == 'ADD':
        modelDict['$M_{D}$'] = MD
        modelDict['$d$'] = d
    elif model == 'stop':
        modelDict['$m_{\tilde t}$'] = mMed
        modelDict['$m_{\tilde \chi_1^0}$'] = mDM
        modelDict['$y_{DM}$'] = yDM

    return modelDict

# ### Define dictionary to store data
def getRecastData(inputFiles,pTcut=150.0,normalize=False):

    if len(inputFiles) > 1:
        print('Combining files:')
        for f in inputFiles:
            print(f)

    modelDict = getModelDict(inputFiles)
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
    # Keep track of yields for each dataset
    cutFlowAll = {ds : 
                    {'Fullsample' : 0.0,
                    'Triggeremulation' : 0.0,
                    '$p_{T}^{miss}>250$GeV' : 0.0, 
                    'Electronveto' : 0.0,
                    'Muonveto' : 0.0, 
                    'Tauveto' : 0.0, 
                    'Bjetveto' : 0.0, 
                    'Photonveto' : 0.0,
                    '$\Delta \phi (jet,p_{T}^{miss})>0.5$ rad' : 0.0,
                    'LeadingAK4jet$p_{T}>100$GeV' : 0.0, 
                    'LeadingAK4jet$\eta<2.4$' : 0.0,            
                    'HCALmitigation(jets)' : 0.0,
                    'HCALmitigation($\phi^{miss}$)' : 0.0}
                for ds in luminosities}


    progressbar = P.ProgressBar(widgets=["Reading %i Events: " %modelDict['Total MC Events'], 
                                P.Percentage(),P.Bar(marker=P.RotatingMarker()), P.ETA()])
    progressbar.maxval = modelDict['Total MC Events']
    progressbar.start()

    ntotal = 0
    for inputFile in inputFiles:
        f = ROOT.TFile(inputFile,'read')
        tree = f.Get("Delphes")
        nevts = tree.GetEntries()
        if normalize:
            norm =nevtsDict[inputFile]/modelDict['Total MC Events']
        else:
            norm = 1.0

        for ievt in range(nevts):    
            
            ntotal += 1
            progressbar.update(ntotal)
            tree.GetEntry(ievt)        

            jets = tree.Jet
            try:
                weightPB = tree.Weight.At(1).Weight
            except:
                weightPB = tree.Weight.At(0).Weight
            weightPB = weightPB*norm
            ns = weightPB*1e3*lumTot # number of signal events
            totalweightPB += weightPB

            missingET = tree.MissingET.At(0)
        #         missingET = tree.GenMissingET.At(0)  # USE REAL MISSING ET!
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
            
            # Apply cut on DM pT to reproduce CMS
            # cut on event generation:
            dmMET = tree.DMMissingET.At(0).MET
            if dmMET < pTcut:
                continue

            # Split event into datasets:
            lumRnd = np.random.uniform(0.,lumTot)
            if lumRnd < luminosities[2016]:
                useDataSet = 2016
            elif lumRnd < luminosities[2016]+luminosities[2017]:
                useDataSet = 2017
            else:
                useDataSet = 2018

            
            cutFlow = cutFlowAll[useDataSet]
            cutFlow['Fullsample'] += ns

            # Apply cuts:
            ## Apply trigger efficiency
            # ns = ns*triggerEff
            if missingET.MET < 120.0:
                continue
            cutFlow['Triggeremulation'] += ns

            ## Cut on MET
            if missingET.MET < minMET: continue              
            cutFlow['$p_{T}^{miss}>250$GeV'] += ns
            ## Veto electrons
            if len(electronList) > nMax_el: continue  
            cutFlow['Electronveto'] += ns
            ## Veto muons
            if len(muonList) > nMax_mu: continue  
            cutFlow['Muonveto'] += ns
            ## Veto tau jets
            if len(taujetList) > nMax_tau: continue  
            cutFlow['Tauveto'] += ns
            ## Veto b jets
            if len(bjetList) > nMax_b: continue  
            cutFlow['Bjetveto'] += ns
            ## Veto photons
            if len(photonList) > nMax_a: continue  
            cutFlow['Photonveto'] += ns
            ## Delta Phi cut
            if deltaPhi < 0.5: continue
            cutFlow['$\Delta \phi (jet,p_{T}^{miss})>0.5$ rad'] += ns
            ## Jet cuts
            if len(jetList) < 1 or jetList[0].PT < pTj1min: continue
            cutFlow['LeadingAK4jet$p_{T}>100$GeV'] += ns
            if abs(jetList[0].Eta) > etamax: continue
            cutFlow['LeadingAK4jet$\eta<2.4$'] += ns
            if useDataSet == 2018 and not passVetoJets2018(jetList):
                continue
            cutFlow['HCALmitigation(jets)'] += ns
            if useDataSet == 2018 and not passVetoPtMiss2018(missingET):
                continue
            cutFlow['HCALmitigation($\phi^{miss}$)'] += ns

            # Store relevant data        
            yields[useDataSet].append(ns)
            metAll[useDataSet].append(missingET.MET)  
            
        f.Close()
    progressbar.finish()

    modelDict['Total xsec-pT150 (pb)'] = 0.0
    # Store total (combined xsec)
    modelDict['Total xsec (pb)'] = totalweightPB
    print('\nCross-section (pb) = %1.3e\n' %totalweightPB)

    for cutFlow in cutFlowAll.values():
        if not cutFlow['Fullsample']:
            continue
        # Get total cross-section after pTDM > 150 cut:
        modelDict['Total xsec-pT150 (pb)'] += cutFlow['Fullsample']/(1e3*lumTot)
        # Normalize cutFlow by FullSample:
        for key,val in cutFlow.items():
            if key == 'Fullsample':
                continue
            cutFlow[key] = val/cutFlow['Fullsample']
        cutFlow['Fullsample'] = 1.0
         

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
    dataDict.update({key : [] for key in list(cutFlowAll.values())[0]})

    for ds in dataDict['Data-takingperiod']:
        # Store common values to all datasets:
        for key,val in modelDict.items():
            dataDict[key].append(val)
        # Store dataset-dependent cutflows:
        cutFlow = cutFlowAll[ds]
        for key,val in cutFlow.items():
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
            + " Store the cutflow and SR bins in a pickle (Pandas DataFrame) file." )
    ap.add_argument('-f', '--inputFile', required=True,nargs='+',
            help='path to the ROOT event file(s) generated by Delphes.', default =[])
    ap.add_argument('-o', '--outputFile', required=False,
            help='path to output file storing the DataFrame with the recasting data.'
                 + 'If not defined, will use the name of the first input file', 
            default = None)
    ap.add_argument('-pt', '--pTcut', required=False,default=150.0,type=float,
            help='Gen level MET cut for computing partial cross-sections.')
    ap.add_argument('-n', '--normalize', required=False,action='store_true',
            help='If set, the input files will be considered to refer to multiple samples of the same process and their weights will be normalized.')
    

    ap.add_argument('-v', '--verbose', default='info',
            help='verbose level (debug, info, warning or error). Default is info')


    t0 = time.time()

    # # Set output file
    args = ap.parse_args()
    inputFiles = args.inputFile
    outputFile = args.outputFile
    if outputFile is None:
        outputFile = inputFiles[0].replace('delphes_events.root','cms_exo_20_004.pcl')

    if os.path.splitext(outputFile)[1] != '.pcl':
        outputFile = os.path.splitext(outputFile)[0] + '.pcl'

    modelDict = getRecastData(inputFiles,args.pTcut,args.normalize)

    # #### Create pandas DataFrame
    df = pd.DataFrame.from_dict(modelDict)

    # ### Save DataFrame to pickle file
    print('Saving to',outputFile)
    df.to_pickle(outputFile)

    print("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
