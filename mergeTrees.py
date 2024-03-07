#!/usr/bin/env python3

# Merge Delphes trees


import os
import progressbar as P
import glob
import multiprocessing

import ROOT
from ROOT import TFile, TTree, TList

delphesDir = os.path.abspath("./DelphesLLP")
os.environ['ROOT_INCLUDE_PATH'] = os.path.join(delphesDir,"external")

ROOT.gSystem.Load(os.path.join(delphesDir,"libDelphes.so"))

ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')


def mergeFileList(pathList):

    treeList = TList()
    pyfilelist = []
    pytreelist = []
    
    for path in pathList:
        inputFile = TFile(path, 'read')
        pyfilelist.append(inputFile) # Make this TFile survive the loop!
        inputTree = inputFile.Get('Delphes')
        pytreelist.append(inputTree) # Make this TTree survive the loop!
        treeList.Add(inputTree)

    outName = pathList[0].replace('_delphes_events.root','_0j1j_delphes_events.root')
    outputFile = TFile(outName, 'recreate')
    outputFile.cd()
    outputTree = TTree.MergeTrees(treeList)
    outputFile.Write()
    outputFile.Close()


if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Combine ROOT trees (from Delphes output) from folders.")
    ap.add_argument('-f', '--inputFolders', required=True,nargs='+',
            help='path to folder containing Delphes ROOT trees.', default =[])    
    

    args = ap.parse_args()
    inputFolders = args.inputFolders

    foldersList = []
    for f in sorted(glob.glob(os.path.join(inputFolders[0],'Events/run_*/*.root'))):
    
        pathList = [f.replace(inputFolders[0],inputF) 
                for inputF in inputFolders]
        foldersList.append(pathList)

    progressbar = P.ProgressBar(widgets=["Reading %i folders: " %len(foldersList), 
                                P.Percentage(),P.Bar(marker=P.RotatingMarker()), P.ETA()])
    progressbar.maxval = len(foldersList)
    progressbar.start()

    ntotal = 0
    children = []
    ncpus =  multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=ncpus)
    for pathList in foldersList:
        
        p = pool.apply_async(mergeFileList, args=(pathList,))                       
        children.append(p)

#     Wait for jobs to finish:
    for p in children:
        p.get()
        ntotal += 1
        progressbar.update(ntotal)
        

    progressbar.finish()