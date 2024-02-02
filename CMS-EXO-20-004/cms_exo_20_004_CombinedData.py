#!/usr/bin/env python3


import pandas as pd
import os,glob


def combineRecastData(files,outputFile):

    # Get files
    if not files:
        print('No valid files found')
    else:
        print('Combining %i files' %len(files))

    allData = pd.read_pickle(files[0])
    for f in files[1:]:
        recastData = pd.read_pickle(f)
        allData = pd.concat((allData,recastData),ignore_index=True)

        
    allColumns = allData.columns.tolist()
    if allData['Coupling'].iloc[0] == 'ADD':
        orderColumns = ['Coupling','Mode','$M_{D}$','$d$','Data-takingperiod']
    elif allData['Coupling'].iloc[0].lower() == 'stop':
        orderColumns = ['Coupling','Mode','$m_{\tilde t}$','$m_{\tilde \chi_1^0}$','Data-takingperiod']
    else:
        orderColumns = ['Coupling','Mode','$m_{med}$','$m_{DM}$','Data-takingperiod']
    allCols = orderColumns[:] + [c for c in allColumns if not c in orderColumns]
    allData = allData[allCols]
    allData.sort_values(orderColumns,inplace=True,
                ascending=[False,False,True,True,False],ignore_index=True)        

    allData.to_pickle(outputFile)



if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser(description=
            "Merge individual DataFrames for model points to a single DataFrame (pickle file).")
    ap.add_argument('-f', '--inputFiles', required=True,nargs='+',
            help='list of pickle files to be merged', default =[])
    ap.add_argument('-o', '--outputFile', required=True, help='output file.')


    args = ap.parse_args()
    inputFiles = args.inputFiles
    outputFile = args.outputFile
    
    combineRecastData(inputFiles,outputFile)







