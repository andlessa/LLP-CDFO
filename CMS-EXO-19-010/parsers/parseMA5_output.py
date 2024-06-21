#!/data/01/lucasmdr/env_py3/bin/python3

from pathlib import Path
import pandas as pd
import numpy as np
import glob
import re
import SAFreader as sfr
import pickle

pd.set_option('display.max_columns', 11)
pd.options.mode.chained_assignment = None #Disable copy warnings

def HL_uL(SR_name):
	HL_sig95obs = {"SR1_2017":0.21146,"SR2_2017":0.03347,"SR3_2017":0.06237,
	"SR1_2018A":0.28914, "SR2_2018A":0.02466, "SR3_2018A":0.02169,
	"SR1_2018B":0.2268497, "SR2_2018B":0.0160657, "SR3_2018B":0.0670251}
	return HL_sig95obs[SR_name]

def luminosity(SR_name):
	lumis = {"SR1_2017":42.0,"SR2_2017":42.0,"SR3_2017":42.0,
	"SR1_2018A":21.0, "SR2_2018A":21.0, "SR3_2018A":21.0,
	"SR1_2018B":39.0, "SR2_2018B":39.0, "SR3_2018B":39.0}
	return lumis[SR_name]

def HL_uL_dR(SR_name):
	HL_sig95obs = {"SR1_2017_dR":0.21146,"SR2_2017_dR":0.03347,"SR3_2017_dR":0.06237,
	"SR1_2018A_dR":0.28914, "SR2_2018A_dR":0.02466, "SR3_2018A_dR":0.02169,
	"SR1_2018B_dR":0.2268497, "SR2_2018B_dR":0.0160657, "SR3_2018B_dR":0.0670251}
	return HL_sig95obs[SR_name]

def luminosity_dR(SR_name):
	lumis = {"SR1_2017_dR":42.0,"SR2_2017_dR":42.0,"SR3_2017_dR":42.0,
	"SR1_2018A_dR":21.0, "SR2_2018A_dR":21.0, "SR3_2018A_dR":21.0,
	"SR1_2018B_dR":39.0, "SR2_2018B_dR":39.0, "SR3_2018B_dR":39.0}
	return lumis[SR_name]

def combine_Eff(df,eff):
	regions = {"SR1_2017":'SR1',"SR1_2018A":'SR1',"SR1_2018B":'SR1',
	"SR2_2017":'SR2', "SR2_2018A":'SR2',"SR2_2018B":'SR2',
	"SR3_2015":'SR3',"SR3_2016A":'SR3',"SR3_2016B":'SR3',"SR3_2017":'SR3',"SR3_2018A":'SR3', "SR3_2018B":'SR3',
	"SR1_2017_dR":'SR1',"SR1_2018A_dR":'SR1',"SR1_2018B_dR":'SR1',
	"SR2_2017_dR":'SR2', "SR2_2018A_dR":'SR2',"SR2_2018B_dR":'SR2',
	"SR3_2015_dR":'SR3',"SR3_2016A_dR":'SR3',"SR3_2016B_dR":'SR3',"SR3_2017_dR":'SR3',"SR3_2018A_dR":'SR3', "SR3_2018B_dR":'SR3'}

	comb_SR = regions[df['SR']]
	if comb_SR in eff.keys():
		eff[comb_SR] += df['Eff']
	else:
		eff[comb_SR] = df['Eff']

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser(description='Generate pickle files with pd.DF cutflows and upper limits from a MA5 Recast mode output')
	parser.add_argument('SAF', metavar='SAF_dir_path', help='Path to the Analysis directory containing MA5 results folders') #required=True if not positional
	parser.add_argument('output', metavar='Output_dir_path', help='Path to the output directory to write the pickle files')
	args = parser.parse_args()

	SAF = args.SAF
	output = args.output
	
	results = [] # Initialize list for UL results and dictionary for Cutflow output
	results_dR = []
	hilumi_results = []
	hilumi_results_dR = []
	
	SRs = {}
	ctau_pat = re.compile(r'\s+1000005\s+([0-9]+\.[0-9]+e.[0-9]+).*wsd3') #Compile regex patterns once to parse parameters from MG5 output for each point
	msb_pat = re.compile(r'\s+1000005\s+([0-9]+\.[0-9]+e.[0-9]+).*')
	mneu1_pat = re.compile(r'\s+1000022\s+([0-9]+\.[0-9]+e.[0-9]+).*')
	xsec_pat = re.compile(r'[0-9]+\s+(.*)\s+[0-9]+.*')
	print("Arguments read: %s and %s" %(SAF, output))
	write_path = Path(args.output) #Initialize path to folder to write results, and create any parents if needed
	write_path.mkdir(parents=True, exist_ok=True)

	for base_path in glob.glob(args.SAF+"/**/", recursive=False):

		in_path=Path(base_path) / "Input" # Initialize paths to Input and Output dirs in the MA5 results folder for the point
		out_path=Path(base_path) / "Output"
		
		param_path = next(Path(next(in_path.glob("*.list")).read_text().strip()).parent.glob("*banner.txt"))       # Get read .list file with path to all input hepmc files
		xsec_path = next(Path(next(in_path.glob("*.list")).read_text().strip()).parent.glob("*merged_xsecs.txt"))  # Assumes 1 dataset (point) per MA5 run, generated by MG5
		xsec = float(re.search(xsec_pat,xsec_path.read_text())[1])                                                 # with 1j to generate *merged_xsecs.txt files
		
		msb=-1
		mneu1=-1
		ctau=-1
		
		with param_path.open('r') as paramfile:
			for line in paramfile:          #Loop over file to match the regex patters line by line, stops once all are matched for performance
				if msb==-1:
					msb_match = re.search(msb_pat, line)
				if mneu1==-1:
					mneu1_match = re.search(mneu1_pat, line)
				if ctau==-1:
					ctau_match = re.search(ctau_pat, line)
				if msb_match is not None:
					msb = float(msb_match[1])
				if mneu1_match is not None:
					mneu1 = float(mneu1_match[1])
				if ctau_match is not None:
					ctau = 6.582e-16/float(ctau_match[1])
				if msb!=-1 and mneu1!=-1 and ctau!=-1:
					break
		
		CLpath = out_path / "SAF" / "CLs_output_summary.dat" #Initialize path to the CLs summary file, and reads it with CSV reader
		CLdf = pd.read_csv(CLpath, header=None,skiprows=1, sep="\s+\|\|\s+|\s+|\|\|", names=["Point", "Analysis", "SR", "Sig95(exp)", "Sig95(obs)", "Eff", "Unc"], engine='python')

		CLaux = CLdf[CLdf["SR"].str.contains(".*(?:2015|2016A|2016B|2017|2018A|2018B)$",regex=True,na=False)]
		CLaux_dR = CLdf[CLdf["SR"].str.contains(".*(?:2015|2016A|2016B|2017|2018A|2018B)_dR$",regex=True,na=False)]

		combinedEff = {}
		combinedEff_dR = {}
		CLaux.apply(combine_Eff,args=(combinedEff,),axis=1) #Computes the overall efficiency for each combined signal region, in order to obtain the upper limit on the production cross section
		CLaux_dR.apply(combine_Eff,args=(combinedEff_dR,),axis=1)

		combinedULs = {'SR1': 0.21531656362673005, 'SR2': 0.06224022369369497, 'SR3':0.05986366949082641}
		combinedRobs = np.array([1000*xsec/(combinedULs[sr] / combinedEff[sr]) if combinedEff[sr]!=0 else 0 for sr in ['SR1', 'SR2', 'SR3']]) #Computes the robs choosing the most stringent among the three combined signal regions
		combinedRobs_dR = np.array([1000*xsec/(combinedULs[sr] / combinedEff_dR[sr]) if combinedEff_dR[sr]!=0 else 0 for sr in ['SR1', 'SR2', 'SR3']])

		results.append({'mLLP':msb, 'mLSP':mneu1, 'deltaM':msb-mneu1, 'tau_ns':ctau, 'Total xsec (pb)':xsec, 'robs':combinedRobs.max()}) #Append dict with extracted results to build DF later
		results_dR.append({'mLLP':msb, 'mLSP':mneu1, 'deltaM':msb-mneu1, 'tau_ns':ctau, 'Total xsec (pb)':xsec, 'robs':combinedRobs_dR.max()}) #Append dict with extracted results to build DF later


		#hilumi_CLaux = CLdf[CLdf["SR"].str.contains(".*(?:2017|2018A|2018B)$",regex=True,na=False)]
		hilumi_CLaux = CLdf[CLdf["SR"].str.contains(".*2018A$",regex=True,na=False)]
		hilumi_CLaux_dR = CLdf[CLdf["SR"].str.contains(".*2018A.*_dR$",regex=True,na=False)]

		hilumi_CLaux['HL_sig95_obs'] = hilumi_CLaux['SR'].apply(HL_uL)
		hilumi_CLaux['lumi'] = hilumi_CLaux['SR'].apply(luminosity)
		hilumi_CLaux['robs'] =1000*xsec*hilumi_CLaux['Eff']*140.4/(hilumi_CLaux['lumi']*hilumi_CLaux['HL_sig95_obs'])
		hilumi_CLaux['Nup']= 1000*hilumi_CLaux['Sig95(obs)']*hilumi_CLaux['Eff']*140.4

		hilumi_CLaux_dR['HL_sig95_obs'] = hilumi_CLaux_dR['SR'].apply(HL_uL_dR)
		hilumi_CLaux_dR['lumi'] = hilumi_CLaux_dR['SR'].apply(luminosity_dR)
		hilumi_CLaux_dR['robs'] =1000*xsec*hilumi_CLaux_dR['Eff']*140.4/(hilumi_CLaux_dR['lumi']*hilumi_CLaux_dR['HL_sig95_obs'])
		hilumi_CLaux_dR['Nup']= 1000*hilumi_CLaux_dR['Sig95(obs)']*hilumi_CLaux_dR['Eff']*140.4

		hilumi_results.append({'mLLP':msb, 'mLSP':mneu1, 'deltaM':msb-mneu1, 'tau_ns':ctau, 'Total xsec (pb)':xsec, 'robs':hilumi_CLaux['robs'].max()}) #Append dict with extracted results to build DF later
		hilumi_results_dR.append({'mLLP':msb, 'mLSP':mneu1, 'deltaM':msb-mneu1, 'tau_ns':ctau, 'Total xsec (pb)':xsec, 'robs':hilumi_CLaux_dR['robs'].max()}) #Append dict with extracted results to build DF later
		
		#Cutflow parsing
		cutflows = glob.glob(base_path+"/**/Cutflows/*.saf", recursive=True) #Find all Cutflows in SAF files

		for file in cutflows:
			cutflow_file = Path(file)
			SRname = cutflow_file.name.replace(".saf","") # Parse SR name from filename
			aux = {'mLLP' : msb, 'mLSP' : mneu1, 'tau_ns' : ctau}    #Initialize dict with the mass values identifying the point

			reader = sfr.read_SAF(cutflow_file) # Read cutflow file and extract cuts and surviving events
			cuts = sfr.get_cuts(reader)
			events = sfr.get_entries(reader)
			weights = [float(value) for value in sfr.get_weights(reader)]

			aux.update(dict(zip(cuts,weights))) # Updates dict with cutflow
			if SRname not in SRs:    # Check if SR is already in the dictionary; if not, initialize it
				SRs[SRname] = []
			SRs[SRname].append(aux) # Append cutflow dict to output list to build DF later

	pickle_path = write_path / "pp2BB1j_cms_exo_19_010.pcl"
	results_df = pd.DataFrame(sorted(results, key=lambda d: d['mLLP'])).to_pickle(pickle_path) # Builds the Exclusion DF and saves it to pickle file in the output folder

	pickle_path_dR = write_path / "pp2BB1j_cms_exo_19_010_dRcut.pcl"
	results_df_dR = pd.DataFrame(sorted(results_dR, key=lambda d: d['mLLP'])).to_pickle(pickle_path_dR) # Builds the Exclusion DF and saves it to pickle file in the output folder

	
	hilumi_pickle_path = write_path / "pp2BB1j_cms_exo_19_010_HL.pcl"
	hilumi_results_df = pd.DataFrame(sorted(hilumi_results, key=lambda d: d['mLLP'])).to_pickle(hilumi_pickle_path)

	hilumi_pickle_path_dR = write_path / "pp2BB1j_cms_exo_19_010_dRcut_HL.pcl"
	hilumi_results_df_dR = pd.DataFrame(sorted(hilumi_results, key=lambda d: d['mLLP'])).to_pickle(hilumi_pickle_path_dR)

	cutflow_folder_path = write_path / "Cutflows"
	cutflow_folder_path.mkdir(exist_ok=True)
	for sr,cutflow in SRs.items():
		pickle_path = cutflow_folder_path / (sr+".pcl")
		pd.DataFrame(sorted(cutflow, key=lambda d: d['mLLP'])).to_pickle(pickle_path) # Builds the Cutflow DF and saves it to pickle file in the output folder for each signal region
	print("Parsing has finished!")

