import xml.etree.ElementTree as ET
import re
import os
from pathlib import Path

def read_SAF(path):
	'''
	Read a SAF file from a pathlib Path, process it to proper XML format and return an ET.root class for parsing
	'''
	with path.open('r') as saf:
		saf_lines = saf.readlines()
		saf_lines.insert(0,"<SAF>\n")
		saf_lines.append("\n</SAF>")
		clean_saf=[line if re.search(r'<.*>',line) else line.replace("&","&amp;").replace(">","&gt;").replace("<","&lt;").replace("'","&apos;").replace('"',"&quot;") for line in saf_lines] 
	with open("tempSAF",'a') as tempfile: # Easiest way to return ET.parser is generating a temp clean file with escaped characters
		for line in clean_saf:
			tempfile.write(line)
	tree = ET.parse("tempSAF")
	root = tree.getroot()
	os.remove("tempSAF")
	return root

def get_cuts(root):
	cuts=[]
	for entry in root.iter("InitialCounter"):
		cuts.append(re.search(r'"(.*)"', entry.text).group(1))
	for entry in root.iter("Counter"):
		cuts.append(re.search(r'"(.*)"', entry.text).group(1))
	return cuts

def get_entries(root):
	nevents=[]
	for entry in root.iter("InitialCounter"):
		nevents.append(re.search(r'(\d+).*#\s+nentries.*', entry.text).group(1))
	for entry in root.iter("Counter"):
		nevents.append(re.search(r'(\d+).*#\s+nentries.*', entry.text).group(1))
	return nevents

def get_weights(root):
	nweights=[]
	for entry in root.iter("InitialCounter"):
		nweights.append(re.search(r'(\d\.\d+e.\d+).*#\s+sum of weights', entry.text).group(1))
	for entry in root.iter("Counter"):
		nweights.append(re.search(r'(\d\.\d+e.\d+).*#\s+sum of weights', entry.text).group(1))
	return nweights
