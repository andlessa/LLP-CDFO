#!/bin/bash

#Simples BASH script to read MG5 results generated from runScanMG5_hepmc.py, and generate .ma5 input files for MadAnalysis5 scripted runs. Might be replaced by python implementation
#Usage: ./generate_ma5.sh MG5_proc_dir output_dir

echo -e "Starting to parse MG5 outputs... "

procDIR=${1:-.};
outpath=${2:-.};
recast_card_path=${3:-../MA5_input/recasting_card.dat}

if ! [[ -d $outpath ]]; then # Check if output directory exists; if not, creates it
	mkdir -p $outpath
	fi

outputs=()

for path in $(for name in $(find ${procDIR} -name '*.hepmc.*' | sort); do readlink -f ${name}; done);do # Finds the paths to all hepmc outputs from MG5 runs
	key=$(basename ${path} | sed -r 's/\_tau.*//g') # Assumes filenames *(%?.f msb)_(%d mn1)_tau_(ctau)*.hepmc.* ; each point is identified by *msb_mn1
	outfile=${outpath}/${key}.ma5 # Generates path to .ma5 output. Make sure there is no point ambiguity in the output format of msb_mn1 when running runScanMG5 beforehand
	if ! [ -f ${outfile} ]; then # Initialize new .ma5 file with default recast options
		echo -e "set main.recast = on" >> ${outfile};
		echo -e "set main.recast.store_root = False" >> ${outfile};
		echo -e "set main.recast.card_path = ${recast_card_path}" >> ${outfile};
		outputs+=($outfile)
	fi
	echo -e "import ${path} as $(basename ${path} | sed -r 's/_tau.*//g' | sed -r 's/\./p/g')" >> $outfile; # Writes a MA5-compliant dataset name replacing . for p
done

for file in ${outputs[@]}; do
	echo -e "submit /data/01/lucasmdr/cms_exo_19_010/outputs/$(basename ${file} | sed -r 's/.ma5//g')" >> ${file}; # Each (name).ma5 will generate a (name) output dir
done
echo -e "Finished!"



