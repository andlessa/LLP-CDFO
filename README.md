# LLP-CDFO

Holds the code and data for reproducing the results in the [Probing conversion-driven freeze-out at the LHC](https://arxiv.org/abs/2404.16086) paper.



## Installation

The following external packages are needed:
  
  * [MadGraph5](https://launchpad.net/mg5amcnlo) with [LHAPDF6](https://lhapdf.hepforge.org/), [Pythia8](https://pythia.org/)
  * [MadAnalysis5](https://launchpad.net/madanalysis5)
  * [DelphesLLP](DelphesLLP.tar.gz): modified version of [Delphes](https://cp3.irmp.ucl.ac.be/projects/delphes) with some additional modules to collect information relevant for recasting LLP searches.
  
Running:

```
./installer.sh
```

Should install the relevant packages.



## Folders and files

Below we describe the main files and folders stored in this repository.

* Folders:

  * [ATLAS-SUSY-2016-08](ATLAS-SUSY-2016-08), [ATLAS-SUSY-2018-13](ATLAS-SUSY-2018-13), [ATLAS-SUSY-2018-42](ATLAS-SUSY-2018-42), [CMS-EXO-19-010](CMS-EXO-19-010) and  [CMS-EXO-20-004](CMS-EXO-20-004): folders containing the relevant recasting code for each analysis (except for CMS-EXO-19-010 which relies on MadAnalysis5)
  * [statisticalTools](statisticalTools): folder with statistical tools used for computing limits and the high luminosity projections
  * [Cards](Cards): relevant cards for event generation
  * [plotting](plotting): folder containing all the plotting notebooks.
  * [results_dataFrames](results_dataFrames): folder with pickle files containing Pandas DataFrames with the recasting results (see below)
  * [eventData](eventData): folder storing Delphes ROOT files for two benchmark points (see [README](eventData/README.md)).

* Files:

  * [CDFOdata_2112_01499v3_Fig9_Good.dat](CDFOdata_2112_01499v3_Fig9_Good.dat): points in the model parameter space used for the scans
  * [DelphesLLP.tar.gz](DelphesLLP.tar.gz): modified Delphes version with some additional modules to collect information relevant for recasting LLP searches. The FastJetFinder module was also modified to compute the summed pT of displaced tracks.
  * [helper.py](helper.py): several convenience functions for dealing with LLP searches.
  * [ParticleData.xml](ParticleData.xml): modified ParticleData file to include R-hadrons.
  * [runScanMG5_hepmc.py](runScanMG5_hepmc.py): scan script for generating HepMC events using MadGraph.
  * [runScanMG5.py](runScanMG5.py): scan script for generating MC events using MadGraph and the DelphesPythia8 interface.
  * [scan_parameters_cdfo_atlas.ini](scan_parameters_cdfo_atlas.ini),[scan_parameters_cdfo_cms.ini](scan_parameters_cdfo_cms.ini) and [scan_parameters_cdfo_hepmc.ini](scan_parameters_cdfo_hepmc.ini): examples for parameter files which can be used as input to the scan scripts (runScanMG5.py or runScanMG5_hepmc.py).
  * [setenv.sh](setenv.sh): convenience script for setting the necessary system variables
  * [WIMPregion_DDconstr_LZ.dat](WIMPregion_DDconstr_LZ.dat): file with direct detection constraints
  * [xsecs_sbottom.csv](xsecs_sbottom.csv): file with NLO cross-sections for sbottoms, which are used to compute k-factors.



## CDFO - LHC Constraints

The recasting of the following searches have been used to compute the LHC constraints:

 * [ATLAS-SUSY-2016-08](https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2016-08/): ATLAS Displaced Vertex plus MET search
 * [ATLAS-SUSY-2018-13](https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-13/): ATLAS Displaced Vertex plus jets search
 * [ATLAS-SUSY-2018-42](https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2018-42/): ATLAS dE/dx (HSCP) search
 * [CMS-EXO-19-010](https://cms-results.web.cern.ch/cms-results/public-results/publications/EXO-19-010/): CMS disappearing track search
 * [CMS-EXO-20-004](https://cms-results.web.cern.ch/cms-results/public-results/publications/EXO-20-004/): CMS jets plus MET prompt search

Each recast folder contains the following structure:

 * `<data folder>`: folder containing HepData files for the analysis and auxiliary code for extracting the relevant information.
 * `<analysis_name>_Recast.py`: main code for running the recasting using as input a Delphes ROOT file.
 * `<analysis_name>_Cutflow.py`: code for obtaining the cutflow using as input a Delphes ROOT file.
 * `<analysis_name>_CombineData.py`: convenience code for merging several pickle files containning the output of the recast code.
 * `<analysis_name>_UpperLimits.py`: main code for computing 95% C.L. limits using as input the pickle file generated by the recast code.
 * `<analysis_name>_UpperLimits_HighLuminosity.py`: main code for computing 95% C.L. limits using the high luminosity projections.
 
The only exception is the recasting of the disappearing track search ([CMS-EXO-19-010](CMS-EXO-19-010)), which is based on the MadAnalysis5 implementation.

### Pipeline for running the recasting

Detailed information for how to run each recasting code can be obtaining running the main recast code with the `--help` flag. Below we schematically list the main steps required for reproducing the results in [Probing conversion-driven freeze-out at the LHC](https://arxiv.org/abs/2404.16086):

 1. Generate MC events and the Delphes ROOT output  for the model points running:

    ```
    runScanMG5.py -p <scan_parameters_file>
    ```
    the output will be a MadGraph5 folder containing Delphes ROOT files for each scan point.


 2. Within a given analysis folder (e.g. [ATLAS-SUSY-2016-08](ATLAS-SUSY-2016-08)) compute the signal yields and efficiencies for the scan points running:
   
    ```
    ./<analysis_name>_Recast.py -f <events_folder>/run_*/*.root
    ```
    
    the output will be a single pickle file for each run containing the results for all the scan points. In case there are multiple run folders for the same model point, the respective events will be combined.

 3. Merge the output pickle files generated in the previous step:
 
    ```
    ./<analysis_name>_CombineData.py -f <events_folder>/run_*/*.pcl -o ../results_dataFrames/<combined_output>
    ```

    the output will be a single pickle file containing the results for all the scan points.

 4. Compute the upper limits using as input the combined output file generated in the previous step:
 
    ```
    ./<analysis_name>_UpperLimits.py -f ../results_dataFrames/<combined_output>
    ```

    the upper limits will be added to the input DataFrame and stored in the same pickle file.


 5. Finally, the generated pickle files can be used to plot the exclusions using the [plots-AllExclusions.ipynb](plotting/plots-AllExclusions.ipynb) notebook.

