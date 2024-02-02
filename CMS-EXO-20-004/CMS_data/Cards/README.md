# CMS Cards

Generation cards provided by CMS for the CMS-SUS-20-004 search.

The cards are for three distinct models:

 * Spin 0 Mediator:
   * scalar scenario
   * pseudo-scalar scenario
 * Spin 1 mediator:
   * axial scenario
   * vectorial scenario
 * t-channel Mediator
 * Leptoquarks
 * ADD
 
 
## Simplified DM and Leptoquark models (MadGraph based)

We provide the following cards:

* run_card.dat: typical madgraph, contains run settings, predefined cuts, etc
* proc_card.dat: typical madgraph, contains the process definition (i.e. "generate" statement)
* customizecards.dat: Specific commands to set parameters for a given signal point
* pythia.cmnd: Pythia shower configuration.

Generally, the logic of the cards is that all parameters are left at the default value in the UFO model, unless explicitly changed in the customizecards file.
If you want to make a parameter scan, leave the proc and run cards as they are and modify the desired parameters as shown in customizecards.

### Spin-1 Mediator

In this case, there are six separate sample types: Monojet, mono-Z and mono-W, each with vector and axial couplings, respectively.
In order to get the full prediction for one parameter point, combine the three samples (monojet, mono-Z and mono-W).

In addition to the typical cards, we also provide a [DMSimp_NLO_cuts.f](./DMSimp_NLO_cuts.f) file that implements a madgraph-level cut on the pt of the DM pair system.
This cut allows to significantly enhance the generation efficiency with respect to the signal selection.


### t-channel Mediator

For each parameter point, we have split the signal generation into two parts:

 1. Final state of two jets + DM pair from mediator pair production (PhiPhiToJJChiChi)

 2. Final state of single jet + DM pair, where the jet is from QCD raidation (JChiChi)

the full signal prediction is simply the sum of the two samples. 
One might also generate them in a single sample for simplicity.



## ADD model (Pythia-based)

We provide the following cards:

* pythia.cmnd: Combined model and shower configuration file. We explicitly mark the places where parameters of the BSM model should be modified for parameter scans.

