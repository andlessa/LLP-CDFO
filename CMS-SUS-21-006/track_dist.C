
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include <stdio.h>      /* printf */
#include <math.h>       /* atan */
#include <complex> 
#else
class ExRootTreeReader;
class ExRootResult;
#endif


//------------------------------------------------------------------------------
double track_length(GenParticle* mother, GenParticle* daughter)
{
	double mx = mother->X, my=mother->Y, mz=mother->Z, dx=daughter->Y, dy=daughter->Y, dz=daughter->Z;
	return sqrt( (mx-dx)*(mx-dx) + (my-dy)*(my-dy) + (mz-dz)*(mz-dz) );
}


void track_dist(const char *inputFile) //REMEMBER TO CHANGE THE NAME TO THE MACRO FILE NAME!!!!!!
{

  gSystem->Load("libDelphes");

  TChain chain("Delphes");  //name of the tree
  chain.Add(inputFile);     // we can add several files this way

    //Create chain of root trees
//  TChain chain("LHEF");
//  TString  rootchain;
 // chain.Add("LHEF");


  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
//  TClonesArray *branchFatJet = treeReader->UseBranch("FatJet");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMet = treeReader->UseBranch("MissingET");
  TClonesArray *branchGenMissingET = treeReader->UseBranch("GenMissingET");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  //test of weights
  //TClonesArray *branchWeight = treeReader->UseBranch("Weight");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");  //AQUI

  //several output files
  TFile f("myoutputs/track_len_dist.root","new");
  TCanvas* c = new TCanvas("c","placeholder");

  vector<const Electron*>   *electrons = new vector<const Electron*>();
  vector<const Muon*>       *muons     = new vector<const Muon*>();
  vector<const MissingET*>  *missinget = new vector<const MissingET*>();
  vector<const Jet*>        *jets      = new vector<const Jet*>();

// Book histograms
//  TH1 *histVetCoord = new TH1F("vertex_coord", "Decay vertices distribution;Time(s)", 100, 0.0, 1000.0);
  TH1 *histTrackLen = new TH1F("track_length", "Chargino track length distribution;Length(mm)", 100, 0.0, 1000.0);

  // Loop over all events

  for(int entry = 0; entry < allEntries; ++entry)
  {
    treeReader->ReadEntry(entry);     // Load selected branches with data from specified event

    Photon *photon;
    Muon *muon;
    Electron *electron;
    MissingET *missinget;
    TLorentzVector lv;
  


    HepMCEvent *event = (HepMCEvent*) branchEvent -> At(0); //AQUI

    //clearing some variables
    electrons -> clear();
    muons     -> clear();


    // loop over generated  particles in the event
    for(Int_t i=0; i < branchParticle->GetEntriesFast(); i++)
      {    
	GenParticle *gen = (GenParticle*) branchParticle->At(i);

	if( gen->PID == 1000024   )
	{       
		GenParticle *daughter = (GenParticle*) branchParticle->At(gen->D1);
		histTrackLen->Fill(track_length(gen,daughter));
	}

      } //end of the loop in the generated particles


  } //end of the loop

histTrackLen->Draw();
c->SaveAs("hist_track_len.png");
f.Write();

cout<< "Reached end of analysis!"<< endl;

}

      


  




