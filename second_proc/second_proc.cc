//g++ -O2 `$ROOTSYS/bin/root-config --libs` -I$ROOTSYS/include second_proc.cc -o second_proc


#include <stdlib.h>
//#include <crtdbg.h>
//  ROOT includes
//
#include "TInterpreter.h"
#include "TObject.h"
#include "TChain.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TText.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TVector3.h"
//
//  C/C++ includes
//
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include <set>
//
#include "NMMEvent.hh"

using namespace std;

//int main(int argc, char* argv[]) 
void second_proc(const string infilename, const string outfilename)
{

  if (!TClass::GetDict("NMMEvent")) {
    gROOT->ProcessLine(".L NMMEvent.cc+");
    //  gInterpreter->GenerateDictionary("NMMEvent","NMMEvent.hh");
  }

  //string infilename = argv[1];
  //string outfilename = argv[2];

  TFile *fin = new TFile(infilename.c_str(),"READ");
  // "melindas_neutrons_PulseTags.root"

  //output ttree
  TFile *fout = new TFile(outfilename.c_str(),"RECREATE");
  TTree * outEventTree = new TTree("Event",
				   "Event tree with the second process");
  NMMEvent outEventFiller;
  cout<<"Building the Branch.\n";
  outEventTree->Branch("NMMEvent","NMMEvent",&outEventFiller,32000,1);
  cout<<"Branch has been built.\n";
  
  TTree *event;
  fin->GetObject("Events",event);
  if (event == 0) {
    delete fin;
    return;
  }
  
  cout<<"Input TTree found.\n";
  vector<double> *pt, *ph, *nf;
  vector<string> *ptc, *ptd, *ptg, *pcv;
  event->SetBranchAddress("mergedPulseTime",&pt);
  event->SetBranchAddress("mergedPulseHeight",&ph);
  event->SetBranchAddress("northFraction",&nf);
  event->SetBranchAddress("pulseTagCap",&ptc);
  event->SetBranchAddress("pulseTagDet",&ptd);
  event->SetBranchAddress("pulseTag",&ptg);
  event->SetBranchAddress("pulseTagCav",&pcv);

  cout<<"Input branches hooked.\n";

  int numEntries = event->GetEntries();
  cout<<"Entry Num = " << numEntries << "\n";
  for (int n=0;n<numEntries;n++) {

    event->GetEntry(n);

    if (n%1000000==0)
      cout<< "Processing event "<< n <<"\n";
    //vector<double> pulseTime = *pt;
    //vector<double> pulseHeight = *ph;
    vector<string> pulseTagCap = *ptc;
    vector<string> pulseTagDet = *ptd;
    vector<string> pulseTag = *ptg;
    vector<string> pulseTagCav = *pcv;

    int M;
    //int mRaw = pulseTag.size();
    int numWaterCaptures = 0;
    int numAcrylicCaptures = 0;
    int numLeadCapturesPromptBKG = 0;
    int numLeadCaptures = 0;
    int numCavBkg = 0;

    bool muonFlag = false;
    bool otherbkgFlag = false;

    NMMEvent ev;
    ev.N = pulseTag.size();
    ev.pulseHeights = *ph;
    ev.pulseTimes = *pt;
    ev.northFraction = *nf;
    
    if (ev.N != 0) {
      // only process an event if not empty
      for (int i = 0; i < ev.N; i++) {
	if (pulseTag[i] == "CapWater") {
	  numWaterCaptures++;
	  ev.pulseTagV2.push_back("CapWater");
	}
	else if(pulseTag[i] == "Unknown") { // should be acrylic captures
	  numAcrylicCaptures++;
	  ev.pulseTagV2.push_back("CapAcrylic");	  
	}
	else if (pulseTagCav[i] == "FromCav") {
	  numCavBkg++;
	  if (pulseTag[i] == "muon") {
	    ev.pulseTagV2.push_back("MuonFromCav");
	    muonFlag = true;
	  }
	  else {
	    ev.pulseTagV2.push_back("OtherFromCav");
	    otherbkgFlag = true;
	  }      
	}
	else {
	  numLeadCapturesPromptBKG++;
	  
	  if (ev.pulseTimes[i] > 1.3){
	    numLeadCaptures++;
	    ev.pulseTagV2.push_back("CapLead");
	  }
	  else {
	    ev.pulseTagV2.push_back("promptX"); // unknown prompt background
	    otherbkgFlag = true;
	  }
	}
      }
      
      ev.M = numWaterCaptures + numAcrylicCaptures + numLeadCaptures;
      
      if (ev.N == ev.M)
	ev.type = "pfn"; // pure fast neutron
      else if (ev.M != 0){
	if(muonFlag)
	  ev.type = "contam_muon";
	else if(otherbkgFlag)
	  ev.type = "contam_other";
	else
	  cout << "Error on contaminated event type!\n";
      }
      else {
	if(muonFlag)
	  ev.type = "bkg_muon";
	else if(otherbkgFlag)
	  ev.type = "bkg_other";
	else
	  cout << "Error on pure bkg event type!\n";
      }
    }
    //output
    if(ev.type != "emptyEvent"){
      outEventFiller = ev;
      outEventTree->Fill();
    }

  }// events

  
  fout = outEventTree->GetCurrentFile();
  fout->Write();
  fout->Close();
  
}
