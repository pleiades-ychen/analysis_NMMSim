/*
 * g++ -O2 `$ROOTSYS/bin/root-config --libs` -I$ROOTSYS/include proc_tagpulses_deteff.cc -o proc_tagpulses_deteff_v2
 */
// modified special version for the study of fast-neutron detection efficiency 
// 07/17/2015
// v2 07/21/2015

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
//#include <crtdbg.h>
//  ROOT includes
//
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
#include <iostream>
#include <set>
//
//  Definitions
//
#define DEBUGGING 0

#define pulseWidth 0.25

using namespace std;

class PMT_Hits {
public:
  int num;
  double *E, *T, *X, *Y;

  PMT_Hits (int);
  ~PMT_Hits ();
  void LinkTree (TTree *);
};

class MC_Truth {
  char name[10];
  char Name[10];
public:
  int num;
  double *X, *Y, *Z, *T;
  double *Px, *Py, *Pz, *E;
  int *TID, *PID, *PDG;

  MC_Truth (const char *, int);
  ~MC_Truth ();
  void LinkTree (TTree *);
};

class CaptureInWater {
  char name[10];
  char Name[10];
public:
  int num;
  double *X, *Y, *Z, *T;

  CaptureInWater (int);
  ~CaptureInWater ( );
  void LinkTree (TTree *);
};
// Add Primary Class for muon track analysis
class Primary {
  char name[10];
  char Name[10];
public:
  int num;
  double *E, *Px, *Py, *Pz;
  double *X, *Y, *Z;
 
  Primary (int);
  ~Primary ();
  void LinkTree (TTree *);
};


struct MC_Track {
  double X, Y, Z, T;
  double Px, Py, Pz, E;
  int TID, PID, PDG;
  int TrigPulse;
};

bool comp_energy(MC_Track a, MC_Track b) { return (a.E > b.E); }
// used at 3rd argument of sorting vector<MC_Track> by energy descendingly
// indicate that the larger one goes first when sorting
bool comp_time(MC_Track a, MC_Track b) { return (a.T < b.T); }
// in order to sort by time ascendingly
double get_mass (int pdg);
double get_cer_thr (int pdg);
void cut_cer_thr(vector<MC_Track> &);
vector<MC_Track> convert_mctruth(const MC_Truth &);
vector<MC_Track> preproc_mctruth(const MC_Truth &);
bool PmtHitsToPhotoElectrons(double photonWavelength);
void collect_hits(vector< vector<double> > &pmtHitTimes, const PMT_Hits &pmtHits);
void get_pmt_pulses(vector<double> &pmtPulseHeights, vector<double> &pmtPulseTime,
		    const vector<double> &pmtHitTimes, const int iPMT);
void get_coin_pulses(vector<double> &, vector<double> &, 
		     const vector< vector<double> >::iterator,
		     const vector< vector<double> >::iterator);
void merge_tank_pulses (vector<double> &, vector<double> &, vector<double> &,
			const vector<double> &, const vector<double> &,
			const vector<double> &, const vector<double> &);
// geometric tools for tagging muon hitting lower panels
class TrackLine {
    double x0,y0,z0; // primary position
    double dx,dy,dz; // direction vector
    
public:
    TrackLine(double,double,double,double,double,double);
    TVector3 GetPoint(double,char);
};

class Rectangle {
  double center[3];
  double side[2];
    
public:
  Rectangle();
  Rectangle(double,double,double,double,double);
  // i = 0,1; for the two dimentions of Rectangle
  // Bottom(i) and Top(i) mean Lower and Upper limits of the rectangle's ranges for each dimensions.
  double Z(){return center[2];}
  bool Inside(double,double);
  double Bottom(int i){return (center[i]-0.5*side[i]);}
  double Top(int i){return (center[i]+0.5*side[i]);}
};

bool Traverse(TrackLine, Rectangle);
////////////////////////////////////////////////////////

int match (double, double); 


//vector<string> GetPulseTag(const vector<double> *, const vector<double> *, 
//			   const MC_Truth *, const MC_Truth *);

int main(int argc, char* argv[])
{

  //  set the random number seed
  //  from /dev/random
  ifstream devrandom("/dev/random");
  int seed;
  devrandom.read( (char*)(&seed), sizeof(int) );
  if(seed < 0) seed = -seed;
  devrandom.close();
  gRandom->SetSeed(seed);
  
  char outfilename[100];
  strcpy(outfilename,argv[1]);
  //strcat(outfilename,"_NewProcessorTest.root");
  strcat(outfilename,"_PulseTags_DetEff.root");
  TFile *fout = new TFile(outfilename,"RECREATE");
  TTree * ot = new TTree("Events","PMT pulse info and MC truth");
  vector<double> * P_mergedPulseHeight;
  vector<double> * P_mergedPulseTime;
  vector<double> * P_northFraction;
  vector<string> * P_pulseTagCap;
  vector<string> * P_pulseTagDet;
  vector<string> * P_pulseTag;
  vector<string> * P_pulseTagCav;
  // for detection efficiency study
  int tm;  // tank muon, 0 or 1
  int lm; // lead muon, 0 or 1
  double lmEn;  // lead muon energy, greatest energy muon entering lead
  int ln;  // lead neutron, 0 or 1
  double lnEn;  // lead neutron energy, greatest energy neutron entering lead
  // v2 detection efficiency study
  // the tag of muon traversing lower panel
  int pm; // panel muon, 0 or 1
  // the tag of clipping muon
  int cm; //clipping muon, 0 or 1 
  
  ot->Branch("mergedPulseHeight",&P_mergedPulseHeight);
  ot->Branch("mergedPulseTime",&P_mergedPulseTime);
  ot->Branch("northFraction",&P_northFraction);
  ot->Branch("pulseTagCap",&P_pulseTagCap);
  ot->Branch("pulseTagDet",&P_pulseTagDet);
  ot->Branch("pulseTag",&P_pulseTag);
  ot->Branch("pulseTagCav",&P_pulseTagCav);
  // for detection efficiency study
  ot->Branch("TankMuon",&tm,"tm/I");
  ot->Branch("LeadMuon",&lm,"lm/I");
  ot->Branch("LeadMuonEnergy",&lmEn,"lm/D");
  ot->Branch("LeadNeutron",&ln,"ln/I");
  ot->Branch("LeadNeutronEnergy",&lnEn,"ln/D");
  // v2 detection efficiency study
  ot->Branch("PanelMuon",&pm,"pm/I");
  ot->Branch("ClipMuon",&cm,"cm/I");
  
  ifstream filelist;
  filelist.open(argv[1]);
  while(!filelist.eof()){
    char rootfilename[100];
    filelist >> rootfilename;
    TFile* fin = new TFile(rootfilename,"READ");
    cout << "Opening " << rootfilename << "\n";
    TTree* event = (TTree*)fin->Get("event");

    PMT_Hits pmtHits(27000000);
    pmtHits.LinkTree(event);
    MC_Truth cavern("cavern",278000);
    cavern.LinkTree(event);
    MC_Truth det("det",36000);
    det.LinkTree(event);
    MC_Truth lead("lead",36000);
    lead.LinkTree(event);
    CaptureInWater cap(1600);
    cap.LinkTree(event);
    Primary primaryMuon(1);
    primaryMuon.LinkTree(event);
    
    int numEntries = event->GetEntries();

    for (int n=0;n<numEntries;n++) {
      event->GetEntry(n);

      vector< vector<double> > pmtHitTimes;
      vector<double> northEastHitTimes;
      vector<double> northWestHitTimes;
      vector<double> southEastHitTimes;
      vector<double> southWestHitTimes;
      pmtHitTimes.push_back(northEastHitTimes);
      pmtHitTimes.push_back(northWestHitTimes);
      pmtHitTimes.push_back(southEastHitTimes);
      pmtHitTimes.push_back(southWestHitTimes);
      
      collect_hits(pmtHitTimes, pmtHits);
      
      vector< vector<double> > pmtPulseHeight;
      vector<double> northEastPulseHeight;
      vector<double> northWestPulseHeight;
      vector<double> southEastPulseHeight;
      vector<double> southWestPulseHeight;
      pmtPulseHeight.push_back(northEastPulseHeight);
      pmtPulseHeight.push_back(northWestPulseHeight);
      pmtPulseHeight.push_back(southEastPulseHeight);
      pmtPulseHeight.push_back(southWestPulseHeight);

      vector< vector<double> > pmtPulseTime;
      vector<double> northEastPulseTime;
      vector<double> northWestPulseTime;
      vector<double> southEastPulseTime;
      vector<double> southWestPulseTime;
      pmtPulseTime.push_back(northEastPulseTime);
      pmtPulseTime.push_back(northWestPulseTime);
      pmtPulseTime.push_back(southEastPulseTime);
      pmtPulseTime.push_back(southWestPulseTime);

      for (int p=0;p<4;p++) {
	if(pmtHitTimes[p].size()){
	  sort(pmtHitTimes[p].begin(),pmtHitTimes[p].end());
	  get_pmt_pulses(pmtPulseHeight[p], pmtPulseTime[p], pmtHitTimes[p],
			 (p+2)%4);
	  // outside this function, PMTs are counted as NE, NW, SE, SW
	  // convert to inside order--SE, SW, NE, NW
	}
      }
      
      vector<double> northSummedResponse;
      vector<double> southSummedResponse;
      vector<double> northCoincidentTime;
      vector<double> southCoincidentTime;
      
      vector< vector<double> >::iterator neph = pmtPulseHeight.begin();
      vector< vector<double> >::iterator nept = pmtPulseTime.begin();
 
      get_coin_pulses(northSummedResponse, northCoincidentTime, neph, nept);
      get_coin_pulses(southSummedResponse, southCoincidentTime, neph+2, nept+2);

      vector<double> mergedPulseHeight;
      vector<double> mergedPulseTime;
      vector<double> northFraction;
      
      merge_tank_pulses(mergedPulseHeight, mergedPulseTime,
			northFraction,
			northSummedResponse, southSummedResponse,
			northCoincidentTime, southCoincidentTime);

      //convert capture data to vector<double> format
      vector<double> captureTime;
      for (int i=0;i<cap.num;i++) {
	captureTime.push_back(cap.T[i]);
      }      
      sort(captureTime.begin(), captureTime.end());
      // reduce cavern, det and lead data by cerenkov threshold and sorted by time
      vector<MC_Track> rcav = preproc_mctruth(cavern);
      vector<MC_Track> rdet = preproc_mctruth(det);
      vector<MC_Track> rlead = convert_mctruth(lead);

      // Up to here, ready to do something
      vector<string> pulseTagCap;
      vector<string> pulseTagDet;
      vector<string> pulseTag;
      vector<string> pulseTagCav;
      /*************************************************
       * check the matching between pulses and captures
       *************************************************
       */
      vector<double>::iterator ip = mergedPulseTime.begin();
      vector<double>::iterator icap = captureTime.begin();
      while (ip != mergedPulseTime.end() && icap != captureTime.end()) {
	// make sure the iterators of pulse time and pulse tag moving together
	// namely pulseTagCap.push_back(tag) and ip++ taking place at same time
	
	if (match(*ip, *icap) == 0) {
	  // do operation for match
	  pulseTagCap.push_back("CapWater");
	  ip++;
	  icap++;
	}
	else if (match(*ip, *icap) < 0) {
	  // if pulse time is smaller
	  pulseTagCap.push_back("NonCapWater");
	  ip++;
	}
	else {
	  icap++;
	}	
      }
      while (ip != mergedPulseTime.end()) {
	pulseTagCap.push_back("NonCapWater");
	ip++;
      }
      /*************************************************/
      
      /******************************************************
       * check the matching between pulses and det particles
       ******************************************************
       */
      ip = mergedPulseTime.begin();
      vector<MC_Track>::iterator idet = rdet.begin();
      vector<MC_Track>::iterator icav = rcav.begin();
      while (ip != mergedPulseTime.end() && idet != rdet.end()) {
	if (match(*ip, (*idet).T) == 0) {
	  if ((*idet).PDG == 13 || (*idet).PDG == -13) {
	    (*idet).TrigPulse = distance(mergedPulseTime.begin(),ip);
	    pulseTagDet.push_back("muon");
	    ip++;
	  }
	  else if ((*idet).PDG == 22) {
	    (*idet).TrigPulse = distance(mergedPulseTime.begin(),ip);
	    pulseTagDet.push_back("gamma");
	    ip++;
	  }
	  else if ((*idet).PDG == 11 || (*idet).PDG == -11) {
	    (*idet).TrigPulse = distance(mergedPulseTime.begin(),ip);
	    pulseTagDet.push_back("e");
	    ip++;
	  }
	  else if ((*idet).PDG == 2212) {
	    (*idet).TrigPulse = distance(mergedPulseTime.begin(),ip);
	    pulseTagDet.push_back("proton");
	    ip++;
	  }
	  //need to check the pdg code histogram to include all interesting particles
	  idet++; 
	  // det particle either the source of pulse, or not interesting (nu, n, etc.)
	}
	else if ( match(*ip, (*idet).T) < 0) {
	  pulseTagDet.push_back("NonDet");
	  ip++;
	}
	else {
	  idet++;
	}	
      }
      while (ip != mergedPulseTime.end()) {
	pulseTagDet.push_back("NonDet");
	ip++;
      }
      /******************************************************/
     
      /* integrate info from pulseTagCap and pulseTagDet
       */
      if (pulseTagCap.size() != pulseTagDet.size())
	cerr << "Error: Tag sizes does not match.\n";
      else {
	for (int i=0;i<pulseTagCap.size();i++) {
	  if (pulseTagCap[i] == "CapWater")
	    pulseTag.push_back("CapWater");
	  else if (pulseTagDet[i] == "NonDet")
	    pulseTag.push_back("Unknown");
	  else 
	    pulseTag.push_back(pulseTagDet[i]);
	}
      }
      /*************************************************/

       //cout<< "size of pulseTag: " << pulseTag.size() <<"\n";
      if (pulseTag.size() != 0) {
	for (int i=0;i<pulseTag.size();i++) {
	  pulseTagCav.push_back("Null");
	}
	for (int i=0;i<rdet.size();i++) {
	  if (rdet[i].TrigPulse > -1) {
	    for (int j=0;j<rcav.size();j++) {
	      if(rdet[i].TID == rcav[j].TID && rdet[i].T > rcav[j].T) {
		rcav[j].TrigPulse = rdet[i].TrigPulse;
		pulseTagCav[rcav[j].TrigPulse] = "FromCav";
	      }
	      
	    }
	  }
	}
      }

      //***************************************
      // tags for fast-neutron efficiency study
      //***************************************
      int TankMuon = 0;
      int LeadMuon = 0;
      double LeadMuonEnergy = 0.;
      int LeadNeutron = 0;
      double LeadNeutronEnergy = 0.;
      for (int i = 0; i < rdet.size(); i++){
	if(rdet[i].PDG == 13 || rdet[i].PDG == -13)
	  TankMuon = 1;
      }
      for (int i = 0; i < rlead.size(); i++){
	if(rlead[i].PDG == 13 || rlead[i].PDG == -13){
	  LeadMuon = 1;
	  LeadMuonEnergy = rlead[i].E > LeadMuonEnergy ? rlead[i].E : LeadMuonEnergy;
	}
	if(rlead[i].PDG == 2112){
	  LeadNeutron = 1;
	  LeadNeutronEnergy = rlead[i].E > LeadNeutronEnergy ? rlead[i].E : LeadNeutronEnergy;
	}
      }
      
      //*** PanelMuon Tag****
      int PanelMuon = 0;
      // Lower Panel Geometry
      double xLeadPos = 5940; // mm
      double yLeadPos = 1984;
      double zLeadPos = -3487;
      Rectangle pan37_1(xLeadPos - 50, yLeadPos - 830, zLeadPos - 650,
			2580, 1480);
      Rectangle pan37_2(xLeadPos - 110, yLeadPos + 670, zLeadPos - 650,
			2540, 1280);
      Rectangle pan38(xLeadPos - 120, yLeadPos, zLeadPos - 1320,
			2940, 3680);
      // Primary Muon
      double p0 = sqrt(primaryMuon.Px[0]*primaryMuon.Px[0]+
		       primaryMuon.Py[0]*primaryMuon.Py[0]+
		       primaryMuon.Pz[0]*primaryMuon.Pz[0]);
      TrackLine muon_track(primaryMuon.X[0],primaryMuon.Y[0],primaryMuon.Z[0],
			   primaryMuon.Px[0]/p0,primaryMuon.Py[0]/p0,
			   primaryMuon.Pz[0]/p0);
      
      bool muTopPanel = Traverse(muon_track,pan37_1) || Traverse(muon_track,pan37_2);
      bool muBottomPanel = Traverse(muon_track,pan38);

      if(muTopPanel && muBottomPanel) PanelMuon = 1;

      //*** ClipMuon Tag ***
      int ClipMuon = 0;
      
      double maxph12;     
      if(mergedPulseHeight.size()>=2)
	maxph12 = mergedPulseHeight[0] > mergedPulseHeight[1] 
	  ? mergedPulseHeight[0] : mergedPulseHeight[1];
      else if(mergedPulseHeight.size()==1)
	maxph12 = mergedPulseHeight[0];
      else
	maxph12 = -999;

      if(TankMuon == 1 && maxph12 > 0 && maxph12 < 200) ClipMuon = 1;

      /*
      // print info for debugging
      if (mergedPulseTime.size() > 5 && mergedPulseTime.size() <11) {
	cout<<"\n";
	cout<<mergedPulseTime.size()<<" pulses.\n";
	for (int i=0;i<mergedPulseTime.size();i++) {
	  cout<< mergedPulseTime[i] *1000. <<", " 
	      << pulseTagCap[i] <<", "
	      << pulseTagDet[i] <<", "
	      << pulseTag[i] <<", "
	      << pulseTagCav[i] <<"\n";
	}
	cout<<captureTime.size()<<" captures.\n";
	for (int i=0;i<captureTime.size();i++) {
	  cout<< captureTime[i] <<"\n";
	}
	cout<<rdet.size()<<" reduced dets.\n";
	for (int i=0;i<rdet.size();i++) {
	  cout<< rdet[i].T  <<", " 
	      << rdet[i].PDG <<", "
	      << rdet[i].TID <<", "
	      << rdet[i].PID 
	      <<"\n";
	}
	cout<<rcav.size()<<" reduced cavs.\n";
	for (int i=0;i<rcav.size();i++) {
	  cout<< rcav[i].T  <<", " 
	      << rcav[i].PDG <<", "
	      << rcav[i].TID <<", "
	      << rcav[i].PID 
	      <<"\n";
	}
      }
      */

      // fill tree
      P_mergedPulseHeight = &mergedPulseHeight;
      P_mergedPulseTime = &mergedPulseTime;
      P_northFraction = &northFraction;
      P_pulseTagCap = &pulseTagCap;
      P_pulseTagDet = &pulseTagDet;
      P_pulseTag = &pulseTag;
      P_pulseTagCav = &pulseTagCav;
      //for detection efficiency study
      tm = TankMuon;
      lm =  LeadMuon;
      lmEn = LeadMuonEnergy;
      ln = LeadNeutron;
      lnEn = LeadNeutronEnergy;
      pm = PanelMuon;
      cm = ClipMuon;
      ot->Fill();
    } // Entries
    
    fin->Close();
    delete fin;
  }// root files
  fout->Write();
  fout->Close();

  /*
  cout<< "Cerenkov threshold for gamma is: " << get_cer_thr(22) <<"\n";
  cout<< "Cerenkov threshold for electron is: " << get_cer_thr(11) <<"\n";
  cout<< "Cerenkov threshold for muon is: " << get_cer_thr(13) <<"\n";
  cout<< "Cerenkov threshold for proton is: " << get_cer_thr(2212) <<"\n";
  cout<< "Cerenkov threshold for pi is: " << get_cer_thr(211) <<"\n";
  */

  return 0;
}

PMT_Hits::PMT_Hits(int len) {
  E = new double[len];
  T = new double[len];
  X = new double[len];
  Y = new double[len];
}

PMT_Hits::~PMT_Hits() {
  delete [] E;
  delete [] T;
  delete [] X;
  delete [] Y;
}

void PMT_Hits::LinkTree (TTree *event) {
  event->SetBranchAddress("numPmtHits",&num);
  event->SetBranchAddress("photonHitEnergy",E);
  event->SetBranchAddress("photonHitT",T);
  event->SetBranchAddress("photonHitX",X);
  event->SetBranchAddress("photonHitY",Y);
}

MC_Truth::MC_Truth(const char *NAME, int len) {
  strcpy(name,NAME);
  strcpy(Name,NAME);
  Name[0] = (char)toupper(NAME[0]);
  
  X = new double[len];
  Y = new double[len];
  Z = new double[len];
  T = new double[len];
  Px = new double[len];
  Py = new double[len];
  Pz = new double[len];
  E = new double[len];
  TID = new int[len];
  PID = new int[len];
  PDG = new int[len];
}

MC_Truth::~MC_Truth() {
  delete [] X;
  delete [] Y;
  delete [] Z;
  delete [] T;
  delete [] Px;
  delete [] Py;
  delete [] Pz;
  delete [] E;
  delete [] TID;
  delete [] PID;
  delete [] PDG;
}

void MC_Truth::LinkTree (TTree *event) {
  event->SetBranchAddress(Form("num%sParticles",Name),&num);
  event->SetBranchAddress(Form("%sParticleX",name),X);
  event->SetBranchAddress(Form("%sParticleY",name),Y);
  event->SetBranchAddress(Form("%sParticleZ",name),Z);
  event->SetBranchAddress(Form("%sParticleT",name),T);
  event->SetBranchAddress(Form("%sParticlePx",name),Px);
  event->SetBranchAddress(Form("%sParticlePy",name),Py);
  event->SetBranchAddress(Form("%sParticlePz",name),Pz);
  event->SetBranchAddress(Form("%sParticleE",name),E);
  event->SetBranchAddress(Form("%sParticleTID",name),TID);
  event->SetBranchAddress(Form("%sParticlePID",name),PID);
  event->SetBranchAddress(Form("%sParticlePDG",name),PDG);
}

CaptureInWater::CaptureInWater(int len)
{
  strcpy(name,"capture");
  strcpy(Name,"Capture");
  X = new double[len];
  Y = new double[len];
  Z = new double[len];
  T = new double[len];
}

CaptureInWater::~CaptureInWater( )
{
  delete [] X;
  delete [] Y;
  delete [] Z;
  delete [] T;
}

void CaptureInWater::LinkTree (TTree *event) {
  event->SetBranchAddress(Form("num%ss",Name),&num);
  event->SetBranchAddress(Form("%sX",name),X);
  event->SetBranchAddress(Form("%sY",name),Y);
  event->SetBranchAddress(Form("%sZ",name),Z);
  event->SetBranchAddress(Form("%sT",name),T);
}

Primary::Primary(int len) {
  strcpy(name,"primary");
  strcpy(Name,"Primary");
  
  E = new double[len];
  Px = new double[len];
  Py = new double[len];
  Pz = new double[len];
  X = new double[len];
  Y = new double[len];
  Z = new double[len];
}

Primary::~Primary() {
  delete [] E;
  delete [] Px;
  delete [] Py;
  delete [] Pz;
  delete [] X;
  delete [] Y;
  delete [] Z;
}

void Primary::LinkTree (TTree *event) {
  event->SetBranchAddress(Form("num%sParticles",Name),&num);
  event->SetBranchAddress(Form("%sParticleEnergy",name),E);
  event->SetBranchAddress(Form("%sParticlePx",name),Px);
  event->SetBranchAddress(Form("%sParticlePy",name),Py);
  event->SetBranchAddress(Form("%sParticlePz",name),Pz);
  event->SetBranchAddress(Form("%sParticleX",name),X);
  event->SetBranchAddress(Form("%sParticleY",name),Y);
  event->SetBranchAddress(Form("%sParticleZ",name),Z);
}

vector<MC_Track> convert_mctruth(const MC_Truth &mcdata)
{
  vector<MC_Track> track_vec;
  for (int i=0;i<mcdata.num;i++){
    MC_Track trk = {mcdata.X[i], mcdata.Y[i], mcdata.Z[i], mcdata.T[i],
		    mcdata.Px[i], mcdata.Py[i], mcdata.Pz[i], mcdata.E[i],
		    mcdata.TID[i], mcdata.PID[i], mcdata.PDG[i], -1};
    track_vec.push_back(trk);
  }
  return track_vec;
}

vector<MC_Track> preproc_mctruth(const MC_Truth &mcdata)
{
  /* 1.Convert MC_Truth data to MC_Track vector
   * 2.Sort tracks by momentum descendingly
   * 3.Erase the tracks having energy smaller than Cerenkov threshold
   * 4.Sort tracks by time, and return reduced track vector data 
  */
  vector<MC_Track> tracks = convert_mctruth(mcdata);
  sort(tracks.begin(),tracks.end(),comp_energy);
  cut_cer_thr(tracks);
  sort(tracks.begin(),tracks.end(),comp_time);

  return tracks;
}

double get_cer_thr (int pdg) {
  // only for charged particles and gamma
  if (abs(pdg) == 2112 || abs(pdg) ==310 || abs(pdg) == 130 ||
      abs(pdg) == 12 || abs(pdg) ==14) {
    cout<< "Please only pass charged particle or gamma to get_cer_thr().\n";
    return 1.E9;
  }
  else if (pdg == 22)
    pdg = 11;
  
  double m = get_mass(pdg);
  double n = 1.333; // water refractive index
  return m*n/sqrt(n*n-1) - m;
}

void cut_cer_thr(vector<MC_Track> &mctrks)
{
// This function erases paticle tracks under Cerenkov threshold
// Pass-in a descendingly energy sorted vector<MC_Track>
// this may help improve the efficiency of erase()
  vector<MC_Track>::iterator aTrack=mctrks.begin();
  for ( ; aTrack!=mctrks.end(); ++aTrack){
    int pdg = (*aTrack).PDG;
    if (abs(pdg)==2112 || abs(pdg)==310 || abs(pdg)==130 || 
	abs(pdg)==12 || abs(pdg)==14 ) {
      (*aTrack).E = -1.;
      continue;
    } 
    double E = (*aTrack).E;
    double thr = get_cer_thr(pdg);
    if (E < thr )
      (*aTrack).E = -1.;
  }
  sort(mctrks.begin(),mctrks.end(),comp_energy);
  for (aTrack=mctrks.begin(); aTrack!=mctrks.end(); ++aTrack)
    if ( (*aTrack).E < 0 ) break;
  mctrks.erase(aTrack,mctrks.end());
}

void collect_hits(vector< vector<double> > &pmtHitTimes, const PMT_Hits &pmtHits)
{
  /* 
   * [0] northEast
   * [1] northWest
   * [2] southEast
   * [3] sorthWest
   */
  double yOffset = 1984.;
  double xOffset = 5940.;
  
  for(int h=0; h<pmtHits.num; h++){
    if(PmtHitsToPhotoElectrons(1240./pmtHits.E[h])){
      if((pmtHits.Y[h]-yOffset)>0.) {	
	if((pmtHits.X[h]-xOffset)>0.)
	  pmtHitTimes[0].push_back(pmtHits.T[h]*0.001);
	else 
	  pmtHitTimes[1].push_back(pmtHits.T[h]*0.001);
      }
      else {					
	if((pmtHits.X[h]-xOffset)>0.)
	  pmtHitTimes[2].push_back(pmtHits.T[h]*0.001);
	else 
	  pmtHitTimes[3].push_back(pmtHits.T[h]*0.001);
      }
    }
  }

}

void get_pmt_pulses(vector<double> &pmtPulseHeight, vector<double> &pmtPulseTime,
		    const vector<double> &pmtHitTimes , const int iPMT)
{     
  /* intend to ensure the pmtHitTimes non-empty and sort it outside this function
   * if(pmtHitTimes.size())
   * sort(pmtHitTimes.begin(), pmtHitTimes.end())
   *********************************************
   * now, inside this function, assume *pmtHitTimes is sorted
   */
  // low-gain settings from Cf252 calibration fit
  //double scaling[4] = {1.516, 1.337, 1.637, 1.532}; // mV/PE
  //double s0 = 2.3; // mV
  //double s1[4] = {2.886, 2.663, 1.768, 2.181}; // sqrt(mV)

  //high-gain settings from Cf252 calibration fit
  //double scaling[4] = {2.279, 1.992, 1.951, 2.047}; // mV/PE
  //double s0 = 2.3; // mV
  //double s1[4] = {1.703, 1.667, 1.95, 1.753}; // sqrt(mV)

  double firstOrderScale = 2.5;
  double secondOrderScale = 0.9;

  double lastPulseStartTime = pmtHitTimes[0];
  double lastPulseHeight = 0.;
  double lastPulseTime = 0.;
  for(int p=0; p<pmtHitTimes.size(); p++){
    if (fabs(pmtHitTimes[p] - lastPulseStartTime) < pulseWidth) {
      lastPulseHeight++;
      lastPulseTime += pmtHitTimes[p];
      
      if (p == pmtHitTimes.size()-1) {
	//double pulseHeight = lastPulseHeight + gRandom->Rndm() - 0.5;
	double pulseHeight = lastPulseHeight;
	// old way:
	pulseHeight *= firstOrderScale;
	pulseHeight = gRandom->Gaus(pulseHeight,secondOrderScale*sqrt(pulseHeight));
	// new way:
	//pulseHeight *= scaling[iPMT];
	//pulseHeight = gRandom->Gaus(pulseHeight,
	//		    sqrt(s0*s0+s1[iPMT]*s1[iPMT]*pulseHeight));

	pmtPulseHeight.push_back(pulseHeight);
	pmtPulseTime.push_back(lastPulseTime/lastPulseHeight);
      }
    }
    else {
      //double pulseHeight = lastPulseHeight + gRandom->Rndm() - 0.5;
      
      double pulseHeight = lastPulseHeight;
      pulseHeight *= firstOrderScale;
      pulseHeight = gRandom->Gaus(pulseHeight,secondOrderScale*sqrt(pulseHeight));
      //pulseHeight *= scaling[iPMT];
      //pulseHeight = gRandom->Gaus(pulseHeight,
      //			  sqrt(s0*s0+s1[iPMT]*s1[iPMT]*pulseHeight));

      pmtPulseHeight.push_back(pulseHeight);
      pmtPulseTime.push_back(lastPulseTime/lastPulseHeight);
      
      lastPulseHeight = 1;
      lastPulseTime = pmtHitTimes[p];
      lastPulseStartTime = pmtHitTimes[p];
    }
  }
  
}

void get_coin_pulses(vector<double> &summed, vector<double> &coin,
		     const vector< vector<double> >::iterator pph,
		     const vector< vector<double> >::iterator ppt)
{
  /*
   * summed: summed pulse heigh over 2 PMTs
   * coin: coincident pulse time between 2 PMTs
   * pph: pmt pulse height
   * ppt: pmt pulse time
   * pph and ppt are iterators starting from the first pmt of a tank
  */
  double threshold = 10.; // mV
  
  int numPulsesPMT1 = (*pph).size();
  int numPulsesPMT2 = (*(pph+1)).size();
  if (numPulsesPMT1 + numPulsesPMT2) {
    for (int t1=0;t1<numPulsesPMT1;t1++) {
      for (int t2=0;t2<numPulsesPMT2;t2++) {
	double eastPH = (*pph)[t1];
	double westPH = (*(pph+1))[t2];
	if ( eastPH>threshold && westPH>threshold
	     && fabs((*ppt)[t1] - (*(ppt+1))[t2])<pulseWidth ) {
	  summed.push_back(eastPH + westPH);
	  coin.push_back( ((*ppt)[t1] + (*(ppt+1))[t2])/2. );
	  break;
	}	
      }
    }
  }

}

void merge_tank_pulses (vector<double> &mergedPH, vector<double> &mergedPT, 
			vector<double> &northFrac,
			const vector<double> &nPH, const vector<double> &sPH,
			const vector<double> &nPT, const vector<double> &sPT)
{
  vector<double>::const_iterator itnpt = nPT.begin();
  vector<double>::const_iterator itspt = sPT.begin();
  vector<double>::const_iterator itnph = nPH.begin();
  vector<double>::const_iterator itsph = sPH.begin();
  
  while(itnpt!=nPT.end() && itspt!=sPT.end()) {
    if ( fabs((*itnpt) - (*itspt)) < pulseWidth ) {
      double NPH = *itnph++;
      double SPH = *itsph++;
      mergedPH.push_back(NPH + SPH);
      mergedPT.push_back( ((*itnpt++) + (*itspt++))/2. );
      northFrac.push_back(NPH / (NPH + SPH));
    }
    else if( (*itnpt) < (*itspt) ) {
      mergedPH.push_back( *itnph++ );
      mergedPT.push_back( *itnpt++ );
      northFrac.push_back( 1. );
    }
    else {
      mergedPH.push_back( *itsph++ );
      mergedPT.push_back( *itspt++ );
      northFrac.push_back( 0. );
    }   
  }

  while (itnpt!=nPT.end()) {
    mergedPH.push_back( *itnph++ );
    mergedPT.push_back( *itnpt++ );
    northFrac.push_back( 1. );
  }
  while (itspt!=sPT.end()) {
    mergedPH.push_back( *itsph++ );
    mergedPT.push_back( *itspt++ );
    northFrac.push_back( 0. );
  }

}

int match (double ptime, double mctime) { 
  // ptime in us and mctime in ns
  double dt = ptime*1000. - mctime;
  if (dt>0. && dt<300)
    return 0;
  else if (dt < 0.)
    return -1;
  else
    return 1;
}

/*
vector<string> GetPulseTag(const vector<double> *pt, const vector<double> *ct, 
			   const MC_Truth *det, const MC_Truth *cav)
{

}
*/

double get_mass (int pdg) {
  // return mass in MeV
  double mass;
  int abs_pdg = abs(pdg);
  switch (abs_pdg) {
  case 11:
    //electron
    mass = 0.511;
    break;
  case 13:
    //muon
    mass = 105.658;
    break;
  case 2212:
    //proton
    mass = 938.272;
    break;
  case 321:
    //K
    mass = 493.677;
    break;
  case 211:
    // pi
    mass = 139.570;
    break;
  case 22:
    //gamma
    mass = 0.;
    break;
  case 2112:
    //neutron
    mass = 939.565;
    break;
  case 310:
    //K_S^0
    mass = 497.614;
    break;
  case 130:
    //K_L^0
    mass = 497.614;
    break;
  case 12:
    //neutrino
    mass = 0.;
    break;
  case 14:
    //neutrino
    mass = 0.;
    break;
  default:
    cerr<<"Unkown particle: "<< pdg <<"!\n";
  }
  return mass;
}

bool PmtHitsToPhotoElectrons(double photonWavelength){
  double wavelength[] = {609.64,599.61,589.63,579.71,569.59,
			 559.57,549.65,539.60,529.69,519.70,
			 509.66,499.80,498.59,496.78,494.79,
			 492.89,491.25,489.95,488.40,486.76,
			 485.46,484.17,482.70,481.15,479.51,
			 477.87,476.07,474.43,472.97,471.51,
			 470.13,468.58,467.38,466.09,464.63,
			 463.35,462.50,461.39,460.45,459.60,
			 458.75,458.16,457.31,456.20,455.44,
			 453.98,452.25,450.08,447.90,445.72,
			 444.04,442.46,441.14,439.81,438.30,
			 437.06,435.90,434.66,433.51,432.52,
			 431.37,430.30,429.32,428.42,427.62,
			 426.55,425.65,424.85,424.04,423.15,
			 422.25,421.53,420.64,419.74,418.85,
			 417.78,416.63,415.65,414.58,413.25,
			 412.26,411.03,409.52,408.10,406.86,
			 405.89,403.95,401.75,399.55,397.36,
			 395.17,392.97,390.79,388.63,385.08,
			 381.53,377.94,374.27,370.77,367.40,
			 364.34,360.58,356.46,352.11,347.74,
			 343.37,340.29,336.64,332.38,328.01,
			 323.71,320.57,318.56,317.17,316.26,
			 315.21,314.25,313.34,312.33,310.46,
			 306.65,302.30,298.08,295.12,292.72,
			 290.63,288.49,286.53,284.61,282.64,
			 280.42,278.37,275.93,271.77,269.15,
			 268.63,268.19,267.98,267.93,267.63,
			 267.24,267.11,266.72,266.32,266.02,
			 265.76,265.06,264.67,264.28,263.89,
			 263.50,263.10,262.71,262.32,261.62,
			 260.97,260.36,259.71,258.97,258.45,
			 257.62,256.83,256.22,255.53,254.66,
			 253.78,252.83,251.61,250.34,247.27,
			 243.20,241.11,239.67,238.36,236.88,
			 235.52,234.22,232.65,229.06,225.74,
			 223.17,220.99};
  
  double quantumEff[] = {0.00751,0.01165,0.01634,0.02339,0.02971,
			 0.04087,0.04988,0.06211,0.07141,0.08714,
			 0.09438,0.10222,0.11050,0.11070,0.11090,
			 0.11300,0.11500,0.11700,0.11988,0.12229,
			 0.12476,0.12476,0.12726,0.12726,0.12983,
			 0.12983,0.13244,0.13244,0.13511,0.13783,
			 0.13783,0.14060,0.14060,0.14343,0.14343,
			 0.14632,0.14632,0.14926,0.14926,0.15227,
			 0.15227,0.15227,0.15227,0.15533,0.15533,
			 0.15846,0.16165,0.16165,0.16490,0.16822,
			 0.16822,0.16822,0.17161,0.17161,0.17506,
			 0.17506,0.17858,0.17858,0.17858,0.17858,
			 0.17858,0.18218,0.18218,0.18218,0.18218,
			 0.18584,0.18584,0.18584,0.18584,0.18584,
			 0.18584,0.18958,0.18958,0.18958,0.18958,
			 0.18958,0.18958,0.19340,0.19340,0.19340,
			 0.19340,0.19340,0.19729,0.19729,0.19729,
			 0.19729,0.19729,0.19729,0.19729,0.19729,
			 0.19729,0.19729,0.19729,0.19729,0.19729,
			 0.19729,0.19729,0.19729,0.19729,0.19729,
			 0.19729,0.19729,0.19340,0.19340,0.18958,
			 0.18958,0.18584,0.18218,0.17506,0.17161,
			 0.16490,0.15533,0.15227,0.14926,0.14926,
			 0.14632,0.14632,0.14060,0.14060,0.13783,
			 0.12229,0.11293,0.08891,0.07000,0.05511,
			 0.04339,0.04339,0.03416,0.02690,0.02117,
			 0.01667,0.01313,0.00814,0.00397,0.00313,
			 0.00000,0.00000,0.00000,0.00000,0.00000,
			 0.00000,0.00000,0.00000,0.00000,0.00000,
			 0.00000,0.00000,0.00000,0.00000,0.00000,
			 0.00000,0.00000,0.00000,0.00000,0.00000,
			 0.00000,0.00000,0.00000,0.00000,0.00000,
			 0.00000,0.00000,0.00000,0.00000,0.00000,
			 0.00000,0.00000,0.00000,0.00000,0.00000,
			 0.00000,0.00000,0.00000,0.00000,0.00000,
			 0.00000,0.00000,0.00000,0.00000,0.00000,
			 0.00000,0.00000};
  
  
  int numEntries = 182;
  double rand = gRandom->Rndm();
  
  //	find closest wavelength entry
  //	use that as the quantum efficiency
  double effDiff = 1000000.;
  double radDiff = 1000000.;
  int quanIndex = 0;
  int radIndex = 0;
  for(int i=0; i<numEntries; i++){
    if(fabs(wavelength[i]-photonWavelength) < effDiff){ 
      effDiff = fabs(wavelength[i]-photonWavelength);
      quanIndex = i;
    }
  }
  
  //	if the random number is less	
  //	than the efficiency, we got a hit
  if(rand > quantumEff[quanIndex]) 
    return false;
  else
    return true;
  
}

// geometric tools for tagging muon hitting lower panels
TrackLine::TrackLine(double posX, double posY, double posZ, double dirX, double dirY, double dirZ){
    x0 = posX;
    y0 = posY;
    z0 = posZ;
    dx = dirX;
    dy = dirY;
    dz = dirZ;
}

TVector3 TrackLine::GetPoint(double coord,char coordName) {
    double x,y,z;
    double ratio;
    if (coordName=='x') {
        x = coord;
        ratio = (x - x0)/dx;
        y = ratio*dy + y0;
        z = ratio*dz + z0;
    } else if(coordName=='y') {
        y = coord;
        ratio = (y - y0)/dy;
        x = ratio*dx + x0;
        z = ratio*dz + z0;
    } else if(coordName=='z'){
        z = coord;
        ratio = (z - z0)/dz;
        x = ratio*dx + x0;
        y = ratio*dy + y0;
    } else {
        cout<<"Please provide coordinate name. Candidates: 'x', 'y', or 'z'.\n";
        return TVector3(x0,y0,z0);
    }
     
    return TVector3(x,y,z);
}

Rectangle::Rectangle(){ }
Rectangle::Rectangle(double cent1,double cent2,double cent3, double w, double h){
  center[0] = cent1; //x
  center[1] = cent2; //y
  center[2] = cent3; //z
  side[0] = w; // x side
  side[1] = h; // y side
}

bool Rectangle::Inside(double a, double b){
    if ( (Bottom(0)<=a && a<=Top(0)) && (Bottom(1)<=b && b<=Top(1)) ) {
        return true;
    } else{
        return false;
    }
}

bool Traverse(TrackLine muon, Rectangle panel){
  
  TVector3 cross_point = muon.GetPoint(panel.Z(),'z');
  return panel.Inside(cross_point.X(),cross_point.Y());
}

/*
// To show PMT hit times, pulse heights, and pulse times, 
// pick up one entry (event), use the following code:

  TH1F h1("h1","pmt hits;time (#mus);",10000,0,10);
  TH1F h2("h2","pulses;time (#mus);",10000,0,10);

// in entry loop, do:
      for (int i=0;i<pmtHitTimes[0].size();i++) {
	h1.Fill(pmtHitTimes[0][i]);
      }

      for (int i=0;i<pmtPulseTime[0].size();i++) {
	int num_ns = (int)(pmtPulseTime[0][i]*1000.);

	h2.SetBinContent(num_ns,pmtPulseHeight[0][i]);
      }

// near the end of main, do:

  TFile *fout = new TFile("PulseAndHitTest.root","RECREATE");
  TCanvas c1("c1","c1");
  h1.SetLineColor(4);
  h2.SetLineColor(2);
  h1.Draw();
  h2.Draw("same"); 
  h1.Write();
  h2.Write();
  c1.Write();
  fout->Close();
  delete fout;
*/
