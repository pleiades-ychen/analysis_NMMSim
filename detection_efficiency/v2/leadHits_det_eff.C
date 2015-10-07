{
  if (!TClass::GetDict("NMMEvent")) {
    gROOT->ProcessLine(".L NMMEvent.cc+");
  }

  double thr = 50.;

  int totLeadHits = 0;
  int numDet = 0;
  int numPassMuonCut = 0;
  TH1F *hMLhit = new TH1F("hMLhit",";Primary Energy (MeV);Events / 10 MeV"
			  ,100,0,1000);
  TH1F *hNLhit = new TH1F("hNLhit",";Primary Energy (MeV);Events / 10 MeV"
			  ,100,0,1000);

  TH1F *hMLdet = new TH1F("hMLdet",";Primary Energy (MeV);Events / 10 MeV"
			  ,100,0,1000);
  TH1F *hNLdet = new TH1F("hNLdet",";Primary Energy (MeV);Events / 10 MeV"
			  ,100,0,1000);

  TFile *infile = new TFile("melindas_neutrons_2ndProc_DetEff.root","READ");
  TTree *event = (TTree*)infile->Get("Event");
  NMMEvent *nmmEvent;
  event->SetBranchAddress("NMMEvent",&nmmEvent);
  int NumEvents = event->GetEntries();
  for (int i = 0; i < NumEvents; i++){
    event->GetEntry(i);

    int N = nmmEvent->N;
    vector<double> ph = nmmEvent->pulseHeights;
    vector<double> pt = nmmEvent->pulseTimes;
    string type = nmmEvent->type;
    int M = nmmEvent->M;
    int TankMuon = nmmEvent-> TankMuon;
    int LeadMuon = nmmEvent->LeadMuon;
    double LeadMuonEnergy = nmmEvent->LeadMuonEnergy;
    int LeadNeutron = nmmEvent->LeadNeutron;
    double LeadNeutronEnergy = nmmEvent->LeadNeutronEnergy;

    if(TankMuon==0 && (LeadMuon==1 || LeadNeutron==1)){
      // for drawing the spectra of lead hits
      // not exactly the "denominator" events 
      if(LeadMuon==1) hMLhit->Fill(LeadMuonEnergy);
      if(LeadNeutron==1) hNLhit->Fill(LeadNeutronEnergy);

      // "denominator" events, add threshold
      if(LeadMuonEnergy>thr || LeadNeutronEnergy>thr){ 
	totLeadHits++;
           
	if (N >= 5) { // detection
	  numDet++;
	  if(ph[0]<200 && ph[1]<200){
	    numPassMuonCut++;
	    if(LeadMuon==1) 
	      hMLdet->Fill(LeadMuonEnergy);
	    if(LeadNeutron==1) 
	      hNLdet->Fill(LeadNeutronEnergy); 
	  }

	}
      }
    }
  }

  // efficiency for N>=5 detection
  double eff = (double)numDet/totLeadHits;
  double err = eff * sqrt(1./numDet + 1./totLeadHits);
  cout<< "Eff = " << eff <<" +/- "<< err <<"\n";
  // efficiency for N>=5 && pass muon cut detection
  double eff2 = (double)numPassMuonCut/totLeadHits;
  double err2 = eff * sqrt(1./numPassMuonCut + 1./totLeadHits);
  cout<< "Eff2 = " << eff2 <<" +/- "<< err2 <<"\n";
  
  hMLhit->SetLineColor(46);
  hMLdet->SetLineColor(2);
  hNLhit->SetLineColor(9);
  hNLdet->SetLineColor(4);

  hMLhit->SetLineWidth(2);
  hMLdet->SetLineWidth(2);
  hNLhit->SetLineWidth(2);
  hNLdet->SetLineWidth(2);

  new TCanvas;
  hNLhit->Draw();
  hNLdet->Draw("same");
  gPad->SetLogy();
  new TCanvas;
  hMLhit->Draw();  
  hMLdet->Draw("same");
  gPad->SetLogy();
}
