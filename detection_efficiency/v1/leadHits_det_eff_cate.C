{
  if (!TClass::GetDict("NMMEvent")) {
    gROOT->ProcessLine(".L NMMEvent.cc+");
  }

  double thr = 50.;
  double MuonCut = 200.;

  int totLeadHits = 0;
  int numTrg = 0;
  int numDet = 0;
  
  TH1F *hhit = new TH1F("hhit","Muon hitting both tank(s) and lead (category 4)"
			";Primary Energy (MeV);Events / 10 GeV",
			100,0,1000000);
  //High-energy neutron above 50 MeV
  //Muon hitting lead
  //Muon hitting both tank(s) and lead (category 3)
  TH1F *htrg = new TH1F("htrg",";Primary Energy (MeV);Events / 10 GeV",
			100,0,1000000);
  TH1F *hdet = new TH1F("hdet",";Primary Energy (MeV);Events / 10 GeV",
			100,0,1000000);

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

    //if(TankMuon==0 &&  LeadNeutron==1 && LeadNeutronEnergy > thr){ // category 1
    //if(TankMuon==0 && LeadMuon==1 && LeadNeutron==0){ // category 2
    //if(TankMuon==1 && LeadMuon==1 && !(LeadNeutron==1 && LeadNeutronEnergy > thr)){ // category 3
    if(TankMuon==1 && LeadMuon==1 && LeadNeutron==1 && LeadNeutronEnergy > thr){ // category 4
      totLeadHits++;
      hhit->Fill(LeadMuonEnergy);   // LeadNeutronEnergy  
      if (N >= 5) { // triggered
	numTrg++;
	htrg->Fill(LeadMuonEnergy); // LeadNeutronEnergy  
	if(ph[0] < MuonCut && ph[1] < MuonCut){ // detected
	  numDet++;
	  hdet->Fill(LeadMuonEnergy); // LeadNeutronEnergy
	}
      }
    }
      
  }

  // efficiency for N>=5 detection
  double eff = (double)numTrg/totLeadHits;
  double err = eff * sqrt(1./numTrg + 1./totLeadHits);
  cout<< "Eff = " << eff <<" +/- "<< err <<"\n";
  // efficiency for N>=5 && pass muon cut detection
  double eff2 = (double)numDet/totLeadHits;
  double err2 = eff * sqrt(1./numDet + 1./totLeadHits);
  cout<< "Eff2 = " << eff2 <<" +/- "<< err2 <<"\n";
  
  hhit->SetLineColor(1);
  htrg->SetLineColor(6);
  hdet->SetLineColor(4);
  //hdet->SetLineColor(2);

  hhit->SetLineWidth(2);
  htrg->SetLineWidth(2);
  hdet->SetLineWidth(2);
  
  new TCanvas;
  hhit->Draw();
  htrg->Draw("same");
  hdet->Draw("same");
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  TLegend * leg = new TLegend(0.6,0.75,0.9,0.9);
  leg->AddEntry(hhit,"Flux hitting lead","l");
  leg->AddEntry(htrg,"Triggered","l");
  leg->AddEntry(hdet,"Detected","l");
  leg->Draw();
}
