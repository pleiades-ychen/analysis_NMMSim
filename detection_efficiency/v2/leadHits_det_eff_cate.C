{
  if (!TClass::GetDict("NMMEvent")) {
    gROOT->ProcessLine(".L NMMEvent.cc+");
  }

  double thr = 0.;
  double MuonCut = 300.;

  int totLeadHits = 0;
  int numTrg = 0;
  int numDet = 0;
  
  TH1F *hhit = new TH1F("hhit","Signal Events (also tagged by lower panels)"
			";Primary Neutron Energy (MeV);Events / MeV",
			1000,0,1000);
  TH1F *htrg = new TH1F("htrg",";Primary Energy (MeV);Events / MeV",
			1000,0,1000);
  TH1F *hdet = new TH1F("hdet",";Primary Energy (MeV);Events / MeV",
			1000,0,1000);

  //TH1F *hph12 = new TH1F("hph12","Signal Events (also tagged by lower panels)"
  //		";Max of First Two Pulses; Counts / mV",
  //		8000,0,8000);

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
    vector<double> nf = nmmEvent->northFraction;
    string type = nmmEvent->type;
    int M = nmmEvent->M;
    int TankMuon = nmmEvent-> TankMuon;
    int LeadMuon = nmmEvent->LeadMuon;
    double LeadMuonEnergy = nmmEvent->LeadMuonEnergy;
    int LeadNeutron = nmmEvent->LeadNeutron;
    double LeadNeutronEnergy = nmmEvent->LeadNeutronEnergy;
    int PanelMuon = nmmEvent-> PanelMuon;
    int ClipMuon = nmmEvent-> ClipMuon;
    
    int Ns = 0; int Nn = 0;
    for(int p = 0; p < nf.size(); p++){
      if(nf[p]<0.5)
	Ns++;
      else
	Nn++;
    }
    if(Ns + Nn != N)
      cout<< "Error: Ns + Nn doesn't match to total.\n";
    //if(PanelMuon==1){// all lower panel firing
    if(TankMuon==0 && LeadMuon==0 && LeadNeutron==1 && LeadNeutronEnergy > thr){ // Joel's Signal
    //if(PanelMuon==1 && TankMuon==0 && LeadMuon==0 && LeadNeutron==1 && LeadNeutronEnergy > thr){ // Ray's Signal
    //if((TankMuon==0 || ClipMuon==1) && LeadMuon==1){ // Joel's Background
    //if(PanelMuon==1 && (TankMuon==0 || ClipMuon==1) && LeadMuon==1){ // Ray's Background
      totLeadHits++;
      //cout<<"LeadNeutronEnergy = "<<LeadNeutronEnergy<<"\n";
      hhit->Fill(LeadNeutronEnergy);   // LeadNeutronEnergy 
      //if (N >= 5) { // normal trigger condition
      if(Ns>=2 && Nn>=2){ // Chiranjibi's trigger condition
	numTrg++;
	//cout<<"Triggered LeadNeutronEnergy = "<<LeadNeutronEnergy<<"\n";
	htrg->Fill(LeadNeutronEnergy); // LeadNeutronEnergy 

	//hph12->Fill(ph[0]>ph[1] ? ph[0]:ph[1]);
	cout << "1st Two Pulse Heights of Triggered: "
	     <<ph[0]<<", "<<ph[1]<<"\n";
 
	if(ph[0] < MuonCut && ph[1] < MuonCut){ // detected
	  cout<<"Detected LeadNeutronEnergy = "<<LeadNeutronEnergy<<"\n";
	  numDet++;
	  hdet->Fill(LeadNeutronEnergy); // LeadNeutronEnergy
	}
      }
    }
      
  }

  // efficiency for triggering
  double eff = (double)numTrg/totLeadHits;
  double err = sqrt(eff*(1-eff)/totLeadHits); // error estimate base on binomial stats
  cout<< "Eff = " << eff <<" +/- "<< err <<"\n";
  // efficiency for triggering && pass muon cut detection
  double eff2 = (double)numDet/totLeadHits;
  double err2 = sqrt(eff2*(1-eff2)/totLeadHits); 
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

  //new TCanvas;
  //hph12->SetLineColor(2);
  //hph12->Draw();
}
