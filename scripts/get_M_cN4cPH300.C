// add cuts, and/or change quantity & title 

{
#include <algorithm>

  if (!TClass::GetDict("NMMEvent")) {
    gROOT->ProcessLine(".L NMMEvent.cc+");
  }

  TH1F * h1 = new TH1F("h1","Multiplicity;Multiplicity;Events",200,0,200);
  TH1F * h2 = new TH1F("h2","Multiplicity;Multiplicity;Events",200,0,200);
  TH1F * h3 = new TH1F("h3","Multiplicity;Multiplicity;Events",200,0,200);
  TH1F * h4 = new TH1F("h4","Multiplicity;Multiplicity;Events",200,0,200);
  TH1F * h5 = new TH1F("h5","Multiplicity;Multiplicity;Events",200,0,200);

  TH1F * hPH = new TH1F("hPH","Other Involved Contaminated Events;Pulse Height (mV);Counts",500,0,500);
  TH1F * hPHc = new TH1F("hPHc",";Pulse Height (mV);Counts",500,0,500);

  TFile *fin = new TFile("melindas_neutrons_2ndProc.root","READ");
  TTree* event = (TTree*) fin->Get("Event");
  NMMEvent *nmm;
  event->SetBranchAddress("NMMEvent",&nmm);

  int NumEntries = event->GetEntries();
  for (int n = 0; n < NumEntries; n++) {
     if (n%10000==0)
      cout<< "Processing event "<< n <<"\n";
    event->GetEntry(n);  
   
    int N = nmm->N;
    if(N==0) continue;
    vector<double> ph = nmm->pulseHeights;
    //double phmax = std::max_element(ph.begin(), ph.end());
    
    bool cPH300 = true;
    for (int i = 0; i  < ph.size(); i++) 
    if (ph[i] > 300) {cPH300 = false; break;}
    
    string type = nmm->type;
    int M = nmm->M;

    if (type=="contam_other"){
      for (int i = 0; i  < ph.size(); i++) 
	hPH->Fill(ph[i]);
    }
    
    if (N>=4 && cPH300){
      if (type=="pfn") h1->Fill(M);
      if (type=="contam_muon") h2->Fill(M);
      if (type=="contam_other") h3->Fill(M);
      if (type=="bkg_muon") h4->Fill(M);
      if (type=="bkg_other") h5->Fill(M);

      if (type=="contam_other"){
	for (int i = 0; i  < ph.size(); i++) 
	  hPHc->Fill(ph[i]);
      }

    }
  }

  new TCanvas;
  h1->SetLineColor(kBlue);
  h1->SetTitle("Multiplicity;Multiplicity;Events");
  h2->SetLineColor(kMagenta);
  h3->SetLineColor(kRed);
  h4->SetLineColor(kGreen);
  h5->SetLineColor(kCyan);
  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");
  //h4->Draw("same");
  //h5->Draw("same");
  gPad->SetLogy();
  
  TLegend leg(0.6,0.65,0.88,0.85);
  leg.AddEntry(h1,"Pure Fast Neutron","l");
  leg.AddEntry(h2,"Muon Involved Contaminated","l");
  leg.AddEntry(h3,"Other Involved Contaminated","l");
  //leg.AddEntry(h4,"Muon Involved Background","l");
  //leg.AddEntry(h5,"Other Involved Background","l");
  leg.Draw();

  cout << "Pure Fast Neutron Events:\t" << h1->Integral(1,201) << "\n";
  cout << "Muon Involved Contaminated Events:\t" << h2->Integral(1,201) << "\n";
  cout << "Other Involved Contaminated Events:\t" << h3->Integral(1,201) << "\n";
  cout << "Muon Involved Background Events:\t" << h4->Integral(1,201) << "\n";
  cout << "Other Involved Background Events:\t" << h5->Integral(1,201) << "\n";

  new TCanvas;
  hPH->SetLineColor(kBlack);
  hPHc->SetLineColor(kRed);
  hPH->Draw();
  hPHc->Draw("same");
  TLegend leg2(0.6,0.65,0.88,0.85);
  leg2.AddEntry(hPH,"Without muon cut","l");
  leg2.AddEntry(hPHc,"With muon cut","l");
  leg2.Draw();
  gPad->SetLogy();
}
