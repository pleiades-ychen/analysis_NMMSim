// add cuts, and/or change quantity & title 

{
#include <algorithm>

  if (!TClass::GetDict("NMMEvent")) {
    gROOT->ProcessLine(".L NMMEvent.cc+");
    gROOT->ProcessLine(".L ph_likelihood.cc+");
  }

  TH1F * h1 = new TH1F("h1","Multiplicity;Multiplicity;Events",200,0,200);
  TH1F * h2 = new TH1F("h2","Multiplicity;Multiplicity;Events",200,0,200);
  TH1F * h3 = new TH1F("h3","Multiplicity;Multiplicity;Events",200,0,200);
  TH1F * h4 = new TH1F("h4","Multiplicity;Multiplicity;Events",200,0,200);
  TH1F * h5 = new TH1F("h5","Multiplicity;Multiplicity;Events",200,0,200);

  TH1F * hSum = new TH1F("hSum","Multiplicity as MC truth;Multiplicity;Events",200,0,200);

  TH1F * h1L = new TH1F("h1L",";Neutron Likelihood Fraction; Events",
			120,-0.1,1.1);
  TH1F * h2L = new TH1F("h2L",";Neutron Likelihood Fraction; Events",
			120,-0.1,1.1);
  TH1F * h3L = new TH1F("h3L",";Neutron Likelihood Fraction; Events",
			120,-0.1,1.1);
  
  //TH1F * hPH = new TH1F("hPH",";Pulse Height (mV);Counts",500,0,500);
  //TH1F * hPHc = new TH1F("hPHc",";Pulse Height (mV);Counts",500,0,500);

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
    
    bool cPH300 = true;
    for (int i = 0; i  < ph.size(); i++) 
      if (ph[i] > 300) {cPH300 = false; break;}
    
    string type = nmm->type;
    int M = nmm->M;

    //if (type=="contam_other"){
    //for (int i = 0; i  < ph.size(); i++) 
    //hPH->Fill(ph[i]);
    //}
    
    if (N>=4 && cPH300){
      
      double LN = ph_likelihood(ph,"n");
      double LG = ph_likelihood(ph,"g");
      double L = LN/(LN + LG);

      if (type=="pfn") h1L->Fill(L);
      if (type=="contam_muon") h2L->Fill(L);
      if (type=="contam_other") h3L->Fill(L);

      if (L>0.9){
	hSum->Fill(M);
	if (type=="pfn") h1->Fill(M); 
	if (type=="contam_muon") h2->Fill(M);
	if (type=="contam_other") h3->Fill(M); 
      
      }

      //if (type=="contam_other"){
      //for (int i = 0; i  < ph.size(); i++) 
      //  hPHc->Fill(ph[i]);
      //}

    }
  }

  new TCanvas;
  hSum->SetLineColor(kBlack);
  hSum->SetLineWidth(2);
  hSum->SetFillStyle(3005);
  hSum->SetFillColor(kBlack);
  h1->SetLineColor(kBlue);
  h1->SetTitle("Multiplicity;Multiplicity;Events");
  h2->SetLineColor(kMagenta);
  h3->SetLineColor(kRed);
  h1->SetLineWidth(2);
  h2->SetLineWidth(2);
  h3->SetLineWidth(2);
  h4->SetLineColor(kGreen);
  h5->SetLineColor(kCyan);
  hSum->Draw();
  h1->Draw("same");
  h2->Draw("same");
  h3->Draw("same");
  //h4->Draw("same");
  //h5->Draw("same");
  gPad->SetLogy();
  
  TLegend leg(0.6,0.65,0.88,0.85);
  leg.AddEntry(h1,"Pure Fast Neutron","l");
  leg.AddEntry(h2,"Muon Involved Contaminated","l");
  leg.AddEntry(h3,"Other Involved Contaminated","l");
  leg.AddEntry(hSum,"All Events in ROI","f");
  //leg.AddEntry(h4,"Muon Involved Background","l");
  //leg.AddEntry(h5,"Other Involved Background","l");
  leg.Draw();

  cout << "Pure Fast Neutron Events:\t" << h1->Integral(1,201) << "\n";
  cout << "Muon Involved Contaminated Events:\t" << h2->Integral(1,201) << "\n";
  cout << "Other Involved Contaminated Events:\t" << h3->Integral(1,201) << "\n";
  cout << "Muon Involved Background Events:\t" << h4->Integral(1,201) << "\n";
  cout << "Other Involved Background Events:\t" << h5->Integral(1,201) << "\n";

  new TCanvas;
  h1L->SetLineColor(kBlue);
  h1L->SetLineWidth(2);
  h2L->SetLineColor(kMagenta);
  h2L->SetLineWidth(2);
  h3L->SetLineColor(kRed);
  h3L->SetLineWidth(2);
  h1L->Draw();
  h3L->Draw("same");
  h2L->Draw("same");
  TLegend leg2(0.6,0.65,0.88,0.85);
  leg2.AddEntry(h1L,"Fast Neutron Events","l");
  leg2.AddEntry(h2L,"Muon Involved Contaminated","l");
  leg2.AddEntry(h3L,"Other Involved Contaminated","l");
  leg2.Draw();
  gPad->SetLogy();
}
