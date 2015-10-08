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

  //TH1F * hPH = new TH1F("hPH","Other Involved Contaminated Events;Pulse Height (mV);Counts",500,0,500);
  TH1F * hPH1 = new TH1F("hPH1","Fast Neutron Events;Pulse Height (mV);Counts",50,0,1000);
  TH1F * hPH2 = new TH1F("hPH2","Muon Involved Contaminated Events;Pulse Height (mV);Counts",50,0,1000);
  TH1F * hPH3 = new TH1F("hPH3","Other Involved Contaminated Events;Pulse Height (mV);Counts",50,0,1000);
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
    if(N<5) continue;
    vector<double> ph = nmm->pulseHeights;
    //find the max
    //double phmax = std::max_element(ph.begin(), ph.end());
    double phmax = 0.;
    for (int i = 0; i < N; i++) {
      if(ph[i] > phmax)
	phmax = ph[i];
    }
    
    string type = nmm->type;
    if (type=="pfn")
      hPH1->Fill(phmax);
    else if (type=="contam_muon")
      hPH2->Fill(phmax);
    else if (type=="contam_other")
      hPH3->Fill(phmax);

    //if(phmax < 1){
    //for(int = 0; i < N; i++)
    //cout<<ph[i]<<", ";
    //cout<< "\n";
    //}

    bool cPH300 = true;
    for (int i = 0; i  < ph.size(); i++) 
    if (ph[i] > 300) {cPH300 = false; break;}
    
    //string type = nmm->type;
    int M = nmm->M;

    //if (type=="contam_other"){
    //for (int i = 0; i  < ph.size(); i++) 
    //hPH->Fill(ph[i]);
    //}
    
    if (N>=5 && cPH300){
      if (type=="pfn") h1->Fill(N);
      if (type=="contam_muon") h2->Fill(N);
      if (type=="contam_other") h3->Fill(N);
      if (type=="bkg_muon") h4->Fill(N);
      if (type=="bkg_other") h5->Fill(N);

      //if (type=="contam_other"){
      //for (int i = 0; i  < ph.size(); i++) 
      //  hPHc->Fill(ph[i]);
      //}

    }
  }

  new TCanvas;
  h1->SetLineColor(kBlue);
  h1->SetTitle("Multiplicity;Multiplicity;Events");
  h2->SetLineColor(kMagenta);
  h3->SetLineColor(kRed);
  h4->SetLineColor(kGreen);
  h5->SetLineColor(kCyan);
  h1->SetLineWidth(2);
  h2->SetLineWidth(2);
  h3->SetLineWidth(2);
  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");
  //h4->Draw("same");
  //h5->Draw("same");
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  
  TLegend leg(0.6,0.65,0.88,0.85);
  leg.AddEntry(h1,"Fast Neutron","l");
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

  /*
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
  gPad->SetGridx();
  gPad->SetGridy();
  */

  new TCanvas;
  hPH1->SetLineColor(kBlue);
  hPH1->SetTitle("Dist. of Max Pulse Heights;Pulse Height (mV);Events");
  hPH2->SetLineColor(kMagenta);
  hPH3->SetLineColor(kRed);
  hPH1->SetLineWidth(2);
  hPH2->SetLineWidth(2);
  hPH3->SetLineWidth(2);
  hPH1->Draw();
  hPH2->Draw("same");
  hPH3->Draw("same");
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  
  TLegend leg3(0.6,0.65,0.88,0.85);
  leg3.AddEntry(hPH1,"Fast Neutron","l");
  leg3.AddEntry(hPH2,"Muon Involved Contaminated","l");
  leg3.AddEntry(hPH3,"Other Involved Contaminated","l");
  leg3.Draw();

}
