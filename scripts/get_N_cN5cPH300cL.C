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

  TH1F * hSum = new TH1F("hSum","Observable Multiplicity;Multiplicity;Events",200,0,200);

  TH1F * hSumM = new TH1F("hSumM","All Events in ROI;Multiplicity;Events",200,0,200);
  // for comparison

  TH1F * h1Sk = new TH1F("h1Sk",";(Observable - True) Multiplicity;Events",
			 32,-1,30);
  TH1F * h2Sk = new TH1F("h2Sk",";(Observable - True) Multiplicity;Events",
			 32,-1,30);
  TH1F * h3Sk = new TH1F("h3Sk",";(Observable - True) Multiplicity;Events",
			 32,-1,30);
  TH1F * hSumSk = new TH1F("hSumSk",
			   "Skewness;(Observable - True) Multiplicity;Events",
			   32,-1,30);

  TH2F *h1SkM = new TH2F("h1SkM","Skewness vs. Observable Multiplicity"
			 ";Observable Multiplicity;Skewness",
			 31,0,30,7,-1,5);
  TH2F *h2SkM = new TH2F("h2SkM","Skewness vs. Observable Multiplicity"
			 ";Observable Multiplicity;Skewness",
			 31,0,30,7,-1,5);
  TH2F *h3SkM = new TH2F("h3SkM","Skewness vs. Observable Multiplicity"
			 ";Observable Multiplicity;Skewness",
			 31,0,30,7,-1,5);
  TH2F *hSumSkM = new TH2F("hSumSkM","Skewness vs. Observable Multiplicity"
			   ";Observable Multiplicity;Skewness",
			   31,0,30,7,-1,5);
  
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
    
    if (N>=5 && cPH300){
      
      double LN = ph_likelihood(ph,"n");
      double LG = ph_likelihood(ph,"g");
      double L = LN/(LN + LG);

      if (L>0.99){
	hSumM->Fill(M);
	hSum->Fill(N);
	if (type=="pfn") h1->Fill(N); 
	if (type=="contam_muon") h2->Fill(N);
	if (type=="contam_other") h3->Fill(N); 

	hSumSk->Fill(N-M);
	if (type=="pfn") h1Sk->Fill(N-M); 
	if (type=="contam_muon") h2Sk->Fill(N-M);
	if (type=="contam_other") h3Sk->Fill(N-M); 

	hSumSkM->Fill(N,N-M);
	if (type=="pfn") h1SkM->Fill(N,N-M); 
	if (type=="contam_muon") h2SkM->Fill(N,N-M);
	if (type=="contam_other") h3SkM->Fill(N,N-M); 

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
  hSumSk->SetLineColor(kBlack);
  hSumSk->SetLineWidth(2);
  hSumSk->SetFillStyle(3005);
  hSumSk->SetFillColor(kBlack);
  h1Sk->SetLineColor(kBlue);
  h1Sk->SetTitle("Multiplicity;Multiplicity;Events");
  h2Sk->SetLineColor(kMagenta);
  h3Sk->SetLineColor(kRed);
  h1Sk->SetLineWidth(2);
  h2Sk->SetLineWidth(2);
  h3Sk->SetLineWidth(2);
  hSumSk->Draw();
  h1Sk->Draw("same");
  h2Sk->Draw("same");
  h3Sk->Draw("same");
  gPad->SetLogy();
  
  TLegend leg2(0.6,0.65,0.88,0.85);
  leg2.AddEntry(h1Sk,"Pure Fast Neutron","l");
  leg2.AddEntry(h2Sk,"Muon Involved Contaminated","l");
  leg2.AddEntry(h3Sk,"Other Involved Contaminated","l");
  leg2.AddEntry(hSumSk,"All Events in ROI","f");
  leg2.Draw();

  //plot 2-D 
  new TCanvas;
  gStyle->SetPalette(1,0);
  h1SkM->Draw("colz");
  new TCanvas;
  h2SkM->Draw("colz");
  new TCanvas;
  h3SkM->Draw("colz");
  new TCanvas;
  hSumSkM->Draw("colz");

  //comparison between observable and true multiplicity
  new TCanvas;
  //hSum->SetLineColor(kRed);
  //hSum->SetFillStyle(0);
  hSumM->SetLineColor(kBlue);
  hSumM->SetLineWidth(2);
  hSum->SetTitle("All Events in ROI;Multiplicity;Events");
  hSum->Draw();
  hSumM->Draw("same");
  TLegend leg3(0.6,0.65,0.88,0.85);
  leg3.AddEntry(hSum,"Observable Multiplicity","l");
  leg3.AddEntry(hSumM,"True Multiplicity","l");
  leg3.Draw();
}
