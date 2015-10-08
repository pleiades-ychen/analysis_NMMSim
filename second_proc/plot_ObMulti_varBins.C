// add cuts, and/or change quantity & title 

{
#include <algorithm>

  if (!TClass::GetDict("NMMEvent")) {
    gROOT->ProcessLine(".L NMMEvent.cc+");
    gROOT->ProcessLine(".L ph_likelihood.cc+");
  }

  // multiplicity from real data
  double xc[10] = {5, 6, 7, 8, 9, 10, 13, 18, 23, 28};
  //double nc[10] = {42.0, 34.0, 30.0, 14.0, 11.0, 8.0, 3.6, 1.2, 1.2, 0.4};
  double nc[10] = {42.0, 34.0, 30.0, 14.0, 11.0, 8.0, 18.0, 6.0, 6.0, 2.0};
  //yet to normalize the wide bins (6,7,8,9)
  double ex[10] = {0.,0.,0.,0.,0.,0.,2.5,2.5,2.5,2.5};
  double ey[10];
  for (int i = 0; i < 6; i++){
    xc[i] += 0.5;
    ey[i] = sqrt(nc[i]);
  }
  for (int i = 6; i < 10; i++){
    xc[i] += 0.5;
    ey[i] = sqrt(nc[i])/5.;
    nc[i] /= 5.; 
  }
  
  TGraphErrors gr(10,xc,nc,ex,ey);

  double xl [16] = {0,1,2,3,4,5,6,7,8,9,10,11,16,21,26,31};
  TH1F * hSum = new TH1F("hSum","Observable Multiplicity;Multiplicity;Events",
			 15,xl);

  TH1F * hSumM = new TH1F("hSumM","All Events in ROI;Multiplicity;Events",
			  15,xl);
  // for comparison

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
    
    if (N>=5 && cPH300){
      
      double LN = ph_likelihood(ph,"n");
      double LG = ph_likelihood(ph,"g");
      double L = LN/(LN + LG);

      if (L>0.99){
	hSumM->Fill(M);
	hSum->Fill(N);
      }

    }
  }

  // normalize the wide bins 
  for (int i = 12; i < 17; i++){
    double newval = hSum->GetBinContent(i)/5.;
    double newvalM = hSumM->GetBinContent(i)/5.;
    hSum->SetBinContent(i,newval);
    hSumM->SetBinContent(i,newvalM);
  }
  //comparison between observable and true multiplicity
  new TCanvas;
  hSum->Scale(158.23/217.36);
  hSumM->Scale(158.23/217.36);
  hSum->SetLineColor(kBlack);
  hSum->SetFillStyle(0);
  hSum->SetLineWidth(2);
  hSumM->SetLineColor(kBlue);
  hSumM->SetLineWidth(2);
  hSum->SetTitle("Multiplicity;Multiplicity;Events / Multiplicity");
  //hSumM->Draw();
  hSum->Draw();
  gr.SetMarkerColor(kRed);
  gr.SetLineColor(kRed);
  gr.SetLineWidth(2);
  gr.SetMarkerStyle(4);
  gr.Draw("Psame");
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  TLegend leg3(0.6,0.65,0.88,0.85);
  leg3.AddEntry(hSum,"Simulation","l");
  //leg3.AddEntry(hSumM,"True Multiplicity","l");
  leg3.AddEntry(&gr,"Data","pl");
  leg3.Draw();
}
