// add cuts, and/or change quantity & title 
// 10/08/2015
// use new muon cut: cut on the first two pulses at ph < 200 mV
// plot likelihood distribution
// count numbers of fast neutron events, 
// for varying Multiplicity
// study M = 4 events (in case of M=5 data):
//   1) 4n + 1g, in which g is indicated by this MC
//   2) 4n in MC, count the number to be used in 4n+1g estimate in real data,
//      in which the g is an accidental gamma from U/Th

{
#include <algorithm>

  if (!TClass::GetDict("NMMEvent")) {
    gROOT->ProcessLine(".L ../second_proc/NMMEvent.cc+");
  }
  gROOT->ProcessLine(".L ph_likelihood.cc+");

  TH1F * h1 = new TH1F("h1","Pulse Height Likelihood; n/(n+g);Events",120,-0.1,1.1);
  TH1F * h2 = new TH1F("h2","Pulse Height Likelihood; n/(n+g);Events",120,-0.1,1.1);
  TH1F * h3 = new TH1F("h3","Pulse Height Likelihood; n/(n+g);Events",120,-0.1,1.1);
  TH1F * h1b = new TH1F("h1b","Pulse Height Likelihood; n/(n+g);Events",120,-0.1,1.1);
  TH1F * h2b = new TH1F("h2b","Pulse Height Likelihood; n/(n+g);Events",120,-0.1,1.1);
  TH1F * h3b = new TH1F("h3b","Pulse Height Likelihood; n/(n+g);Events",120,-0.1,1.1);

  TH1F * hBkgMu = new TH1F("hBkgMu","Pulse due to background;Pulse Height (mV);Counts / 20 mV",15,0,300);
  TH1F * hBkgG = new TH1F("hBkgG","Pulse due to background;Pulse Height (mV);Counts / 20 mV",15,0,300);


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
    int M = nmm->M;
    if(M!=4) continue; // change this to work on different multiplicities
    vector<double> ph = nmm->pulseHeights;
    //double phmax = 0.;
    //for (int i = 0; i < N; i++) {
    //if(ph[i] > phmax)
    //phmax = ph[i];
    //}
    double ph12max = 0.;
    ph12max = ph[0] > ph[1] ? ph[0] : ph[1];
    
    string type = nmm->type;
    vector<string> tag = nmm->pulseTagV2;
   
    //if(phmax < 1){
    //for(int = 0; i < N; i++)
    //cout<<ph[i]<<", ";
    //cout<< "\n";
    //}

    //bool cPH300 = true;
    //for (int i = 0; i  < ph.size(); i++) 
    //if (ph[i] > 300) {cPH300 = false; break;} 
    
    if (ph12max < 200){
      //now modify to plot L and pulse heights
      //
      double nllh ,gllh;
      nllh = ph_likelihood(ph,"n","l");
      gllh = ph_likelihood(ph,"g","l");
      double L = nllh / (nllh + gllh);
      
      if (N==M){
	if (type=="pfn") h1->Fill(L);
	if (type=="contam_muon") h2->Fill(L);
	if (type=="contam_other") h3->Fill(L);
      }
      else if(N==(M+1)){
	if (type=="pfn") h1b->Fill(L);
	if (type=="contam_muon") h2b->Fill(L);
	if (type=="contam_other") h3b->Fill(L);
      }

      for (int i = 0; i < ph.size(); i++){
	if (tag[i] == "MuonFromCav"){
	  hBkgMu->Fill(ph[i]);
	  cout << "Position of the muon bkg pulse: " << i+1 << "\n";
	}

	if (tag[i] == "OtherFromCav" ||
	    tag[i] == "promptX"){
	  hBkgG->Fill(ph[i]);
	  cout << "Position of the gamma-like bkg pulse: " << i+1 << "\n";
	}
      }

      //if (type=="contam_other"){
      //for (int i = 0; i  < ph.size(); i++) 
      //  hPHc->Fill(ph[i]);
      //}

    }
  }

  double ScFactor = 158.23/217.36;

  new TCanvas;
  h1->SetLineColor(kBlue);
  h1->SetTitle("Pulse Height Likelihood; n/(n+g);Events");
  h2->SetLineColor(kMagenta);
  h3->SetLineColor(kRed);
  h1->SetLineWidth(2);
  h2->SetLineWidth(2);
  h3->SetLineWidth(2);
  h1->Draw();
  h2->Draw("same");
  h3->Draw("same");
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();  
  TLegend leg(0.6,0.65,0.88,0.85);
  leg.AddEntry(h1,"Fast Neutron","l");
  leg.AddEntry(h2,"Muon Involved Contaminated","l");
  leg.AddEntry(h3,"Other Involved Contaminated","l");
  leg.Draw();

  cout << "Events with N = M:\n";
  cout << "Pure Fast Neutron Events:\t" << h1->Integral() << "+/-" 
       << sqrt(h1->Integral()) <<"\n";
  cout << "Normalized to 158 live days:\t" << ScFactor*(h1->Integral()) << "+/-" 
       << ScFactor*sqrt(h1->Integral()) <<"\n";
  //cout << "Muon Involved Contaminated Events:\t" << h2->Integral() << "\n";
  //cout << "Other Involved Contaminated Events:\t" << h3->Integral() << "\n";

  new TCanvas;
  h1b->SetLineColor(kBlue);
  h1b->SetTitle("Pulse Height Likelihood; n/(n+g);Events");
  h2b->SetLineColor(kMagenta);
  h3b->SetLineColor(kRed);
  h1b->SetLineWidth(2);
  h2b->SetLineWidth(2);
  h3b->SetLineWidth(2);
  h1b->Draw();
  h2b->Draw("same");
  h3b->Draw("same");
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();  
  TLegend leg_1b(0.6,0.65,0.88,0.85);
  leg_1b.AddEntry(h1b,"Fast Neutron","l");
  leg_1b.AddEntry(h2b,"Muon Involved Contaminated","l");
  leg_1b.AddEntry(h3b,"Other Involved Contaminated","l");
  leg_1b.Draw();

  cout << "Events with N = M + 1:\n";
  //cout << "Pure Fast Neutron Events:\t" << h1b->Integral() << "\n";
  cout << "Muon Involved Contaminated Events:\t" << h2b->Integral() << "\n";
  cout << "Other Involved Contaminated Events:\t" << h3b->Integral() << "\n";
  cout << "Total 4n+1contam events:\t" << (h2b->Integral() + h3b->Integral()) << "+/-"
       << sqrt(h2b->Integral() + h3b->Integral()) <<"\n";
  cout << "Normalized to 158 live days:\t" << ScFactor*(h2b->Integral() + h3b->Integral()) << "+/-"
       << ScFactor*sqrt(h2b->Integral() + h3b->Integral()) <<"\n";

  new TCanvas;
  hBkgMu->SetLineColor(kMagenta);
  hBkgG->SetLineColor(kRed);
  hBkgG->Draw();
  hBkgMu->Draw("same");
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

  /*
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
  */
}
