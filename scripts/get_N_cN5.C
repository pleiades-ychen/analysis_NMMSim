// add cuts, and/or change quantity & title 

{
  if (!TClass::GetDict("NMMEvent")) {
    gROOT->ProcessLine(".L NMMEvent.cc+");
  }
  TFile *fin = new TFile("melindas_neutrons_2ndProc.root","READ");
  
  TTree* event = (TTree*) fin->Get("Event");
  new TCanvas;
  
  event->Draw("NMMEvent.N>>h1(200,0,200)",
	      "NMMEvent.type==\"pfn\""
	      "&& NMMEvent.N>=5");
  event->Draw("NMMEvent.N>>h2(200,0,200)",
	      "NMMEvent.type==\"contam_muon\""
	      "&& NMMEvent.N>=5");
  event->Draw("NMMEvent.N>>h3(200,0,200)",
	      "NMMEvent.type==\"contam_other\""
	      "&& NMMEvent.N>=5");
  event->Draw("NMMEvent.N>>h4(200,0,200)",
	      "NMMEvent.type==\"bkg_muon\""
	      "&& NMMEvent.N>=5");
  event->Draw("NMMEvent.N>>h5(200,0,200)",
	      "NMMEvent.type==\"bkg_other\""
	      "&& NMMEvent.N>=5");
  TH1F *h1 = (TH1F*) gDirectory->Get("h1");
  TH1F *h2 = (TH1F*) gDirectory->Get("h2");
  TH1F *h3 = (TH1F*) gDirectory->Get("h3");
  TH1F *h4 = (TH1F*) gDirectory->Get("h4");
  TH1F *h5 = (TH1F*) gDirectory->Get("h5");
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
  gPad->SetGridx();
  gPad->SetGridy();
  
  TLegend leg(0.6,0.65,0.88,0.85);
  leg.AddEntry(h1,"Fast Neutron Events","l");
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

}
