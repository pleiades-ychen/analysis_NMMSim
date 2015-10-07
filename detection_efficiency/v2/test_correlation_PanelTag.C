{
  TH2F * corr = new TH2F("corr","Correlation between"
			 " LeadMuon Tag and PanelMuon Tag;PanelMuon;LeadMuon"
			 ,3,0,3,3,0,3);
  
  TFile * fin = new TFile("melindas_neutrons_PulseTags_DetEff.root","READ");
  TTree * events = (TTree*)fin->Get("Events");

  int lm;
  int pm;
  events->SetBranchAddress("LeadMuon",&lm);
  events->SetBranchAddress("PanelMuon",&pm);
  
  int num = events->GetEntries();
  for(int i = 0; i < num; i++){
    events->GetEntry(i);
    corr->Fill(pm,lm);

  }

  new TCanvas;
  corr->Draw("colz");
  gPad->SetLogz();
}
