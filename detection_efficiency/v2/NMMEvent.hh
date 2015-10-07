// special version for detection efficiency study

#include <fstream>
#include <stdio.h>
#include <math.h>
#include <string>
#include <vector>
#include <iostream>
#include "TObject.h"

using namespace std;

class NMMEvent : public TObject
{
public:
  int N; // Pulse Number = observable M
  vector<double> pulseHeights;
  vector<double> pulseTimes;
  vector<double> northFraction;
  vector<string> pulseTagV2;
  string type;
  int M; // MC True M
  int TankMuon;
  int LeadMuon;
  double LeadMuonEnergy;
  int LeadNeutron;
  double LeadNeutronEnergy;
  int PanelMuon;
  int ClipMuon;
  NMMEvent();
  ~NMMEvent() { }
  ClassDef(NMMEvent, 1); 
};
