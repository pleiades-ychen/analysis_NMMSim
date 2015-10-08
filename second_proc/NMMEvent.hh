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
  NMMEvent();
  ~NMMEvent() { }
  ClassDef(NMMEvent, 1); 
};
