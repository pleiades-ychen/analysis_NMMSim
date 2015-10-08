#include <iostream.h>
#include "TObject.h"

class NMMEvent : public TObject
{
public:
  int N; // Pulse Number = observable M
  vector<double> pulseHeights;
  vector<double> pulseTimes;
  vector<string> pulseTagV2;
  string type;
  int M; // MC True M
  NMMEvent();
  ~NMMEvent();
  ClassDef(NMMEvent, 1); 
};
