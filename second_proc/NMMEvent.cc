#include "NMMEvent.hh"

#if !defined(__CINT__)
ClassImp(NMMEvent);
#endif

NMMEvent::NMMEvent()
{
  N = 0;
  type = "emptyEvent";
  M = 0;
}
