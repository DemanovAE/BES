//
// The StFemtoXi class holds indices of the StFemtoV0 and StFemtoTrack arrays
// which correspond to Xi daughters: Xi^- -> Lambda + pi^-
//

// C++ headers
#include <iostream>

// FemtoDst headers
#include "StFemtoXi.h"

ClassImp(StFemtoXi)

//________________
StFemtoXi::StFemtoXi() : TObject(), mV0Id(0), mTrackId(0) {
  /* empty */
}

//________________
StFemtoXi::StFemtoXi(const StFemtoXi &xi) :
  TObject(), mV0Id( xi.mV0Id ), mTrackId( xi.mTrackId ) {
  /* empty */
}

//________________
void StFemtoXi::print() {
  std::cout << "Xi daughters. V0: " << mV0Id << " Track: " << mTrackId
            << std::endl;
}
