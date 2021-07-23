//
// The StFemtoV0 class holds indices of positive and
// negative daughters in StFemtoTrack array
//

// C++ headers
#include <iostream>

// FemtoDst headers
#include "StFemtoV0.h"

ClassImp(StFemtoV0)

//________________
StFemtoV0::StFemtoV0() : TObject(), mPosId(0), mNegId(0) {
  /* empty */
}

//________________
StFemtoV0::StFemtoV0(const StFemtoV0 &v0) :
  TObject(), mPosId(v0.mPosId), mNegId(v0.mNegId) {
  /* empty */
}

//________________
void StFemtoV0::print() {
  std::cout << "V0 daughters. PosId: " << mPosId << " NegId: "
            << mNegId << std::endl;
}
