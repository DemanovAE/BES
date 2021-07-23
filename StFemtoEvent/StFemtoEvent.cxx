//
// The StFemtoEvent holds information about the current event
//

// C++ headers
#include <limits>
#include <iostream>

// FemtoDst headers
#include "StFemtoEvent.h"

// ROOT headers
#include "TMath.h"

ClassImp(StFemtoEvent)

//_________________
StFemtoEvent::StFemtoEvent() :
  mEventId(0), mRunId(0), mRefMult(0), mRefMultCorr(0), mRefMultCorrWeight(0),
  mRefMultPos(0), mRefMult2(0), mRefMult2Pos(0), mGRefMult(0), mCent16(-1),
  mNumberOfPrimaryTracks(0), mNumberOfGlobalTracks(0),
  mNumberOfBTofHit(0), mNumberOfTofMatched(0),mNumberOfBEMCMatched(0),
  mMagField(0), mVertexPositionX(0), mVertexPositionY(0), mVertexPositionZ(0),
  mVpdVz(0), mPrimaryVertexRanking(0), mNumberOfPrimaryVertices(0),
  mTransverseSphericity(-1), mTransverseSphericity2(-1),
  mZdcSumAdcEast(0), mZdcSumAdcWest(0),
  mZdcSmdEastHorizontal{}, mZdcSmdEastVertical{}, mZdcSmdWestHorizontal{}, mZdcSmdWestVertical{},
  mBbcAdcEast{}, mBbcAdcWest{} { // Default constructor

  if( !mTriggerIds.empty() ) {
    mTriggerIds.clear();
  }

}

//_________________
StFemtoEvent::StFemtoEvent(const StFemtoEvent& ev) {
  // Copy constructor
  mEventId = ev.mEventId;
  mRunId = ev.mRunId;
  mRefMult = ev.mRefMult;
  mRefMultCorr = ev.mRefMultCorr;
  mRefMultCorrWeight = ev.mRefMultCorrWeight;
  mRefMultPos = ev.mRefMultPos;
  mGRefMult = ev.mGRefMult;
  mCent16 = ev.mCent16;
  mNumberOfPrimaryTracks = ev.mNumberOfPrimaryTracks;
  mNumberOfGlobalTracks = ev.mNumberOfGlobalTracks;
  mNumberOfBTofHit = ev.mNumberOfBTofHit;
  mNumberOfTofMatched = ev.mNumberOfTofMatched;
  mNumberOfBEMCMatched = ev.mNumberOfBEMCMatched;
  mMagField = ev.mMagField;
  mVertexPositionX = ev.mVertexPositionX;
  mVertexPositionY = ev.mVertexPositionY;
  mVertexPositionZ = ev.mVertexPositionZ;
  mVpdVz = ev.mVpdVz;
  mPrimaryVertexRanking = ev.mPrimaryVertexRanking;
  mNumberOfPrimaryVertices = ev.mNumberOfPrimaryVertices;
  mTransverseSphericity = ev.mTransverseSphericity;
  mTransverseSphericity2 = ev.mTransverseSphericity2;

  mZdcSumAdcEast = ev.mZdcSumAdcEast;
  mZdcSumAdcWest = ev.mZdcSumAdcWest;
  for(int iModule=0; iModule<8; iModule++) {
    mZdcSmdEastHorizontal[iModule] = ev.mZdcSmdEastHorizontal[iModule];
    mZdcSmdEastVertical[iModule] = ev.mZdcSmdEastVertical[iModule];
    mZdcSmdWestHorizontal[iModule] = ev.mZdcSmdWestHorizontal[iModule];
    mZdcSmdWestVertical[iModule] = ev.mZdcSmdWestVertical[iModule];
  }

  for(int iModule=0; iModule<24; iModule++) {
    mBbcAdcEast[iModule] = ev.mBbcAdcEast[iModule];
    mBbcAdcWest[iModule] = ev.mBbcAdcWest[iModule];
  }

  for(unsigned int iIter=0; ev.mTriggerIds.size(); iIter++) {
    mTriggerIds.at(iIter) = ev.mTriggerIds.at(iIter);
  }
}

//_________________
Bool_t StFemtoEvent::isTrigger(const unsigned int& id) const {
  // Check if id is in the list of triggers
  return std::find(mTriggerIds.begin(), mTriggerIds.end(), id) != mTriggerIds.end();
}

//_________________
void StFemtoEvent::setPrimaryVertexPosition(const Float_t& x, const Float_t& y,
                                            const Float_t& z) {
  // Set primary vertex position
  mVertexPositionX = x;
  mVertexPositionY = y;
  mVertexPositionZ = z;
}

//_________________
void StFemtoEvent::setPrimaryVertexRanking(const Float_t& rank) {
  // Set primary vertex ranking according to the StFemtoEvent scheme
  if(rank < -1e9)       { mPrimaryVertexRanking = -10; }
  else if(rank < -1e6)  { mPrimaryVertexRanking = -9;  }
  else if(rank < -1e3)  { mPrimaryVertexRanking = -8;  }
  else if(rank < -1e2)  { mPrimaryVertexRanking = -7;  }
  else if(rank < -10)   { mPrimaryVertexRanking = -6;  }
  else if(rank < -5)    { mPrimaryVertexRanking = -5;  }
  else if(rank < -4)    { mPrimaryVertexRanking = -4;  }
  else if(rank < -3)    { mPrimaryVertexRanking = -3;  }
  else if(rank < -2)    { mPrimaryVertexRanking = -2;  }
  else if(rank < -1)    { mPrimaryVertexRanking = -1;  }
  else if(rank < 1)     { mPrimaryVertexRanking = 0;   }
  else if(rank < 2)     { mPrimaryVertexRanking = 1;   }
  else if(rank < 3)     { mPrimaryVertexRanking = 2;   }
  else if(rank < 4)     { mPrimaryVertexRanking = 3;   }
  else if(rank < 5)     { mPrimaryVertexRanking = 4;   }
  else if(rank < 10)    { mPrimaryVertexRanking = 5;   }
  else if(rank < 1e2)   { mPrimaryVertexRanking = 6;   }
  else if(rank < 1e3)   { mPrimaryVertexRanking = 7;   }
  else if(rank < 1e6)   { mPrimaryVertexRanking = 8;   }
  else if(rank < 1e9)   { mPrimaryVertexRanking = 9;   }
  else { mPrimaryVertexRanking = 10; }
}

//_________________
void StFemtoEvent::setRefMult(const Int_t& mult) {
  // Set refMult
  if( mult <= 0 ) {
    mRefMult = 0;
  }
  else {
    mRefMult = ( mult>std::numeric_limits<unsigned short>::max() ?
		 std::numeric_limits<unsigned short>::max() : (UShort_t)mult );
  }
}

//_________________
void StFemtoEvent::setRefMultPos(const Int_t& mult) {
  // Set refMult estimated by positive particles
  if( mult <= 0) {
    mRefMultPos = 0;
  }
  else {
    mRefMultPos = ( mult>std::numeric_limits<unsigned short>::max() ?
		    std::numeric_limits<unsigned short>::max() : (UShort_t)mult );
  }
}

//_________________
void StFemtoEvent::setRefMult2(const Int_t& mult) {
  // Set refMult2
  if( mult <= 0 ) {
    mRefMult2 = 0;
  }
  else {
    mRefMult2 = ( mult>std::numeric_limits<unsigned short>::max() ?
		  std::numeric_limits<unsigned short>::max() : (UShort_t)mult );
  }
}

//_________________
void StFemtoEvent::setRefMult2Pos(const Int_t& mult) {
  // Set refMult2 estimated by positive particles
  if( mult <= 0) {
    mRefMult2Pos = 0;
  }
  else {
    mRefMult2Pos = ( mult>std::numeric_limits<unsigned short>::max() ?
		     std::numeric_limits<unsigned short>::max() : (UShort_t)mult );
  }
}

//_________________
void StFemtoEvent::setRefMultCorr(const Float_t& mult) {
  // Set refMultCorr(ected) that comes from StRefMultCorr
  if( mult <= 0 ) {
    mRefMultCorr = 0;
  }
  else {
    mRefMultCorr = (  (mult * 10.)>std::numeric_limits<unsigned short>::max()  ?
		      std::numeric_limits<unsigned short>::max() :
		      (UShort_t)TMath::Nint( mult * 10. ) );
  }
}

//_________________
void StFemtoEvent::setNumberOfTofMatched(const Int_t& mult) {
  // Set number of primary tracks that matched TOF
  if( mult <= 0 ) {
    mNumberOfTofMatched = 0;
  }
  else {
    mNumberOfTofMatched = ( mult>std::numeric_limits<unsigned short>::max() ?
			    std::numeric_limits<unsigned short>::max() : (UShort_t)mult );
  }
}

//_________________
void StFemtoEvent::setNumberOfBEMCMatched(const Int_t& mult) {
  // Set number of primary tracks that matched BEMC
  if( mult <= 0 ) {
    mNumberOfBEMCMatched = 0;
  }
  else {
    mNumberOfBEMCMatched = ( mult>std::numeric_limits<unsigned short>::max() ?
			     std::numeric_limits<unsigned short>::max() : (UShort_t)mult );
  }
}

//_________________
void StFemtoEvent::setRefMultCorrWeight(const Float_t& w) {
  // Set refMultCorrWeight obtained from StRefMultCorr
  if( w<=0 ) {
    mRefMultCorrWeight = 0;
  }
  else {
    mRefMultCorrWeight = ( (w * 10000.) > std::numeric_limits<unsigned short>::max() ?
			   std::numeric_limits<unsigned short>::max() :
			   (UShort_t)TMath::Nint( w * 10000. ) );
  }
}

//_________________
void StFemtoEvent::setGRefMult(const Int_t& mult) {
  // Set gRefMult
  if( mult <= 0 ) {
    mGRefMult = 0;
  }
  else {
    mGRefMult = ( mult>std::numeric_limits<unsigned short>::max() ?
		  std::numeric_limits<unsigned short>::max() : (UShort_t)mult );
  }
}

//_________________
void StFemtoEvent::setNumberOfPrimaryTracks(const Int_t& mult) {
  // Set number of primary track reconstructed in the event
  if( mult <= 0 ) {
    mNumberOfPrimaryTracks = 0;
  }
  else {
    mNumberOfPrimaryTracks = ( mult>std::numeric_limits<unsigned short>::max() ?
			       std::numeric_limits<unsigned short>::max() :
			       (UShort_t)mult );
  }
}

//_________________
void StFemtoEvent::setNumberOfGlobalTracks(const Int_t& mult) {
  // Set number of global tracks reconstructed in the event
  if( mult <= 0 ) {
    mNumberOfGlobalTracks = 0;
  }
  else {
    mNumberOfGlobalTracks = ( mult>std::numeric_limits<unsigned short>::max() ?
			      std::numeric_limits<unsigned short>::max() :
			      (UShort_t)mult );
  }
}

//_________________
void StFemtoEvent::setNumberOfPrimaryVertices(const Int_t& mult) {
  // Set number of primary vertices reconstructed in the event
  if( mult <= 0 ) {
    mNumberOfPrimaryVertices = 0;
  }
  else {
    mNumberOfPrimaryVertices = ( mult>std::numeric_limits<unsigned short>::max() ?
				 std::numeric_limits<unsigned short>::max() :
				 (UShort_t)mult );
  }
}

//_________________
void StFemtoEvent::setNumberOfBTofHit(const Int_t& mult) {
  // Set number of hits in TOF (trays)
  if( mult <= 0 ) {
    mNumberOfBTofHit = 0;
  }
  else {
    mNumberOfBTofHit = ( mult>std::numeric_limits<unsigned short>::max() ?
			 std::numeric_limits<unsigned short>::max() :
			 (UShort_t)mult );
  }
}

//_________________
void StFemtoEvent::setTransverseSphericity(const Float_t& sph) {
  // Set transverse sphericity estimated in |eta|<0.5
  if(sph<0. || sph>1.) {
    mTransverseSphericity = -1;
  }
  else {
    mTransverseSphericity = (Char_t)(sph * 100.);
  }
}

//_________________
void StFemtoEvent::setTransverseSphericity2(const Float_t& sph) {
  // Set transverse sphericity estimated in |eta|<1
  if(sph<0. || sph>1.) {
    mTransverseSphericity2 = -1;
  }
  else {
    mTransverseSphericity2 = (UChar_t)(sph * 100.);
  }
}
//_________________
Int_t StFemtoEvent::cent9() const {
  // Return centality esitmated for 9 bins
  if(mCent16 == 15) {
    return 8;  //0-5%
  }
  else if(mCent16 == 14) {
    return 7;  //5-10%
  }
  else if( (mCent16 == 13) || (mCent16 == 12) ) {
    return 6;  //10-20%
  }
  else if( (mCent16 == 11) || (mCent16 == 10) ) {
    return 5;  //20-30%
  }
  else if( (mCent16 == 9) || (mCent16 == 8) ) {
    return 4;  //30-40%
  }
  else if( (mCent16 == 7) || (mCent16 == 6) ) {
    return 3;  //40-50%
  }
  else if( (mCent16 == 5) || (mCent16 == 4) ) {
    return 2;  //50-60%
  }
  else if( (mCent16 == 3) || (mCent16 == 2) ) {
    return 1;  //60-70%
  }
  else if( (mCent16 == 1) || (mCent16 == 0) ) {
    return 0;  //70-80%
  }
  else {
    return -1;
  }
}

//_________________
void StFemtoEvent::setZdcSumAdcEast(const Float_t& x) {
  // Set sum of ADC for east ZDC
  if( x < 0) {
    mZdcSumAdcEast = 0;
  }
  else {
    mZdcSumAdcEast = ( x>std::numeric_limits<unsigned short>::max() ?
		       std::numeric_limits<unsigned short>::max() :
		       (UShort_t)TMath::Nint( x ) );
  }
}

//_________________
void StFemtoEvent::setZdcSumAdcWest(const Float_t& x) {
  // Set sum of ADC for west ZDC
  if( x < 0) {
    mZdcSumAdcWest = 0;
  }
  else {
    mZdcSumAdcWest = ( x>std::numeric_limits<unsigned short>::max() ?
		       std::numeric_limits<unsigned short>::max() :
		       (UShort_t)TMath::Nint( x ) );
  }
}

//_________________
void StFemtoEvent::setZdcSmdEastHorizontal(const Int_t& iMod, const Float_t& x) {
  // Set ADC for east ZDC i-th horizontal strip
  if( x < 0 ) {
    mZdcSmdEastHorizontal[iMod] = 0;
  }
  else {
    mZdcSmdEastHorizontal[iMod] = ( x>std::numeric_limits<unsigned short>::max() ?
				    std::numeric_limits<unsigned short>::max() :
				    (UShort_t)TMath::Nint( x ) );
  }
}

//_________________
void StFemtoEvent::setZdcSmdEastVertical(const Int_t& iMod, const Float_t& x) {
  // Set ADC for east ZDC i-th vertical strip
  if( x < 0 ) {
    mZdcSmdEastVertical[iMod] = 0;
  }
  else {
    mZdcSmdEastVertical[iMod] = ( x>std::numeric_limits<unsigned short>::max() ?
				  std::numeric_limits<unsigned short>::max() :
				  (UShort_t)TMath::Nint( x ) );
  }
}

//_________________
void StFemtoEvent::setZdcSmdWestHorizontal(const Int_t& iMod, const Float_t& x) {
  // Set ADC for west ZDC i-th horizontal strip
  if( x < 0 ) {
    mZdcSmdWestHorizontal[iMod] = 0;
  }
  else {
    mZdcSmdWestHorizontal[iMod] = ( x>std::numeric_limits<unsigned short>::max() ?
				    std::numeric_limits<unsigned short>::max() :
				    (UShort_t)TMath::Nint( x ) );
  }
}

//_________________
void StFemtoEvent::setZdcSmdWestVertical(const Int_t& iMod, const Float_t& x) {
  // Set ADC for west ZDC i-th vertical strip
  if( x < 0 ) {
    mZdcSmdWestVertical[iMod] = 0;
  }
  else {
    mZdcSmdWestVertical[iMod] = ( x>std::numeric_limits<unsigned short>::max() ?
				  std::numeric_limits<unsigned short>::max() :
				  (UShort_t)TMath::Nint( x ) );
  }
}

//_________________
void StFemtoEvent::setBbcAdcEast(const Int_t& iMod, const UShort_t& x) {
  // Set ADC for east BBC i-th PMT box
  if( x < 0 ) {
    mBbcAdcEast[iMod] = 0;
  }
  else {
    mBbcAdcEast[iMod] = ( x>std::numeric_limits<unsigned short>::max() ?
			  std::numeric_limits<unsigned short>::max() :
			  (UShort_t)x );
  }
}

//_________________
void StFemtoEvent::setBbcAdcWest(const Int_t& iMod, const UShort_t& x) {
  // Set ADC for west BBC i-th PMT box
  if( x < 0 ) {
    mBbcAdcWest[iMod] = 0;
  }
  else {
    mBbcAdcWest[iMod] = ( x>std::numeric_limits<unsigned short>::max() ?
			  std::numeric_limits<unsigned short>::max() :
			  (UShort_t)x );
  }
}

//_________________
void StFemtoEvent::setTriggerId(const UInt_t& id) {
  // Add trigger id to the list if it is not there yet

  // Check if the trigger id is already in the mTriggerIds vector
  if ( std::find(mTriggerIds.begin(),
                 mTriggerIds.end(),
                 id) != mTriggerIds.end() ) {
    mTriggerIds.push_back(id);
  }
}

//_________________
void StFemtoEvent::addTriggerId(const UInt_t& id) {
  // Add trigger id to the list if it is not there yet
  setTriggerId(id);
}

//_________________
void StFemtoEvent::print() {
  // Print basic event information
  std::cout << "eventId: " << eventId() << " runId: " << runId() << " year: " << year() << " day: " << day()
	    << std::endl
	    << "refMult: " << refMult() << " refMult2: " << refMult2() << " refMultCorr: " << refMultCorr()
	    << " gRefMult: " << gRefMult() << std::endl
	    << "cent9: " << cent9() << " cent16: " << cent16() << std::endl
	    << "# bTofHit: " << numberOfBTofHit() << " # primaryTracks: " << numberOfPrimaryTracks()
	    << "# globalTracks: " << numberOfGlobalTracks() << " # tofMatched: " << numberOfTofMatched()
	    << std::endl
	    << "bField: " << magneticField() << " transverseSph: " << transverseSphericity()
	    << " transverseSph2: " << transverseSphericity2() << std::endl
	    << "Vertex position x,y,z: " << primaryVertex().X() << " , " << primaryVertex().Y()
	    << " , " << primaryVertex().Z() << std::endl
	    << "vpdVz: " << vpdVz() << " ranking: " << ranking() << std::endl
	    << "trigger Ids:";
  for(UInt_t iTrg=0; triggerIds().size(); iTrg++) {
    std::cout << " " << triggerIds().at(iTrg);
  }
  std::cout << std::endl;
}
