//
// The StFemtoTrack class stores track information
//

// C++ headers
#include <limits>
#include <cmath>
#include <iostream>

// FemtoDst headers
#include "StFemtoTrack.h"

// ROOT headers
#include "TMath.h"

ClassImp(StFemtoTrack)

//_________________
StFemtoTrack::StFemtoTrack() : TObject(),
  mId(0), mNHits(0), mNHitsPoss(0), mNHitsDedx(0),
  mDedx(0), mMap{}, mChi2(0),
  mNSigmaElectron(-60), mNSigmaPion(-60), mNSigmaKaon(-60), mNSigmaProton(-60),
  mBeta(0),
  mPrimMomX(0), mPrimMomY(0), mPrimMomZ(0),
  mGlobMomX(0), mGlobMomY(0), mGlobMomZ(0),
  mOriginX(0), mOriginY(0), mOriginZ(0) { // Default constructor
  /* empty */
}

//_________________
StFemtoTrack::StFemtoTrack(const StFemtoTrack &t) : TObject() {
  // Copy constructor
  mId = t.mId;
  mNHits = t.mNHits;
  mNHitsPoss = t.mNHitsPoss;
  mNHitsDedx = t.mNHitsDedx;
  mDedx = t.mDedx;
  mMap[0] = t.mMap[0];
  mMap[1] = t.mMap[1];
  mChi2 = t.mChi2;
  mNSigmaElectron = t.mNSigmaElectron;
  mNSigmaPion = t.mNSigmaPion;
  mNSigmaKaon = t.mNSigmaKaon;
  mNSigmaProton = t.mNSigmaProton;
  mBeta = t.mBeta;
  mPrimMomX = t.mPrimMomX;
  mPrimMomY = t.mPrimMomY;
  mPrimMomZ = t.mPrimMomZ;
  mGlobMomX = t.mGlobMomX;
  mGlobMomY = t.mGlobMomY;
  mGlobMomZ = t.mGlobMomZ;
  mOriginX = t.mOriginX;
  mOriginY = t.mOriginY;
  mOriginZ = t.mOriginZ;
}

//_________________
Float_t StFemtoTrack::massSqr() const {
  // Calculate squared mass of the track via TOF beta measurement
  Float_t mSqr = -999.;
  if( isTofTrack() && isPrimary() ) {
    mSqr = p() * p() * ( invBeta() * invBeta() - 1.);
  }

  return mSqr;
}

//_________________
void StFemtoTrack::setDedx(const Float_t& dEdx) {
  // Set dE/dx value
  if( (dEdx * 1e9) >= std::numeric_limits<unsigned short>::max() ) {
    mDedx = std::numeric_limits<unsigned short>::max();
  }
  else {
    mDedx = (UShort_t)(dEdx * 1e9);
  }
}

//_________________
void StFemtoTrack::setDedx(const Double_t& dEdx) {
  // Set dE/dx value
  if( (dEdx * 1e9) >= std::numeric_limits<unsigned short>::max() ) {
    mDedx =  std::numeric_limits<unsigned short>::max();
  }
  else {
    mDedx = (UShort_t)(dEdx * 1e9);
  }
}

//_________________
Float_t StFemtoTrack::gDCAz(const Float_t& pVtxZ) const {
  // Calculate DCA z of the global track to pVtx
  return (mOriginZ - pVtxZ);
}

//_________________
Float_t StFemtoTrack::gDCAxy(const Float_t& pVtxX, const Float_t& pVtxY) const {
  // Calculate DCA perp of the global track to pVtx
  return TMath::Sqrt( (mOriginX - pVtxX) * (mOriginX - pVtxX) +
		                  (mOriginY - pVtxY) * (mOriginY - pVtxY) );
}

//_________________
Float_t StFemtoTrack::gDCA(const Float_t& pVtxX, const Float_t& pVtxY,
                           const Float_t& pVtxZ) const {
  // Calculate DCA of the global track to pVtx
  return TMath::Sqrt( gDCAxy(pVtxX, pVtxY) * gDCAxy(pVtxX, pVtxY) +
		                  gDCAz(pVtxZ) * gDCAz(pVtxZ) );
}

//_________________
TVector3 StFemtoTrack::gDCA(TVector3 point) const {
  // Return DCA value between the point of DCA to isPrimary
  // vertex and the point
  return ( origin() - point );
}

//_________________
void StFemtoTrack::print() {
  std::cout << "Track info: " << std::endl
	    << "id: " << id() << " type: " << type() << " isPrimary: " << isPrimary()
	    << " flag: " << flag() << " charge: " << charge() << std::endl
	    << "nHits: " << nHits() << " nHitsPoss: " << nHitsPoss() << std::endl
	    << "ns(e): " << nSigmaElectron() << " ns(pi): " << nSigmaPion()
	    << " ns(K): " << nSigmaKaon() << " ns(p): " << nSigmaProton() << std::endl
	    << "dEdx: " << dEdx() << " map0: " << map0() << " map1: " << map1() << chi2() << std::endl
	    << "isTofTrack: " << isTofTrack() << " beta: " << beta() << " m^2: " << massSqr() << std::endl
	    << "pMomentum z,y,z: " << pMom().X() << " , " << pMom().Y() << " , " << pMom().Z() << std::endl
	    << "origin x,y,z: " << origin().X() << " , " << origin().Y() << " , " << origin().Z() << std::endl
	    << "gMomentum x,y,z: " << gMom().X() << " , " << gMom().Y() << " , " << gMom().Z() << std::endl
      << std::endl;
}

//_________________
void StFemtoTrack::setNSigmaElectron(const Float_t& ns) {
  // Set nSigma(electron)
  mNSigmaElectron = ( TMath::Abs(ns * 1000.) > std::numeric_limits<short>::max() ?
		      ( (ns>0) ?  std::numeric_limits<short>::max() :
			std::numeric_limits<short>::min() ) :
		      (Short_t)(ns * 1000.) );
}

//_________________
void StFemtoTrack::setNSigmaPion(const Float_t& ns) {
  // Set nSigma(pion)
  mNSigmaPion = ( TMath::Abs(ns * 1000.) > std::numeric_limits<short>::max() ?
		  ( (ns>0) ?  std::numeric_limits<short>::max() :
		    std::numeric_limits<short>::min() ) :
		  (Short_t)(ns * 1000.) );
}

//_________________
void StFemtoTrack::setNSigmaKaon(const Float_t& ns) {
  // Set nSigma(kaon)
  mNSigmaKaon = ( TMath::Abs(ns * 1000.) > std::numeric_limits<short>::max() ?
		  ( (ns>0) ?  std::numeric_limits<short>::max() :
		    std::numeric_limits<short>::min() ) :
		  (Short_t)(ns * 1000.) );
}

//_________________
void StFemtoTrack::setNSigmaProton(const Float_t& ns) {
  // Set nSigma(proton)
  mNSigmaProton = ( TMath::Abs(ns * 1000.) > std::numeric_limits<short>::max() ?
		    ( (ns>0) ?  std::numeric_limits<short>::max() :
		      std::numeric_limits<short>::min() ) :
		    (Short_t)(ns * 1000.) );
}

//_________________
void StFemtoTrack::setChi2(const Float_t& chi2) {
  // Set chi2 value
  if( chi2<=0 ) {
    mChi2 = 0;
  }
  else {
    mChi2 = ( (chi2 * 10.) > std::numeric_limits<unsigned char>::max() ?
	      std::numeric_limits<unsigned char>::max() :
	      (UChar_t)(chi2 * 10.) );
  }
}

//_________________
void StFemtoTrack::setBeta(const Float_t& beta) {
  // Set TOF beta
  if( beta <= 0 ) {
    mBeta = 0;
  }
  else {
    mBeta = ( (beta * 20000.) > std::numeric_limits<unsigned short>::max() ?
	      std::numeric_limits<unsigned short>::max() :
	      (UShort_t)(beta * 20000.) );
  }
}
