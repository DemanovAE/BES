/**
 * \class StFemtoTrack
 * \brief Holds track information
 *
 * The StFemtoTrack class describes track (primary or global). The class
 * is designed in the both standalone and specific environment modes.
 * For the former one needs _VANILLA_ROOT_ global variable
 * to be defined at the compilation and run times.
 *
 * \author Grigory Nigmatkulov; e-mail: nigmatkulov@gmail.com
 * \date June 28, 2018
 */

#ifndef StFemtoTrack_h
#define StFemtoTrack_h

// ROOT headers
#include "TObject.h"
#include "TVector3.h"
#include "TBufferFile.h"
#include "TMath.h"

// FemtoDst headers
#include "StFemtoHelix.h"
#include "StFemtoPhysicalHelix.h"
#ifdef _VANILLA_ROOT_
#include "SystemOfUnits.h"
#include "PhysicalConstants.h"
#else
#include "StarClassLibrary/SystemOfUnits.h"
#include "StarClassLibrary/PhysicalConstants.h"
#endif

//_________________
class StFemtoTrack : public TObject {

 public:
  /// Default constructor
  StFemtoTrack();
  /// Copy constructor
  StFemtoTrack(const StFemtoTrack &copy);
  /// Destructor
  virtual ~StFemtoTrack() { /* empty */}
  /// Print track information
  void print();

  //
  // Setters
  //

  /// Set track unique ID
  void setId(const Short_t& id)        { mId = (UShort_t)id; }
  /// Set number of hits with next encoding: nHits * charge
  void setNHits(const Char_t& nHits)   { mNHits = nHits; }
  /// Set number of hits with next encoding: nHits * charge
  void setNHits(const Short_t& nHits)  { mNHits = (Char_t)nHits; }
  /// Set number of hits with next encoding: nHits * charge
  void setNHits(const Int_t& nHits)    { mNHits = (Char_t)nHits; }
  /// Set possible number of hits in TPC
  void setNHitsPoss(const UShort_t& nHits)  { mNHitsPoss = (UChar_t)nHits; }
  /// Set number of hits that were used in dE/dx estimated
  void setNHitsDedx(const UShort_t& nHits)  { mNHitsDedx = (UChar_t)nHits; }
  /// Set chi2 of the track reconstruction
  void setChi2(const Float_t& chi2);
  /// Set nSigma(electron)
  void setNSigmaElectron(const Float_t& ns);
  /// Set nSigma(pion)
  void setNSigmaPion(const Float_t& ns);
  /// Set nSigma(kaon)
  void setNSigmaKaon(const Float_t& ns);
  /// Set nSigma(proton)
  void setNSigmaProton(const Float_t& ns);
  /// Set dE/dx
  void setDedx(const Float_t& dEdx);
  /// Set dE/dx
  void setDedx(const Double_t& dEdx);
  /// Set the first (out of two) word of topology map
  void setMap0(const Int_t& map0)    { mMap[0] = (UInt_t)map0; }
  /// Set the second (out of two) word of topology map
  void setMap1(const Int_t& map1)    { mMap[1] = (UInt_t)map1; }
  /// Set topology map
  void setTopologyMap(const Int_t& i, const Int_t& map) { if (i==0 || i==1) { mMap[i]=(UInt_t)map; } }
  /// Set X component of primary track momentum
  void setPx(const Float_t& px)      { mPrimMomX = px; }
  /// Set Y component of primary track momentum
  void setPy(const Float_t& py)      { mPrimMomY = py; }
  /// Set Z component of primary track momentum
  void setPz(const Float_t& pz)      { mPrimMomZ = pz; }
  /// Set three-momentum of the primary track (px,py,pz)
  void setP(const Float_t& px, const Float_t& py, const Float_t& pz)
  { setPrimaryMomentum( px, py, pz ); }
  /// Set three-momentum of the primary track (px,py,pz)
  void setPrimaryMomentum(const Float_t& px, const Float_t& py, const Float_t& pz)
  {mPrimMomX = px; mPrimMomY = py; mPrimMomZ = pz; }
  /// Set position (x,y,z) of the track at DCA to primary vertex
  void setOriginXYZ(const Float_t& x, const Float_t& y, const Float_t& z)
  { mOriginX = x; mOriginY = y; mOriginZ = z; }
  /// Set X position of the track origin at DCA to primary vertex
  void setOriginX(const Float_t& x)  { mOriginX = x; }
  /// Set Y position of the track origin at DCA to primary vertex
  void setOriginY(const Float_t& y)  { mOriginY = y; }
  /// Set Z position of the track origin at DCA to primary vertex
  void setOriginZ(const Float_t& z)  { mOriginZ = z; }
  /// Set momentum of the global track (px,py,pz) at DCA to primary vertex
  void setGlobalP(const Float_t& px, const Float_t& py, const Float_t& pz)
  { mGlobMomX = px; mGlobMomY = py; mGlobMomZ = pz; }
  /// Set TOF beta
  void setBeta(const Float_t& beta);

  //
  // Getters
  //

  /// Return track unique ID
  Short_t   id() const              { return (Short_t)mId; }
  /// Return track type: 0 - global; 1 - primary
  Short_t   type() const            { return ( isPrimary() ) ? 1 : 0; }
  /// Check if the track is primary
  Bool_t    isPrimary() const       { return (pMom().Mag()>0) ? true : false; }
  /// Return track flag. Not store it but 300 represents
  /// that track was reconstructed in TPC. Need it for historical reasons
  Short_t   flag() const            { return 300; }
  /// Return number of hits
  UShort_t  nHits() const           { return (UShort_t)TMath::Abs(mNHits); }
  /// Return number of hits used in a track fit
  UShort_t  nHitsFit() const        { return nHits(); }
  /// Return possible number of hits in TPC
  UShort_t  nHitsPoss() const       { return (UShort_t)mNHitsPoss; }
  /// Return number of hits used in dE/dx estimation
  UShort_t  nHitsDedx() const       { return (UShort_t)mNHitsDedx; }
  /// Return nSigma(electron)
  Float_t   nSigmaElectron() const  { return (Float_t)mNSigmaElectron / 1000.; }
  /// Return nSigma(pion)
  Float_t   nSigmaPion() const      { return (Float_t)mNSigmaPion / 1000.; }
  /// Return nSigma(kaon)
  Float_t   nSigmaKaon() const      { return (Float_t)mNSigmaKaon / 1000.; }
  /// Return nSigma(proton)
  Float_t   nSigmaProton() const    { return (Float_t)mNSigmaProton / 1000.; }
  /// Dummy method that returns probability to be an electron
  Float_t   pidProbElectron() const { return 0.5; }
  /// Dummy method that returns probability to be an pion
  Float_t   pidProbPion() const     { return 0.5; }
  /// Dummy method that returns probability to be an kaon
  Float_t   pidProbKaon() const     { return 0.5; }
  /// Dummy method that returns probability to be an proton
  Float_t   pidProbProton() const   { return 0.5; }
  /// Return dE/dx of the track in [GeV * cm]
  Double_t  dEdx() const            { return (Double_t)mDedx * 1e-9; }
  /// Return the first word of topology map
  UInt_t    map0() const            { return (UInt_t)mMap[0]; }
  /// Return the second word of topology map
  UInt_t    map1() const            { return (UInt_t)mMap[1]; }
  /// Return the word of topology map
  UInt_t    map(const Int_t& word)  { return (word==0) ? mMap[0] : mMap[1]; }
  /// Return chi2 of the track
  Float_t   chi2() const            { return (Float_t)mChi2 / 10.; }
  /// Return charge of the track
  Short_t   charge() const          { return (mNHits>0) ? +1 : -1; }
  /// Return origin position (x,y,z) of the global track at DCA to the primary vertex
  TVector3  origin() const          { return TVector3( mOriginX, mOriginY, mOriginZ); }
  /// Return transverse distance between (pointX,pointY) and (originX,originY)
  Float_t   gDCAxy(const Float_t& pointX, const Float_t& pointY) const;
  /// Return longitudinal distance between (pointZ) and (originZ)
  Float_t   gDCAz(const Float_t& pointZ) const;
  /// Return distance between point (pointX,pointY,pointZ) and origin (originX,originY,originZ)
  Float_t   gDCA(const Float_t& pVtxX, const Float_t& pVtxY, const Float_t& pVtxZ) const;
  /// Return 3D distance between point (pointX,pointY,pointZ) and origin (originX,originY,originZ)
  TVector3  gDCA(TVector3 point) const;

  /// Return 3-momentum of the primary track. (0,0,0) if the track is global
  TVector3  pMom() const       { return TVector3(mPrimMomX, mPrimMomY, mPrimMomZ); }
  /// Return magnitude of primary track 3-momentum. 0 if the track is global
  Float_t   p() const          { return ptot(); }
  /// Return magnitude of primary track 3-momentum. 0 if the track is global
  Float_t   ptot() const       { return ( isPrimary() ) ? pMom().Mag() : 0.; }
  /// Return squared magnitude of primary track 3-momentum. 0 if the track is global
  Float_t   ptot2() const      { return ( isPrimary() ) ? pMom().Mag2() : 0.; }
  /// Return transverse momentum of primary track. 0 if the track is global
  Float_t   pt() const         { return ( isPrimary() ) ? pMom().Perp() : 0.; }
  /// Return squared transverse momentum of primary track. 0 if the track is global
  Float_t   pt2() const        { return ( isPrimary() ) ? pMom().Perp2() : 0.; }
  /// Return pseudorapidity of primary track. 0 if the track is global
  Float_t   eta() const        { return ( isPrimary() ) ? pMom().PseudoRapidity() : 0.; }
  /// Return azimuthal angle of primary track. 0 if the track is global
  Float_t   phi() const        { return ( isPrimary() ) ? pMom().Phi() : 0.; }

  /// Return the 3-momentum of global track at the dca point to the
  /// point (usually it is the primary vertex). B - magnetic field
  /// from FemtoEvent::bField()
  TVector3  gMom() const       { return TVector3(mGlobMomX, mGlobMomY, mGlobMomZ); }
  /// Return magnitude of global track 3-momentum
  Float_t   gP() const         { return gPtot(); }
  /// Return magnitude of global track 3-momentum
  Float_t   gPtot() const      { return gMom().Mag(); }
  /// Return squared magnitude of global track 3-momentum
  Float_t   gPtot2() const     { return gMom().Mag2(); }
  /// Return transverse momentum of global track
  Float_t   gPt() const        { return gMom().Perp(); }
  /// Return squared transverse momentum of global track
  Float_t   gPt2() const       { return gMom().Perp2(); }
  /// Return pseudorapidity of global track
  Float_t   gEta() const       { return gMom().PseudoRapidity(); }
  /// Return beta measured in TOF
  Float_t   beta() const       { return ( isTofTrack() ) ? (Float_t)mBeta / 20000. : -999.; }
  /// Return 1./beta
  Float_t   invBeta() const    { return ( isTofTrack() ) ? 1./beta() : -999.; }
  /// Return 1./beta^2
  Float_t   invBeta2() const   { return ( isTofTrack() ) ? invBeta()*invBeta() : -999.; }
  /// Return squared mass estimated by TOF beta
  Float_t   massSqr() const;
  /// Check if track had signal in TOF
  Bool_t    isTofTrack() const { return (mBeta>0) ? true : false; }

  /// Global momentum at point of DCA to pVtx, B should be in kilogauss
  TVector3 gMom(const TVector3 pVtx, const Float_t B) const;
  /// Return helix of the global track
  StFemtoPhysicalHelix helix(const Float_t& B) const;

 private:

  /// Track unique ID
  UShort_t  mId;
  /// Number of hits * charge
  Char_t    mNHits;
  /// Number of possible hits
  UChar_t   mNHitsPoss;
  /// Number of hits used for dE/dx estimation
  UChar_t   mNHitsDedx;
  /// dEdx * 10^9 (in eV/cm)
  UShort_t mDedx;
  /// Topology map (word0 and word1)
  UInt_t mMap[2];
  /// chi2 of the track reconstruction (covnersion factor = *10)
  UChar_t mChi2;
  /// nSigma(electron) via dE/dx (conversion factor = *1000)
  Short_t mNSigmaElectron;
  /// nSigma(pion) via dE/dx (conversion factor = *1000)
  Short_t mNSigmaPion;
  /// nSigma(kaon) via dE/dx (conversion factor = *1000)
  Short_t mNSigmaKaon;
  /// nSigma(proton) via dE/dx (conversion factor = *1000)
  Short_t mNSigmaProton;
  /// Beta in TOF (0 if no signal) (conversion factor = *20000)
  UShort_t mBeta;

  /// X component of track momentum fitted to the primary vertex
  Float_t mPrimMomX;
  /// Y component of track momentum fitted to the primary vertex
  Float_t mPrimMomY;
  /// Z component of track momentum fitted to the primary vertex
  Float_t mPrimMomZ;

  /// X component of global track momentum at DCA to primary vertex
  Float_t mGlobMomX;
  /// Y component of global track momentum at DCA to primary vertex
  Float_t mGlobMomY;
  /// Z component of global track momentum at DCA to primary vertex
  Float_t mGlobMomZ;
  /// X position of global track at DCA to primary vertex
  Float_t mOriginX;
  /// Y position of global track at DCA to primary vertex
  Float_t mOriginY;
  /// Z position of global track at DCA to primary vertex
  Float_t mOriginZ;

  ClassDef(StFemtoTrack, 2)
};

//_________________
inline TVector3 StFemtoTrack::gMom(const TVector3 pVtx, const Float_t B) const {
  StFemtoPhysicalHelix gHelix = helix(B);
  return gHelix.momentumAt( gHelix.pathLength(pVtx), B * kilogauss );
}

//_________________
inline StFemtoPhysicalHelix StFemtoTrack::helix(const Float_t& B) const {
  return StFemtoPhysicalHelix( TVector3(mGlobMomX, mGlobMomY, mGlobMomZ),
			       TVector3(mOriginX, mOriginY, mOriginZ),
			       B * kilogauss,
			       static_cast<float>( charge() ) );
}

#endif // StFemtoTrack_h
