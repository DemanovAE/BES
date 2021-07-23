/**
 * \class StFemtoEvent
 * \brief Holds parameters of the current event
 *
 * The StFemtoEvent class stores parameters of event. In order to work
 * in a standalone mode one should define _VANILLA_ROOT_ global variable
 * at the compilation and run times.
 *
 * \author Grigory Nigmatklov, Nikita Ermakov; e-mail: nigmatkulov@gmail.com
 * \date June 28, 2018
 */

#ifndef StFemtoEvent_h
#define StFemtoEvent_h

// C++ headers
#include <vector>
#include <limits>

// ROOT headers
#include "TObject.h"
#include "TVector3.h"

//_________________
class StFemtoEvent : public TObject {

 public:
  /// Default constructor
  StFemtoEvent();
  /// Copy constructor
  StFemtoEvent(const StFemtoEvent &copy);
  /// Destructor
  virtual ~StFemtoEvent()    { /* emtpy */ }
  /// Print some event information
  void print();

  //
  // Setters
  //

  /// Set event ID
  void setEventId(const Int_t& id)       { mEventId = id; }
  /// Set run ID
  void setRunId(const Int_t& id)         { mRunId = id; }
  /// Set fill ID
  void setFillId(Int_t id)               { mFillId = (UShort_t)id; }
  /// Set fill ID
  void setFillId(Float_t id)             { mFillId = (id > 0) ? (UShort_t)id : 0; }
  /// Set time stamp
  void setTime(Int_t time)               { mTime = time; }
  /// Set refMult
  void setRefMult(const Int_t& m);
  /// Set refMult2
  void setRefMult2(const Int_t& m);
  /// Set refMultCorr(ected) from StRefMult
  void setRefMultCorr(const Float_t& rf);
  /// Set refMultCorrWeight from StRefMult
  void setRefMultCorrWeight(const Float_t& w);
  /// Set refMult of positive particles
  void setRefMultPos(const Int_t& m);
  /// Set refMult2 of positive particles
  void setRefMult2Pos(const Int_t& m);
  /// Set gRefMult
  void setGRefMult(const Int_t& m);
  /// Set cent16
  void setCent16(const Int_t& c)        { mCent16 = (Char_t)c; }
  /// Set number of bTofHits (number of hits in TOF trays)
  void setNumberOfBTofHit(const Int_t& m);
  /// Set number of tracks that matched TOF
  void setNumberOfTofMatched(const Int_t& m);
  /// Set number of tracks that matched BEMC
  void setNumberOfBEMCMatched(const Int_t& m);
  /// Set number of primary tracks in event
  void setNumberOfPrimaryTracks(const Int_t& m);
  /// Set number of global tracks in event
  void setNumberOfGlobalTracks(const Int_t& m);
  /// Set east ZDC ADC sum
  void setZdcSumAdcEast(const Float_t& x);
  /// Set west ZDC ADC sum
  void setZdcSumAdcWest(const Float_t& x);
  /// Set ADC of east ZDC horizontal i-th strip
  void setZdcSmdEastHorizontal(const Int_t& i, const Float_t& x);
  /// Set ADC of east ZDC vertical i-th strip
  void setZdcSmdEastVertical(const Int_t& i, const Float_t& x);
  /// Set ADC of west ZDC horizontal i-th strip
  void setZdcSmdWestHorizontal(const Int_t& i, const Float_t& x);
  /// Set ADC of west ZDC vertical i-th strip
  void setZdcSmdWestVertical(const Int_t& i, const Float_t& x);
  /// Set ADC of east BBC i-th PMT
  void setBbcAdcEast(const Int_t& i, const UShort_t& x);
  /// Set ADC of west BBC i-th PMT
  void setBbcAdcWest(const Int_t& i, const UShort_t& x);
  /// Set magnetic field
  void setMagneticField(const Float_t& f)  { mMagField = f; }
  /// Set magnetic field
  void setBField(const Float_t& f)         { setMagneticField( f ); }
  /// Set primary vertex position (x,y,z)
  void setPrimaryVertexPosition(const Float_t& x, const Float_t& y, const Float_t& z);
  /// Set x position of primary vertex
  void setPrimaryVertexX(const Float_t& x) { mVertexPositionX = x; }
  /// Set y position of primary vertex
  void setPrimaryVertexY(const Float_t& y) { mVertexPositionY = y; }
  /// Set z position of primary vertex
  void setPrimaryVertexZ(const Float_t& z) { mVertexPositionZ = z; }
  /// Set number of primary vertices in event
  void setNumberOfPrimaryVertices(const Int_t& mult);
  /// Set transverse sphericity |eta|<0.5
  void setTransverseSphericity(const Float_t& sph);
  /// Set transverse sphericity |eta|<1
  void setTransverseSphericity2(const Float_t& sph);
  /// Set z position of primary vertex estimated by VPD
  void setVpdVz(const Float_t& vz)         { mVpdVz = vz; }
  /// Set primary vertex ranking
  void setPrimaryVertexRanking(const Float_t& ranking);
  /// Add trigger
  void addTriggerId(const UInt_t& id);
  /// Add trigger
  void setTriggerId(const UInt_t& id);
  /// Add list (std::vector) of triggers
  void setTriggerIds(std::vector<unsigned int> ids) { mTriggerIds = ids; }

  /// Set ZDC coincidence rate
  void setZDCx(Float_t zdcCoinRate)         { mZDCx = (UInt_t)zdcCoinRate; }
  /// Set BBC coincidence rate
  void setBBCx(Float_t bbcCoinRate)         { mBBCx = (UInt_t)bbcCoinRate; }
  /// Set background rate
  void setBackgroundRate(Float_t bckgRate)  { mBackgroundRate = (Float_t)bckgRate; }
  /// Set "blue"-beam background rate
  void setBbcBlueBackgroundRate(Float_t bbcBlueBckgRate)
  { mBbcBlueBackgroundRate = (Float_t)bbcBlueBckgRate; }
  /// Set "yellow"-beam background rate
  void setBbcYellowBackgroundRate(Float_t bbcYellowBckgRate)
  { mBbcYellowBackgroundRate = (Float_t)bbcYellowBckgRate; }
  /// Set east BBC rate
  void setBbcEastRate(Float_t bbcEastRate)  { mBbcEastRate = (Float_t)bbcEastRate; }
  /// Set west BBC rate
  void setBbcWestRate(Float_t bbcWestRate)  { mBbcWestRate = (Float_t)bbcWestRate; }
  /// Set east ZDC rate
  void setZdcEastRate(Float_t zdcEastRate)  { mZdcEastRate = (Float_t)zdcEastRate; }
  /// Set west ZDC rate
  void setZdcWestRate(Float_t zdcWestRate)  { mZdcWestRate = (Float_t)zdcWestRate; }

  //
  // Getters
  //

  /// Return event ID
  Int_t    eventId() const          { return mEventId; }
  /// Return run ID
  Int_t    runId() const            { return mRunId; }
  /// Return fill ID
  Int_t    fillId() const           { return (Int_t)mFillId; }
  /// Return time stamp
  Int_t    time() const             { return mTime; }
  /// Return year
  Int_t    year() const             { return ( (mRunId / 1000000) - 1 + 2000 ); }
  /// Return day
  Int_t    day() const              { return ( (mRunId % 1000000) / 1000 ); }
  /// Return refMult
  Int_t    refMult() const          { return (Int_t)mRefMult; }
  /// Return refMult2
  Int_t    refMult2() const         { return (Int_t)mRefMult2;}
  /// Return refMultCorr(ected) from StRefMultCorr
  Float_t  refMultCorr() const      { return (Float_t)mRefMultCorr * 0.1; }
  /// Return refMultCorrWeight from StRefMultCorr
  Float_t  refMultCorrWeight() const{ return (Float_t)mRefMultCorrWeight * 0.0001; }
  /// Return gRefMult (refMult estimated by global tracks)
  UInt_t   gRefMult() const         { return (UInt_t)mGRefMult; }
  /// Return centrality for 9 bins:
  /// 8 = 0-5%
  /// 7 = 5-10%
  /// 6 = 10-20%
  /// 5 = 20-30%
  /// 4 = 30-40%
  /// 3 = 40-50%
  /// 2 = 50-60%
  /// 1 = 60-70%
  /// 0 = 70-80%
  Int_t    cent9() const;
  /// Return centrality for 16 bins:
  /// 15 = 0-5%
  /// 14 = 5-10%
  /// 13 = 10-15%
  /// 12 = 15-20%
  /// 11 = 20-25%
  /// 10 = 25-30%
  ///  9 = 30-35%
  ///  8 = 35-40%
  ///  7 = 40-45%
  ///  6 = 45-50%
  ///  5 = 50-55%
  ///  4 = 55-60%
  ///  3 = 60-65%
  ///  2 = 65-70%
  ///  1 = 70-75%
  ///  0 = 75-80%
  Int_t    cent16() const           { return (Int_t)mCent16; }
  /// Return refMult estimated by positive tracks
  Int_t    refMultPos() const       { return (Int_t)mRefMultPos; }
  /// Return refMult estimated by negative tracks
  Int_t    refMultNeg() const       { return ( refMult() - refMultPos() ); }
  /// Return refMult2 estimated by positive tracks
  Int_t    refMult2Pos() const      { return (Int_t)mRefMult2Pos; }
  /// Return refMult2 estimated by negative tracks
  Int_t    refMult2Neg() const      { return ( refMult2() - refMult2Pos() ); }

  /// Return sum of ADC for east ZDC
  Float_t  zdcSumAdcEast() const    { return (Float_t)mZdcSumAdcEast; }
  /// Return sum of ADC for west ZDC
  Float_t  zdcSumAdcWest() const    { return (Float_t)mZdcSumAdcWest; }
  /// Return ADC for i-th horizontal strip in east ZDC
  Float_t  zdcSmdEastHorizontal(const Int_t& i) const
  { return (Float_t)mZdcSmdEastHorizontal[i]; }
  /// Return ADC for i-th vertical strip in east ZDC
  Float_t  zdcSmdEastVertical(const Int_t& i) const
  { return (Float_t)mZdcSmdEastVertical[i]; }
  /// Return ADC for i-th horizontal strip in west ZDC
  Float_t  zdcSmdWestHorizontal(const Int_t& i) const
  { return (Float_t)mZdcSmdWestHorizontal[i]; }
  /// Return ADC for i-th vertical strip in west ZDC
  Float_t  zdcSmdWestVertical(const Int_t& i) const
  { return (Float_t)mZdcSmdWestVertical[i]; }
  /// Return ADC for i-th PMT box of east BBC
  UShort_t bbcAdcEast(const Int_t& i) const
  { return mBbcAdcEast[i]; }
  /// Return ADC for i-th PMT box of west BBC
  UShort_t bbcAdcWest(const Int_t& i) const
  { return mBbcAdcWest[i]; }
  /// Return number of hits in TOF (in trays)
  UShort_t numberOfBTofHit() const       { return mNumberOfBTofHit; }
  /// Return number of primary tracks
  UShort_t numberOfPrimaryTracks() const { return mNumberOfPrimaryTracks; }
  /// Return number of global tracks
  UShort_t numberOfGlobalTracks() const  { return mNumberOfGlobalTracks; }
  /// Return number of tracks that matched TOF
  UShort_t numberOfTofMatched() const    { return mNumberOfTofMatched; }
  /// Return number of tracks that matched BEMC
  UShort_t numberOfBEMCMatched() const   { return mNumberOfBEMCMatched; }
  /// Return magnetic field
  Float_t  magneticField() const         { return mMagField; }
  /// Return postion of primary vertex (x,y,z)
  TVector3 primaryVertex() const { return TVector3(mVertexPositionX, mVertexPositionY, mVertexPositionZ); }
  /// Return postion of primary vertex (x,y,z)
  TVector3 vertexPosition() const        { return primaryVertex(); }
  /// Return number of reconstructed primary vertices in event
  UInt_t   numberOfPrimaryVertices() const { return (UInt_t)mNumberOfPrimaryVertices; }
  /// Return transverse sphericity |eta|<0.5
  Float_t  transverseSphericity() const  { return (Float_t)mTransverseSphericity * 0.01; }
  /// Return transverse sphericity |eta|<1
  Float_t  transverseSphericity2() const { return (Float_t)mTransverseSphericity2 * 0.01; }
  /// Return radial position of primary vertex
  Float_t  vertexPositionR() const       { return primaryVertex().Perp(); }
  /// Return z postion of primary vertex estimated by VPD
  Float_t  vpdVz() const                 { return mVpdVz; }
  /// Return list (std::vector) of trigger IDs
  std::vector<unsigned int> triggerIds() const { return mTriggerIds; }
  /// Check if trigger ID is in trigger list
  Bool_t   isTrigger(const unsigned int& id) const;
  /// Return primary vertex ranking. FemtoDst encoding:
  /// -10 = x < -10e9
  ///  -9 = -10e9 <= x < -10e6
  ///  -8 = -10e6 <= x < -10e3;
  ///  -7 = -10e3 <= x < -10e2
  ///  -6 = -10e2 <= x < -10
  ///  -5 = -10 <= x < -5
  ///  -4 = -5 <= x < -4
  ///  -3 = -4 <= x < -3
  ///  -2 = -3 <= x < -2
  ///  -1 = -2 <= x < -1
  ///   0 = -1 <= x < 1
  ///   1 = 1 <= x < 2
  ///   2 = 2 <= x < 3
  ///   3 = 3 <= x < 4
  ///   4 = 4 <= x < 5
  ///   5 = 5 <= x < 10
  ///   6 = 10 <= x < 10e2
  ///   7 = 10e2 <= x < 10e3
  ///   8 = 10e3 <= x < 10e6
  ///   9 = 10e6 <= x < 10e9
  ///  10 = 10e9 <= x
  Float_t  ranking() const
  { return (Float_t)mPrimaryVertexRanking; }

  /// Return ZDC coincidence rate
  Float_t  ZDCx() const                  { return (Float_t)mZDCx; }
  /// Return BBC coincidence rate
  Float_t  BBCx() const                  { return (Float_t)mBBCx; }
  /// Return background rate
  Float_t  backgroundRate() const        { return mBackgroundRate; }
  /// Return "blue"-beam background rate
  Float_t  bbcBlueBackgroundRate() const { return mBbcBlueBackgroundRate; }
  /// Return "yellow"-beam background rate
  Float_t  bbcYellowBackgroundRate() const { return mBbcYellowBackgroundRate; }
  /// Return east BBC rate
  Float_t  bbcEastRate() const           { return mBbcEastRate; }
  /// Return west BBC rate
  Float_t  bbcWestRate() const           { return mBbcWestRate; }
  /// Return east ZDC rate
  Float_t  zdcEastRate() const           { return mZdcEastRate; }
  /// Return west ZDC rate
  Float_t  zdcWestRate() const           { return mZdcWestRate; }

 private:

  /// Event ID
  Int_t    mEventId;
  /// Run ID
  Int_t    mRunId;

  /// refMult
  UShort_t mRefMult;
  /// mRefMultCorr * 10
  UShort_t mRefMultCorr;
  /// mRefMultCorrWight * 10000
  UShort_t mRefMultCorrWeight;
  // n positive for refMult
  UShort_t mRefMultPos;
  /// refMult for |eta| < 1.0
  UShort_t mRefMult2;
  /// n positive for refMult2
  UShort_t mRefMult2Pos;
  /// mGRefMult
  UShort_t mGRefMult;
  /// Centrality for 16 bins (-1 if not defined)
  Char_t   mCent16;

  /// Number of primary tracks in muEvent
  UShort_t mNumberOfPrimaryTracks;
  /// Number of global tracks in muEvent
  UShort_t mNumberOfGlobalTracks;
  /// btofTrayMultiplicity
  UShort_t mNumberOfBTofHit;
  /// Number of TOF matched tracks
  UShort_t mNumberOfTofMatched;
  /// Number of BEMC matched tracks
  UShort_t mNumberOfBEMCMatched;

  /// Magnetic field
  Float_t  mMagField;
  /// Primary vertex position
  Float_t  mVertexPositionX, mVertexPositionY, mVertexPositionZ;
  /// VpdVz
  Float_t  mVpdVz;

  /// Primary vertex ranking according to detector-based weighting procedure. FemtoDst encoding:
  /// -10 = x < -10e9
  ///  -9 = -10e9 <= x < -10e6
  ///  -8 = -10e6 <= x < -10e3;
  ///  -7 = -10e3 <= x < -10e2
  ///  -6 = -10e2 <= x < -10
  ///  -5 = -10 <= x < -5
  ///  -4 = -5 <= x < -4
  ///  -3 = -4 <= x < -3
  ///  -2 = -3 <= x < -2
  ///  -1 = -2 <= x < -1
  ///   0 = -1 <= x < 1
  ///   1 = 1 <= x < 2
  ///   2 = 2 <= x < 3
  ///   3 = 3 <= x < 4
  ///   4 = 4 <= x < 5
  ///   5 = 5 <= x < 10
  ///   6 = 10 <= x < 10e2
  ///   7 = 10e2 <= x < 10e3
  ///   8 = 10e3 <= x < 10e6
  ///   9 = 10e6 <= x < 10e9
  ///  10 = 10e9 <= x
  Char_t   mPrimaryVertexRanking;
  /// Number of primary vertices in the muEvent
  UChar_t  mNumberOfPrimaryVertices;
  /// Transverse sphericity: pTrk with p>0.15 GeV, |eta|<0.5, nHits>10
  Char_t  mTransverseSphericity;
  /// Transverse sphericity: pTrk with p>0.15 GeV, |eta|<1, nHits>10
  Char_t  mTransverseSphericity2;

  /// List (std::vector) of triggers that were fired
  /// for the current event
  std::vector<unsigned int> mTriggerIds;

  /// Sum of east ZDC ADC
  UShort_t mZdcSumAdcEast;
  /// Sum of west ZDC ADC
  UShort_t mZdcSumAdcWest;
  /// ADC of east ZDC i-th horizontal strip
  UShort_t mZdcSmdEastHorizontal[8];
  /// ADC of east ZDC i-th vertical strip
  UShort_t mZdcSmdEastVertical[8];
  /// ADC of west ZDC i-th horizontal strip
  UShort_t mZdcSmdWestHorizontal[8];
  /// ADC of west ZDC i-th vertical strip
  UShort_t mZdcSmdWestVertical[8];

  /// ADC of east BBC i-th PMT box
  UShort_t mBbcAdcEast[24];
  /// ADC of west BBC i-th PMT box
  UShort_t mBbcAdcWest[24];

  /// Fill number
  UShort_t mFillId;
  /// To set timestamp for St_db_Maker. This is what StMuDstMaker used to
  /// GetEvtHddr()->SetGMTime(cast into unsigned int).
  Int_t    mTime;

  /// ZDC coincidende rate
  UInt_t   mZDCx;
  /// BBC coincidende rate
  UInt_t   mBBCx;
  /// Background rate
  Float_t  mBackgroundRate;
  /// "Blue" beam rate measured in BBC
  Float_t  mBbcBlueBackgroundRate;
  /// "Yellow" beam rate measured in BBC
  Float_t  mBbcYellowBackgroundRate;
  /// East BBC rate
  Float_t  mBbcEastRate;
  /// West BBC rate
  Float_t  mBbcWestRate;
  /// East ZDC rate
  Float_t  mZdcEastRate;
  /// West ZDC rate
  Float_t  mZdcWestRate;

  ClassDef(StFemtoEvent, 3)
};

#endif // StFemtoEvent_h
