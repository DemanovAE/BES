/**
 * \class StFemtoXi
 * \brief The class keeps pointers to the V0 and track
 * in the StFemtoV0 and StFemtoTrack arrays
 *
 * The StFemtoXi holds indices of V0 and track in the
 * StFemtoV0 and StFemtoTrack arrays
 *
 * \author Egor Alpatov, Grigory Nigmatkulov; e-mail: nigmatkulov@gmail.com
 * \date November 20, 2018
 */

#ifndef StFemtoXi_h
#define StFemtoXi_h

// ROOT headers
#include "TObject.h"

//________________
class StFemtoXi : public TObject {

 public:
  /// Default costructor
  StFemtoXi();
  /// Copy constructor
  StFemtoXi(const StFemtoXi &xi);
  /// Constructor that takes indices of V0 and track
  StFemtoXi(const Int_t& v0Id, const Int_t& trkId)
  { mV0Id = v0Id; mTrackId = trkId; }
  /// Destructor
  virtual ~StFemtoXi()                  { /* empty */}
  /// Print Xi information
  void print();

  //
  // Setters
  //

  void setV0Id(const Int_t& v0Id)       { mV0Id = v0Id; }
  void setTrackId(const Int_t& trackId) { mTrackId = trackId; }

  //
  // Getters
  //

  Int_t v0Id() const                    { return mV0Id; }
  Int_t trackId() const                 { return mTrackId; }

 private:

  /// Index of Lambda in the StFemtoV0 array
  UShort_t mV0Id;
  /// Index of bachelor in the StFemtoTrack array
  UShort_t mTrackId;

  ClassDef(StFemtoXi, 1);
};

#endif // #define StFemtoXi_h
