/**
 * \class StFemtoV0
 * \brief Keeps V0 daugheter IDs
 *
 * The class keeps indices of the positive and negative
 * daughters from StFemtoTrack StFemtoTrack array
 *
 * \author Egor Alpatov, Grigory Nigmatkulov; e-mail: nigmatkulov@gmail.com
 * \date November 20, 2018
 */

#ifndef StFemtoV0_h
#define StFemtoV0_h

// ROOT headers
#include "TObject.h"

//________________
class StFemtoV0 : public TObject {

 public:
  /// Default costructor
  StFemtoV0();
  /// Copy constructor
  StFemtoV0(const StFemtoV0 &v0);
  /// Constructor that takes indices of positive and negative daughters
  /// in the StFemtoTrack array
  StFemtoV0(const Int_t& posId, const Int_t& negId)
  { mPosId = posId; mNegId = negId; }
  /// Destructor
  virtual ~StFemtoV0()            { /* empty */ }
  /// Print V0 information
  void print();

  //
  // Setters
  //

  /// Set index of positive track in StFemtoTrack array
  void setPosId(const Int_t& pos) { mPosId = pos; }
  /// Set index of negative track in StFemtoTrack array
  void setNegId(const Int_t& neg) { mNegId = neg; }

  //
  // Getters
  //

  /// Return index of positive track in StFemtoTrack array
  Int_t posDaughterId() const     { return mPosId; }
  /// Return index of negative track in StFemtoTrack array
  Int_t negDaughterId() const     { return mNegId; }

 private:

	/// Index of positive track in StFemtoTrack array
  UShort_t mPosId;
  /// Index of negative track in StFemtoTrack array
  UShort_t mNegId;

  ClassDef(StFemtoV0, 1);
};

#endif // #define StFemtoV0_h
