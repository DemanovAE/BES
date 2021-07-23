/**
 * \class StFemtoDst
 * \brief The class holds FemtoDst structure
 *
 * The class keeps the structure of FemtoDst and allows
 * one to get event, track, V0 or Xi information
 *
 * \author Grigory Nigmatkulov; e-mail: nigmatkulov@gmail.com
 * \date July 15, 2018
 */

#ifndef StFemtoDst_h
#define StFemtoDst_h

// ROOT headers
#include "TClonesArray.h"

// FemtoDst headers
#include "StFemtoArrays.h"

// Forward declarations
class StFemtoEvent;
class StFemtoTrack;
class StFemtoEpdHit;
class StFemtoV0;
class StFemtoXi;

//_________________
class StFemtoDst {

 public:
  /// Default constructor
  StFemtoDst()  {/* emtpy */}
  /// Destructor
  ~StFemtoDst() {/* emtpy*/}

  /// Set the pointers to the TClonesArrays
  static void set(TClonesArray** array);
  /// Reset the pointers to the TClonesArrays to 0
  static void unset();
  /// Return pointer to the n-th TClonesArray
  static TClonesArray* femtoArray(Int_t type) { return femtoArrays[type]; }
  /// Return pointer to current StFemtoEvent (class holding the event wise information)
  static StFemtoEvent* event() { return (StFemtoEvent*)femtoArrays[StFemtoArrays::Event]->UncheckedAt(0); }
  /// Return pointer to the i-th track
  static StFemtoTrack* track(Int_t i) { return (StFemtoTrack*)femtoArrays[StFemtoArrays::Track]->UncheckedAt(i); }
  /// Return pointer to the i-th V0
  static StFemtoV0* v0(Int_t i) { return (StFemtoV0*)femtoArrays[StFemtoArrays::V0]->UncheckedAt(i); }
  /// REturn pointer to the i-th Xi
  static StFemtoXi* xi(Int_t i) { return (StFemtoXi*)femtoArrays[StFemtoArrays::Xi]->UncheckedAt(i); }
  /// Return pointer to i-th epd hit
  static StFemtoEpdHit* epdHit(Int_t i) { return (StFemtoEpdHit*)femtoArrays[StFemtoArrays::EpdHit]->UncheckedAt(i); }

  /// Return number of tracks in the femto arrays
  static UInt_t numberOfTracks() { return femtoArrays[StFemtoArrays::Track]->GetEntriesFast(); }
  /// Return number of V0s in the femto femto array
  static UInt_t numberOfV0s() { return femtoArrays[StFemtoArrays::V0]->GetEntriesFast(); }
  /// Return number of Xis in the femto array
  static UInt_t numberOfXis() { return femtoArrays[StFemtoArrays::Xi]->GetEntriesFast(); }
  /// Return number of EPD hits
  static UInt_t numberOfEpdHits() { return femtoArrays[StFemtoArrays::EpdHit]->GetEntriesFast(); }

  /// Print information
  void print() const;
  /// Print track information
  static void printTracks();
  /// Print V0 information
  static void printV0s();
  /// Print Xi information
  static void printXis();

 private:

  /// Array of TClonesArrays
  static TClonesArray** femtoArrays;
};

#endif // StFemtoDst_h
