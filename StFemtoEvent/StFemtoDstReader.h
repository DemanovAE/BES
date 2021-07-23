/**
 * \class StFemtoDstReader
 * \brief Allows to read femtoDst file or a list of files
 *
 * This class allows to read femtoDst.root file or a list of files
 * that contain femtoDst and sets up pointers to the femtoDst, and
 * certain TClonesArrays that keep Event, Track, etc...
 * One can also turn on or off certain branches using the
 * SetStatus method.
 *
 * \author Grigory Nigmatkulov; e-mail: nigmatkulov@gmail.com
 * \date May 28, 2018
 */

#ifndef StFemtoDstReader_h
#define StFemtoDstReader_h

// ROOT headers
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TClonesArray.h"

// FemtoDst headers
#include "StFemtoDst.h"
#include "StFemtoEvent.h"
#include "StFemtoArrays.h"

//_________________
class StFemtoDstReader : public TObject {

 public:
  /// Constructor that takes either femtoDst file or file that
  /// contains a list of femtoDst.root files
  StFemtoDstReader(const Char_t* inFileName);
  /// Destructor
  ~StFemtoDstReader();

  /// Return a pointer to femtoDst (return NULL if no dst is found)
  StFemtoDst *femtoDst()    { return mFemtoDst; }
  /// Return pointer to the chain of .femtoDst.root files
  TChain *chain()           { return mChain; }
  /// Return pointer to the current TTree
  TTree *tree()             { return mTree; }

  /// Set enable/disable branch matching when reading femtoDst
  void SetStatus(const Char_t* branchNameRegex, Int_t enable);

  /// Calls openRead()
  void Init();
  /// Read next event in the chain
  Bool_t readFemtoEvent(Long64_t iEvent);
  /// Close files and finilize
  void Finish();

 private:

  /// Name of the inputfile (or of the inputfiles.list)
  TString mInputFileName;

  /// Turn off streamers
  void streamerOff();

  /// Create TClonesArray of femto classes and set them to femtoDst
  void createArrays();
  /// Clear all TClonesArrays with femtoDst classes
  void clearArrays();
  /// Set adresses of femtoArrays and their statuses (enable/disable) to chain
  void setBranchAddresses(TChain *chain);

  /// Pointer to the input/output femtoDst structure
  StFemtoDst *mFemtoDst;
  /// Pointer to the chain
  TChain *mChain;
  /// Pointer to the current tree
  TTree *mTree;

  /// Event counter
  Int_t mEventCounter;

  /// Pointer to the TClonesArray with the data
  TClonesArray *mFemtoArrays[StFemtoArrays::NAllFemtoArrays];
  /// Status of the TClonesArray
  Char_t        mStatusArrays[StFemtoArrays::NAllFemtoArrays];

  ClassDef(StFemtoDstReader, 0)
};

#endif // StFemtoDstReader_h
