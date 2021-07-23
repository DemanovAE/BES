//
// Allows to read femtoDst file or a list of femtoDst files
//

// C++ headers
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <assert.h>

// FemtoDst headers
#include "StFemtoMessMgr.h"
#include "StFemtoDstReader.h"
#include "StFemtoEvent.h"
#include "StFemtoTrack.h"
#include "StFemtoV0.h"
#include "StFemtoXi.h"
#include "StFemtoEpdHit.h"
#include "StFemtoArrays.h"
#include "StFemtoDst.h"

// ROOT headers
#include "TRegexp.h"

ClassImp(StFemtoDstReader)

//_________________
StFemtoDstReader::StFemtoDstReader(const Char_t* inFileName) :
  mFemtoDst(new StFemtoDst()), mChain(NULL), mTree(NULL),
  mEventCounter(0), mFemtoArrays{}, mStatusArrays{} {
  // Constructor
  streamerOff();
  createArrays();
  std::fill_n(mStatusArrays, sizeof(mStatusArrays) / sizeof(mStatusArrays[0]), 1);
  mInputFileName = inFileName;
}

//_________________
StFemtoDstReader::~StFemtoDstReader() {
  // Destructor
  if(mChain) {
    delete mChain;
  }
  if(mFemtoDst) {
    delete mFemtoDst;
  }
}

//_________________
void StFemtoDstReader::clearArrays() {
  // Clear arrays
  for(Int_t iArr=0; iArr<StFemtoArrays::NAllFemtoArrays; iArr++) {
    mFemtoArrays[iArr]->Clear();
  }
}

//_________________
void StFemtoDstReader::SetStatus(const Char_t *branchNameRegex, Int_t enable) {
  // Set branch status
  if(strncmp(branchNameRegex, "St", 2) == 0) {
    // Ignore first "St"
    branchNameRegex += 2;
  }

  TRegexp re(branchNameRegex, 1);
  for(Int_t iArr=0; iArr<StFemtoArrays::NAllFemtoArrays; iArr++) {
    Ssiz_t len;
    if(re.Index(StFemtoArrays::femtoArrayNames[iArr], &len) < 0) continue;
    LOG_INFO << "StFemtoDstMaker::SetStatus " << enable
	     << " to " << StFemtoArrays::femtoArrayNames[iArr] << endm;
    mStatusArrays[iArr] = enable;
  }

  setBranchAddresses(mChain);
}

//_________________
void StFemtoDstReader::setBranchAddresses(TChain *chain) {
  // Set addresses of branches listed in femtoArrays
  if (!chain) return;
  chain->SetBranchStatus("*", 0);
  TString ts;
  for (Int_t i = 0; i < StFemtoArrays::NAllFemtoArrays; ++i) {
    if (mStatusArrays[i] == 0) continue;
    char const* bname = StFemtoArrays::femtoArrayNames[i];
    TBranch* tb = chain->GetBranch(bname);
    if (!tb) {
      LOG_WARN << "setBranchAddress: Branch name " << bname << " does not exist!" << endm;
      continue;
    }
    ts = bname;
    ts += "*";
    chain->SetBranchStatus(ts, 1);
    chain->SetBranchAddress(bname, mFemtoArrays + i);
    assert(tb->GetAddress() == (char*)(mFemtoArrays + i));
  }
  mTree = mChain->GetTree();
}

//_________________
void StFemtoDstReader::streamerOff() {
  // This is to to save space on the file. No need for TObject bits for this structure.
  // see: https://root.cern.ch/doc/master/classTClass.html#a606b0442d6fec4b1cd52f43bca73aa51
  StFemtoEvent::Class()->IgnoreTObjectStreamer();
  StFemtoTrack::Class()->IgnoreTObjectStreamer();
  StFemtoEpdHit::Class()->IgnoreTObjectStreamer();
}

//_________________
void StFemtoDstReader::createArrays() {
  // Create arrays
  for(Int_t iArr=0; iArr<StFemtoArrays::NAllFemtoArrays; iArr++) {
    mFemtoArrays[iArr] = new TClonesArray(StFemtoArrays::femtoArrayTypes[iArr],
					  StFemtoArrays::femtoArraySizes[iArr]);
  }
  mFemtoDst->set(mFemtoArrays);
}

//_________________
void StFemtoDstReader::Finish() {
  // Finish and clean everything
  if(mChain) {
    delete mChain;
  }
  mChain = NULL;
}

//_________________
void StFemtoDstReader::Init() {
  // FemtoDst initialization
  if(!mChain) {
    mChain = new TChain("FemtoDst");
  }

  std::string const dirFile = mInputFileName.Data();

  if( dirFile.find(".list") != std::string::npos ||
      dirFile.find(".lis") != std::string::npos ) {

    std::ifstream inputStream( dirFile.c_str() );

    if(!inputStream) {
      LOG_ERROR << "ERROR: Cannot open list file " << dirFile << endm;
    }

    Int_t nFile = 0;
    std::string file;
    while(getline(inputStream, file)) {
      if(file.find(".femtoDst.root") != std::string::npos) {
        TFile* ftmp = TFile::Open(file.c_str());
        if(ftmp && !ftmp->IsZombie() && ftmp->GetNkeys()) {
          LOG_INFO << " Read in femtoDst file " << file << endm;
          mChain->Add(file.c_str());
          ++nFile;
        } //if(ftmp && !ftmp->IsZombie() && ftmp->GetNkeys())

        if (ftmp) {
	        ftmp->Close();
	      } // if (ftmp)
      } //if(file.find(".femtoDst.root") != std::string::npos)
    } //while (getline(inputStream, file))

    LOG_INFO << " Total " << nFile << " files have been read in. " << endm;
  } //if(dirFile.find(".list") != std::string::npos || dirFile.find(".lis" != string::npos))
  else if(dirFile.find(".femtoDst.root") != std::string::npos) {
    mChain->Add( dirFile.c_str() );
  }
  else {
    LOG_WARN << " No good input file to read ... " << endm;
  }

  if(mChain) {
    setBranchAddresses(mChain);
    mChain->SetCacheSize(50e6);
    mChain->AddBranchToCache("*");
    mFemtoDst->set(mFemtoArrays);
  }
}

//_________________
Bool_t StFemtoDstReader::readFemtoEvent(Long64_t iEvent) {
  // Read femtoEvent
  Int_t mStatusRead = true; // true - okay, false - nothing to read

  if (!mChain) {
    LOG_WARN << " No input files ... ! EXIT" << endm;
    mStatusRead = false;
    return mStatusRead;
  }

  Int_t bytes = mChain->GetEntry(mEventCounter++);
  Int_t nCycles = 0;
  while( bytes <= 0) {
    if( mEventCounter >= mChain->GetEntriesFast() ) {
    }

    LOG_WARN << "Encountered invalid entry or I/O error while reading event "
	     << mEventCounter << " from \"" << mChain->GetName() << "\" input tree\n";
    bytes = mChain->GetEntry(mEventCounter++);
    nCycles++;
    LOG_WARN << "Not input has been found for: " << nCycles << " times" << endm;
    if(nCycles >= 10) {
      LOG_ERROR << "Terminating StFemtoDstReader::ReadProcess(Long64_t) after "
		<< nCycles << " times!" << endm;
      mStatusRead = false;
      break;
    }
  }
  return mStatusRead;
}
