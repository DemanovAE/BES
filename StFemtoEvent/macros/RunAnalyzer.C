/**
 * \brief Helper macros for analyzing name.femtoDst.root files using FemtoDstAnalyzer.C
 *
 *  This macros takes inFileName argument with a femtoDst.root file
 *  or with a list of files (name.lis or name.list). It sets _VANILLA_ROOT_
 *  (necessary for standalone mode), loads pre-compiled libStFemtoDst.so
 *  (from StFemtoEvent), compiles and executes a text
 *  FemtoDstAnalyzer.C macro with passing inFileName to it, and
 *  cleans up the directory from the compilation products at the end.
 *
 *  Some details:
 *    inFileName - is a name of name.femtoDst.root file or a name
 *                 of a name.lis(t) files that contains a list of
 *                 name1.femtoDst.root files.
 *    NOTE: inFileName should contain either /absolutePath/inFileName
 *          or /relative2currentDir/inFileName
 *  It is assumed that FemtoDstAnalyzer.C is placed in the same
 *  directory where the RunAnalyzer.C is stored.
 *
 * \author Grigory Nigmatkulov
 * \date July 5, 2018
 */

// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"

const Char_t* defaultFileName = "/mnt/pool/rhic/2/nigmatkulov/femtoDst/dau/200gev/st_physics_17140035_raw_sum_3.femtoDst.root";

//_________________
void RunAnalyzer(const Char_t *inFileName = defaultFileName,
		 const Char_t *oFileName = "oTest.root") {
  // Next line is not needed if you are not running in a standalone mode
  gROOT->ProcessLine("#define _VANILLA_ROOT_");
  gSystem->Load("$MYSOFT/StFemtoEvent/libStFemtoDst.so");
  TString str;
  str = ".x FemtoDstAnalyzer.C+(\"";
  str += inFileName;
  str += "\",\"";
  str += oFileName;
  str += "\")";
  gROOT->ProcessLine( str.Data() );
  // Next line should be commented if you run in a batch mode
  //gROOT->ProcessLine(".!rm -f FemtoDstAnalyzer_C* ");
}
