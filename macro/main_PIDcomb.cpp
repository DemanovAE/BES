#include <iostream>

#include <TROOT.h>
#include <TStopwatch.h>
#include <TFile.h>
#include <TString.h>

#include "PIDcomb.C"
//#include "functions.C"

int main(int argc, char** argv)
{
  TString inFileName, outFileName, WorkMode, Energy;
  TStopwatch timer;
  if (argc < 9)
  {
    std::cerr << "./PIDcomb -i INPUTFILE -o OUTPUTFILE -m WORKMODE -g ENERGY" << std::endl;
    return 10;
  }

  for (int i=1;i<argc;i++)
  {
    if (std::string(argv[i]) != "-i" &&
        std::string(argv[i]) != "-o" &&
        std::string(argv[i]) != "-m" &&
        std::string(argv[i]) != "-g" )
    {
      std::cerr << "\nUnknown parameter: " << argv[i] << std::endl;
      return 11;
    }
    else
    {
      if (std::string(argv[i]) == "-i" && i != argc-1)
      {
        inFileName = TString(argv[++i]);
      }
      if (std::string(argv[i]) == "-i" && i == argc-1)
      {
        std::cerr << "\nInput file was not specified!" << std::endl;
        return 12;
      }
      if (std::string(argv[i]) == "-o" && i != argc-1)
      {
        outFileName = TString(argv[++i]);
      }
      if (std::string(argv[i]) == "-o" && i == argc-1)
      {
        std::cerr << "\nOutput file was not specified!" << std::endl;
        return 13;
      }
      if (std::string(argv[i]) == "-m" && i != argc-1)
      {
        WorkMode = TString(argv[++i]);
      }
      if (std::string(argv[i]) == "-m" && i == argc-1)
      {
        std::cerr << "\nmode not specified!" << std::endl;
  return 14;
      }
      if (std::string(argv[i]) == "-g" && i != argc-1)
      {
        Energy = TString(argv[++i]);
      }
      if (std::string(argv[i]) == "-g" && i == argc-1)
      {
        std::cerr << "\nenergy not specified!" << std::endl;
  return 15;
      }
    }
  }
  if (inFileName == "" || outFileName == "" || WorkMode == "" || Energy == "")
  {
    std::cerr << "\nInput/Output file has not been set properly!" << std::endl;
    return 16;
  }

  PIDcomb(inFileName.Data(),outFileName.Data(),WorkMode.Data(),Energy.Atoi());

  return 0;
}
