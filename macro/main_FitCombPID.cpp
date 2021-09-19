#include <iostream>

#include <TROOT.h>
#include <TStopwatch.h>
#include <TFile.h>
#include <TString.h>

#include "combPID_Gaus.C"
//#include "functions.C"

int main(int argc, char** argv)
{
  TString inFileName, inFitName, outFileName, WorkMode, Energy, Pt;
  TStopwatch timer;
  if (argc < 13)
  {
    std::cerr << "./combPID_Gaus -i INPUTFILE -f INFITFILE -o OUTPUTFILE -g ENERGY -m WORKMODE -p PT" << std::endl;
    return 14;
  }

  for (int i=1;i<argc;i++)
  {
    if (std::string(argv[i]) != "-i" &&
        std::string(argv[i]) != "-f" &&
        std::string(argv[i]) != "-o" &&
        std::string(argv[i]) != "-g" &&
        std::string(argv[i]) != "-m" &&
        std::string(argv[i]) != "-p" )
    {
      std::cerr << "\nUnknown parameter: " << argv[i] << std::endl;
      return 15;
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
	return 16;
      }
      if (std::string(argv[i]) == "-f" && i != argc-1)
      {
        inFitName = TString(argv[++i]);
      }
      if (std::string(argv[i]) == "-f" && i == argc-1)
      {
        std::cerr << "\nfit file was not specified!" << std::endl;
	return 17;
      }
  if (std::string(argv[i]) == "-o" && i != argc-1)
      {
        outFileName = TString(argv[++i]);
      }
      if (std::string(argv[i]) == "-o" && i == argc-1)
      {
        std::cerr << "\nOutput file was not specified!" << std::endl;
  return 19;
      }
      if (std::string(argv[i]) == "-g" && i != argc-1)
      {
        Energy = TString(argv[++i]);
      }
      if (std::string(argv[i]) == "-g" && i == argc-1)
      {
        std::cerr << "\nenergy not specified!" << std::endl;
  return 20;
      }
      if (std::string(argv[i]) == "-m" && i != argc-1)
      {
        WorkMode = TString(argv[++i]);
      }
      if (std::string(argv[i]) == "-m" && i == argc-1)
      {
        std::cerr << "\npt bin not specified!" << std::endl;
  return 21;
      }
      if (std::string(argv[i]) == "-p" && i != argc-1)
      {
        Pt = TString(argv[++i]);
      }
      if (std::string(argv[i]) == "-p" && i == argc-1)
      {
        std::cerr << "\nenergy not specified!" << std::endl;
  return 22;
      }
    }
  }
  if (inFileName == "" || inFitName == "" || outFileName == "" || Energy == "" || WorkMode == "" || Pt == "")
  {
    std::cerr << "\nInput/Output file has not been set properly!" << std::endl;
    return 23;
  }

  combPID_Gaus(inFileName.Data(),inFitName.Data(), outFileName.Data(), Energy.Atoi(), WorkMode.Atoi(), Pt.Atoi());

  return 0;
}

