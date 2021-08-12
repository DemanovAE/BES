/**
 * \brief Example of how to read a file (list of files) using StFemtoEvent classes
 *
 * RunFemtoDstAnalyzer.C is an example of reading FemtoDst format.
 * One can use either FemtoDst file or a list of femtoDst files (inFile.lis or
 * inFile.list) as an input, and preform physics analysis
 *
 * \author Grigory Nigmatkulov, Povarov Alexey, Demanov Alexandr
 * \date May 29, 2018
 */

// This is needed for calling standalone classes
#define _VANILLA_ROOT_

// C++ headers
#include <string>
#include <vector>
#include <iostream>
#include <fstream> 
#include <map>

// ROOT headers
#include <TStopwatch.h>
#include "TProfile.h"
#include "TProfile2D.h"
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TMath.h"
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TVector.h"
#include "TVector2.h"
#include "TF2.h"
#include "TF1.h"

// FemtoDst headers

#include "/scratch2/demanov/STAR/BES/StFemtoEvent/StFemtoDstReader.h"
#include "/scratch2/demanov/STAR/BES/StFemtoEvent/StFemtoDst.h"
#include "/scratch2/demanov/STAR/BES/StFemtoEvent/StFemtoEvent.h"
#include "/scratch2/demanov/STAR/BES/StFemtoEvent/StFemtoTrack.h"
#include "/scratch2/demanov/STAR/BES/StFemtoEvent/StFemtoV0.h"
#include "/scratch2/demanov/STAR/BES/StFemtoEvent/StFemtoXi.h"

// Constant
#include "/scratch2/demanov/STAR/BES/macro/Constants.h"

// Load libraries (for ROOT_VERSTION_CODE >= 393215)
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
R__LOAD_LIBRARY(/scratch2/demanov/STAR/BES/StFemtoEvent/libStFemtoDst.so)
#endif

// inFile - is a name of name.FemtoDst.root file or a name
//          of a name.lis(t) files, that contains a list of
//          name1.FemtoDst.root, name2.FemtoDst.root, ... files

// Used function
Bool_t isGoodEvent(StFemtoEvent *const &event, const Int_t _energy, const Bool_t BadRunIdKeyFlowStage);
Bool_t isGoodTrackEP(StFemtoEvent *const &event, StFemtoTrack *const &track, const Int_t _energy);
Bool_t isGoodTrackFlow(StFemtoEvent *const &event, StFemtoTrack *const &track, const Int_t _energy);
TVector2 CalculateQvector(Int_t harmonic, StFemtoTrack *const &track);
Int_t GetEtaDirection(StFemtoTrack *const &track);
Int_t GetBinVtxZ(StFemtoEvent *const &event, const Int_t _energy);
Int_t GetBinPtRange(StFemtoTrack *const &track);
Int_t GetBinEta(StFemtoTrack *const &track);
Float_t GetRapidity(StFemtoTrack *const &track, Int_t particle);
Float_t GetMass(Int_t particle);
Double_t GetWeight(StFemtoTrack *const &track);
int PID_TPC_TOF(StFemtoTrack *const &track, const Int_t _energy);
int GetCharge(StFemtoTrack *const &track);
Bool_t TofMatchedCut(StFemtoDst *const &dst, Int_t cutTofMatched);

Double_t GetNSigmaM2(Double_t x_pt, Double_t mean, Double_t sigma);
int PID_nSigma(StFemtoTrack *const &track, const Int_t _energy, Int_t charge, TF1 *funMean[][3], TF1 *funSigma[][3], Double_t nSigma, Double_t nSigmaAntSimm);
int PID_CombPID_lowPt(StFemtoTrack *const &track, const Int_t _energy, Int_t charge, TF1 *funMean[][3], TF1 *funSigma[][3], Double_t nSigma);
int PID_CombPIDfix95(StFemtoTrack *const &track, const Int_t _energy, Double_t x, Double_t y, Int_t charge, Int_t pt, TF1 *funMean[][3], TF1 *funSigma[][3], Double_t nSigma);


//TVector2 CalcNewXY(StFemtoTrack *const &track, Int_t pt_bin_number, Int_t chPar );

TVector2 CalcNewXY(StFemtoTrack *const &track, Double_t width_nsigma, Double_t width_m2, 
                                               Double_t mean_nsigma_Kaon, Double_t mean_m2_Kaon,
                                               Double_t mean_nsigma_Pion, Double_t mean_m2_Pion);


//_________________
void FemtoDstAnalyzer_PID(const Char_t *inFile = "st_physics_12150008_raw_4030001.femtoDst.root",
                          const Char_t *outFileName = "oTest.root",
                          const Char_t *mode = "QAmode",
                          const Int_t energy = 39) {

  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

  #if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
    gSystem->Load("/scratch2/demanov/STAR/BES/StFemtoEvent/libStFemtoDst.so");
  #endif

  StFemtoDstReader* femtoReader = new StFemtoDstReader(inFile);
  femtoReader->Init();

  TStopwatch timer1;
  timer1.Start();

  // This is a way if you want to spead up IO
  std::cout << "Explicit read status for some branches" << std::endl;
  femtoReader->SetStatus("*",0);
  femtoReader->SetStatus("Event",1);
  femtoReader->SetStatus("Track",1);
  std::cout << "Status has been set" << std::endl;

  std::cout << "Now I know what to read, Master!" << std::endl;

  Int_t RunIdRange = RunIdMax.at(energy) - RunIdMin.at(energy);
  
  Int_t cent;
  Int_t RunID;

  Double_t weight;
  Double_t Phi;
  Double_t pt;
  Double_t v2;
  Double_t v3;
  Double_t sinPsi2 = 0.;
  Double_t cosPsi2 = 0.;
  Double_t sinPsi3 = 0.;
  Double_t cosPsi3 = 0.;
  Double_t dPsi2 = 0.;
  Double_t dPsi3 = 0.;

  Double_t Psi2TPC[2][nEtaGap];
  Double_t Psi3TPC[2][nEtaGap];

  Int_t mult[2][nEtaGap]; 
  TVector2 Q2vecTPC[2][nEtaGap]; 
  TVector2 Q3vecTPC[2][nEtaGap];
  TVector2 MeanQ2;
  TVector2 MeanQ3;

  Int_t halfTPC=10; //0 or 1 (10 - bad)

  // Mode works
  Bool_t mode_raw = false; 
  Bool_t mode_rec = false; 
  Bool_t mode_flow = false;

  Bool_t check_QvectorAndPsi = false;
  Bool_t check_TofMatch = false;
  Bool_t check_histo = false;
  
  if( strncmp(mode, "raw",3) == 0){
    mode_raw = true;
  }
  if( strncmp(mode, "rec",3) == 0){
    mode_rec = true;
  }
  if( strncmp(mode, "flow",4) == 0){
    mode_flow = true;
  }

  //create file 
  TFile *outFile = new TFile(outFileName, "RECREATE");
  TFile *FileRec, *FileFlow, *FileFitM2;
  outFile->cd();
  
  //histogram for Q-vectors and event planes with eta-gap and without error
  TH1D *h_Qx2[2][nEtaGap][nBinCent][nBinVtxZ_PID];
  TH1D *h_Qy2[2][nEtaGap][nBinCent][nBinVtxZ_PID];
  TH1D *h_Qx3[2][nEtaGap][nBinCent][nBinVtxZ_PID];
  TH1D *h_Qy3[2][nEtaGap][nBinCent][nBinVtxZ_PID];
  TH1D *h_Psi2[2][nEtaGap][nBinCent][nBinVtxZ_PID];
  TH1D *h_Psi3[2][nEtaGap][nBinCent][nBinVtxZ_PID];
  
  TH1D *h_PtAndP;
  TH1D *h_MultBeforeCut;
  TH1D *h_nHitsdEdx;
  TH1D *h_MultAfterCutTof;
  TH1D *h_MultAfterCutTrigger;
  TH1D *h_TofMatch;
  TH1D *h_MultWithTrigger;
  TH2D *h_BTofMatchedVsRefMult;

  TH2D *h2_nSigmaM2_all[2][3];
  TH2D *h2_nSigma[2][3][30];
  TH2D *h2_nSigma2[2][3][30];

  //***************TProfile for read in file **//
  TProfile2D *tp_read_profile;
  double binCont;
  double binEntr;

  //***** map for recentering **********************************//
  std::map<Double_t,Double_t> mp_Qx2[2][nEtaGap][nBinVtxZ_PID][nBinCent];
  std::map<Double_t,Double_t> mp_Qy2[2][nEtaGap][nBinVtxZ_PID][nBinCent];
  std::map<Double_t,Double_t> mp_Qx3[2][nEtaGap][nBinVtxZ_PID][nBinCent];
  std::map<Double_t,Double_t> mp_Qy3[2][nEtaGap][nBinVtxZ_PID][nBinCent];

  //***** map for flattening **********************************//
  std::map<Double_t,Double_t> mp_sinPsi2[2][nEtaGap][4][nBinVtxZ_PID][nBinCent];
  std::map<Double_t,Double_t> mp_cosPsi2[2][nEtaGap][4][nBinVtxZ_PID][nBinCent];
  std::map<Double_t,Double_t> mp_sinPsi3[2][nEtaGap][4][nBinVtxZ_PID][nBinCent];
  std::map<Double_t,Double_t> mp_cosPsi3[2][nEtaGap][4][nBinVtxZ_PID][nBinCent];

  //*************** FOR Q-VECTOR***************//
  TProfile2D *tp2_fill_Qx2[2][nEtaGap][nBinVtxZ_PID]; 
  TProfile2D *tp2_fill_Qx3[2][nEtaGap][nBinVtxZ_PID];
  TProfile2D *tp2_fill_Qy2[2][nEtaGap][nBinVtxZ_PID];
  TProfile2D *tp2_fill_Qy3[2][nEtaGap][nBinVtxZ_PID];
  
  //*************** FOR COS AND SIN***************//
  TProfile2D *tp2_sinPsi2East[nEtaGap][4][nBinVtxZ_PID];
  TProfile2D *tp2_cosPsi2East[nEtaGap][4][nBinVtxZ_PID];
  TProfile2D *tp2_sinPsi3East[nEtaGap][4][nBinVtxZ_PID];
  TProfile2D *tp2_cosPsi3East[nEtaGap][4][nBinVtxZ_PID];

  TProfile2D *tp2_sinPsi2West[nEtaGap][4][nBinVtxZ_PID];
  TProfile2D *tp2_cosPsi2West[nEtaGap][4][nBinVtxZ_PID];
  TProfile2D *tp2_sinPsi3West[nEtaGap][4][nBinVtxZ_PID];
  TProfile2D *tp2_cosPsi3West[nEtaGap][4][nBinVtxZ_PID];
  
  //*****Resolution^2 for Psi2 and Psi3*****// 
  TProfile *tp_SqRes2[nEtaGap];
  TProfile *tp_SqRes3[nEtaGap];

  //*****Flows*****//

  TProfile2D *tp2_v2_CombPID[nEtaGap][2][3];
  TProfile2D *tp2_v3_CombPID[nEtaGap][2][3];

  TProfile *tp2_v2_CombPID_cent[nEtaGap][2][3];
  TProfile *tp2_v3_CombPID_cent[nEtaGap][2][3];

  TProfile2D *tp2_v2_TPC_TOF[nEtaGap][2][3];
  TProfile2D *tp2_v3_TPC_TOF[nEtaGap][2][3];

  TProfile2D *tp2_v2_TPC_TOF_NoWeight[nEtaGap][2][3];
  TProfile2D *tp2_v3_TPC_TOF_NoWeight[nEtaGap][2][3];

  TProfile *tp_v2_TPC_TOF_cent[nEtaGap][2][3];
  TProfile *tp_v3_TPC_TOF_cent[nEtaGap][2][3];
  
  TProfile2D *tp2_v2_TPC_TOF_eta[nEtaGap][2][3];
  TProfile2D *tp2_v3_TPC_TOF_eta[nEtaGap][2][3];

  TProfile2D *tp2_v2_nSigmaPID[nEtaGap][2][3];
  TProfile2D *tp2_v3_nSigmaPID[nEtaGap][2][3];

  TProfile *tp2_v2_nSigmaPID_cent[nEtaGap][2][3];
  TProfile *tp2_v3_nSigmaPID_cent[nEtaGap][2][3];

  //********mean pt *******//
  TProfile2D *tp2_meanPt_CombPID[nEtaGap][2][3];
  TProfile2D *tp2_meanPt_nSigmaPID[nEtaGap][2][3];
  TProfile2D *tp2_meanPt_TPC_TOF[nEtaGap][2][3];

  //********* Comb PID ******//
  TF1 *tf1_fitFun;
  TF2 *tf2_fitFun;
  TF2 *tf2_FitNewXY[2][30];
  Double_t meanNSigma_xy[2][3][30] = {0.};
  Double_t sigmaNSigma_xy[2][3][30] = {0.};
  Double_t meanM2_xy[2][3][30] = {0.};
  Double_t sigmaM2_xy[2][3][30] = {0.};

  TF1 *tf1_FitMeanM2[2][3];
  TF1 *tf1_FitSigmaM2[2][3];

  if(check_QvectorAndPsi == true){
    // check histograms
    for(Int_t l = 0; l < 2; l++) {
      //loop by cent                                                                                                  
      for(Int_t c = 0; c < nBinCent; c++) {
        //loop by eta-gap 
        for(Int_t i = 0; i < nEtaGap; i++) {
          //loop by VtxZ bins
          for(Int_t z=0; z < nBinVtxZ_PID; z++){
            if(mode_raw == true || mode_rec == true ){
              //historam of Q-vectors for eta-gap 
              h_Qx2[l][i][c][z] = new TH1D(Form("h_Qx2%s%scent%iVtxZ%i",direction[l],NameEtaPID[i],c,z),Form("Q_{x} for #psi_{2} %s %s cent %i VtxZ %i;Q_{x}; Entries",direction[l],NameEtaPID[i],c,z),300,-1.5,1.5);
              h_Qy2[l][i][c][z] = new TH1D(Form("h_Qy2%s%scent%iVtxZ%i",direction[l],NameEtaPID[i],c,z),Form("Q_{y} for #psi_{2} %s %s cent %i VtxZ %i;Q_{y}; Entries",direction[l],NameEtaPID[i],c,z),300,-1.5,1.5);
              h_Qx3[l][i][c][z] = new TH1D(Form("h_Qx3%s%scent%iVtxZ%i",direction[l],NameEtaPID[i],c,z),Form("Q_{x} for #psi_{3} %s %s cent %i VtxZ %i;Q_{x}; Entries",direction[l],NameEtaPID[i],c,z),300,-1.5,1.5);
              h_Qy3[l][i][c][z] = new TH1D(Form("h_Qy3%s%scent%iVtxZ%i",direction[l],NameEtaPID[i],c,z),Form("Q_{y} for #psi_{3} %s %s cent %i VtxZ %i;Q_{y}; Entries",direction[l],NameEtaPID[i],c,z),300,-1.5,1.5);
            }
            //historam of event planes for east eta-gap 
            h_Psi2[l][i][c][z] = new TH1D(Form("h_Psi2%s%scent%iVtxZ%i",direction[l],NameEtaPID[i],c,z),Form("#psi_{2} %s %s cent %i;#psi_{2}; VtxZ %i",direction[l],NameEtaPID[i],c,z),100,-0.05,3.2); 
            h_Psi3[l][i][c][z] = new TH1D(Form("h_Psi3%s%scent%iVtxZ%i",direction[l],NameEtaPID[i],c,z),Form("#psi_{3} %s %s cent %i;#psi_{3}; VtxZ %i",direction[l],NameEtaPID[i],c,z),100,-0.05,2.15);
          }// loop by VtxZ
        }// loop by n  
      }// loop by cent
    }// loop by direstion 
  }
  
  if(mode_raw==true){
    
    //h_PtAndP = new TH1D("h_PtAndP","Reference multiplicity before cut ;p_T - p;Entries", 2000, -2, 2);
    h_nHitsdEdx = new TH1D("h_nHitsdEdx","nHits dE/dx ;nHits dE/dx;Entries", 100, -0.5, 99.5);
    h_MultBeforeCut = new TH1D("h_MultBeforeCut","Reference multiplicity before cut ;RefMult;Entries", 1000, -0.5, 999.5);
    h_MultAfterCutTof = new TH1D("h_MultAfterCutTof","Reference multiplicity After cut tof matched;RefMult;Entries", 1000, -0.5, 999.5);
    h_MultAfterCutTrigger = new TH1D("h_MultAfterCutTrigger","Reference multiplicity After cut trigger;RefMult;Entries", 1000, -0.5, 999.5);
    h_BTofMatchedVsRefMult = new TH2D("h_BTofMatchedVsRefMult","TOF-matched tracks vs. refMult ;refMult;TOF-matched", 600, -0.5, 599.5, 1000, -0.5, 999.5);
    h_TofMatch = new TH1D("h_TofMatch","Number of TOF-matched tracks ;bTofMatched;Entries",1000, -0.5, 999.5);
    h_MultWithTrigger = new TH1D("h_MultWithTrigger","Reference multiplicity with trigger ;RefMult;Entries",1000, -0.5, 999.5);

    // loop by direction 
    for(Int_t l = 0; l < 2; l++) {
       //loop by systematics value
      for(Int_t i = 0; i < nEtaGap; i++) {
        //loop by VtxZ
        for(Int_t z = 0; z < nBinVtxZ_PID; z++){
          //TProfile2D for recentering 
          tp2_fill_Qx2[l][i][z] = new TProfile2D(Form("tp2_Qx2%s%sVtxZ%i",direction[l],NameEtaPID[i],z),Form("<Q_{x}> for #psi_{2} %s %s VtxZ %i ;RunID;cent",direction[l],NameEtaPID[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
          tp2_fill_Qy2[l][i][z] = new TProfile2D(Form("tp2_Qy2%s%sVtxZ%i",direction[l],NameEtaPID[i],z),Form("<Q_{y}> for #psi_{2} %s %s VtxZ %i ;RunID;cent",direction[l],NameEtaPID[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
          tp2_fill_Qx3[l][i][z] = new TProfile2D(Form("tp2_Qx3%s%sVtxZ%i",direction[l],NameEtaPID[i],z),Form("<Q_{x}> for #psi_{3} %s %s VtxZ %i ;RunID;cent",direction[l],NameEtaPID[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
          tp2_fill_Qy3[l][i][z] = new TProfile2D(Form("tp2_Qy3%s%sVtxZ%i",direction[l],NameEtaPID[i],z),Form("<Q_{y}> for #psi_{3} %s %s VtxZ %i ;RunID;cent",direction[l],NameEtaPID[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
        }// loop by VtxZ
      }//loop by systematics value (Eta)
    }// loop by direction
  }

  if( mode_rec==true ) { 

  //loop by systematics value
    for(Int_t i = 0; i < nEtaGap; i++) {
      for(Int_t j = 0; j < 4; j++) {
        for(Int_t z=0; z<nBinVtxZ_PID; z++){
          tp2_sinPsi2East[i][j][z] = new TProfile2D(Form("tp2_%isinPsi2%s%sVtxZ%i",j+1,direction[0],NameEtaPID[i],z),Form("<sin(%i*2#psi_{2})> %s %s VtxZ %i ",j+1,direction[0],NameEtaPID[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
          tp2_cosPsi2East[i][j][z] = new TProfile2D(Form("tp2_%icosPsi2%s%sVtxZ%i",j+1,direction[0],NameEtaPID[i],z),Form("<cos(%i*2#psi_{2})> %s %s VtxZ %i ",j+1,direction[0],NameEtaPID[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
          tp2_sinPsi3East[i][j][z] = new TProfile2D(Form("tp2_%isinPsi3%s%sVtxZ%i",j+1,direction[0],NameEtaPID[i],z),Form("<sin(%i*3#psi_{3})> %s %s VtxZ %i ",j+1,direction[0],NameEtaPID[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
          tp2_cosPsi3East[i][j][z] = new TProfile2D(Form("tp2_%icosPsi3%s%sVtxZ%i",j+1,direction[0],NameEtaPID[i],z),Form("<cos(%i*3#psi_{3})> %s %s VtxZ %i ",j+1,direction[0],NameEtaPID[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);          
        
          tp2_sinPsi2West[i][j][z] = new TProfile2D(Form("tp2_%isinPsi2%s%sVtxZ%i",j+1,direction[1],NameEtaPID[i],z),Form("<sin(%i*2#psi_{2})> %s %s VtxZ %i ",j+1,direction[1],NameEtaPID[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
          tp2_cosPsi2West[i][j][z] = new TProfile2D(Form("tp2_%icosPsi2%s%sVtxZ%i",j+1,direction[1],NameEtaPID[i],z),Form("<cos(%i*2#psi_{2})> %s %s VtxZ %i ",j+1,direction[1],NameEtaPID[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
          tp2_sinPsi3West[i][j][z] = new TProfile2D(Form("tp2_%isinPsi3%s%sVtxZ%i",j+1,direction[1],NameEtaPID[i],z),Form("<sin(%i*3#psi_{3})> %s %s VtxZ %i ",j+1,direction[1],NameEtaPID[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
          tp2_cosPsi3West[i][j][z] = new TProfile2D(Form("tp2_%icosPsi3%s%sVtxZ%i",j+1,direction[1],NameEtaPID[i],z),Form("<cos(%i*3#psi_{3})> %s %s VtxZ %i ",j+1,direction[1],NameEtaPID[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);          
        }
      }
    }

    FileRec = new TFile(Form("%s/OUT/%iGeV/raw_%iGeV_PID.root",path,energy,energy),"READ");
    FileRec->cd();

    for(Int_t dir = 0; dir < 2; dir++) {
      for(Int_t i = 0; i < nEtaGap; i++) {
        for(Int_t z=0; z<nBinVtxZ_PID; z++){
          
          tp_read_profile = (TProfile2D*) FileRec -> Get( Form("tp2_Qx2%s%sVtxZ%i",direction[dir],NameEtaPID[i],z) );
          for(Int_t c = 0; c < nBinCent; c++){
            for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
              binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
              binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
              if (binCont == 0. || binEntr == 0.) continue;
              mp_Qx2[dir][i][z][c][(Double_t)run] = binCont;
            }
          }
          delete tp_read_profile;    

          tp_read_profile = (TProfile2D*) FileRec -> Get( Form("tp2_Qy2%s%sVtxZ%i",direction[dir],NameEtaPID[i],z) );
          for(Int_t c = 0; c < nBinCent; c++){
            for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
              binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
              binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
              if (binCont == 0. || binEntr == 0.) continue;
              mp_Qy2[dir][i][z][c][(Double_t)run] = binCont;
            }
          }
          delete tp_read_profile;  

          tp_read_profile = (TProfile2D*) FileRec -> Get( Form("tp2_Qx3%s%sVtxZ%i",direction[dir],NameEtaPID[i],z) );
          for(Int_t c = 0; c < nBinCent; c++){
            for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
              binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
              binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
              if (binCont == 0. || binEntr == 0.) continue;
              mp_Qx3[dir][i][z][c][(Double_t)run] = binCont;
            }
          }
          delete tp_read_profile;  

          tp_read_profile = (TProfile2D*) FileRec -> Get( Form("tp2_Qy3%s%sVtxZ%i",direction[dir],NameEtaPID[i],z) );
          for(Int_t c = 0; c < nBinCent; c++){
            for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
              binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
              binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
              if (binCont == 0. || binEntr == 0.) continue;
              mp_Qy3[dir][i][z][c][(Double_t)run] = binCont;
            }
          }
          delete tp_read_profile; 
        }
      }
    }

    FileRec->Close();
  }// if( mode_rec==true )

  if( mode_flow==true ) {

    for(Int_t i = 0; i < nEtaGap; i++) {
      
      tp_SqRes2[i] = new TProfile(Form("tp_SqRes2TPC%s",NameEtaPID[i]),Form("Resolution^{2} for v_{2} %s ",NameEtaPID[i]),nBinCent,0,nBinCent);
      tp_SqRes3[i] = new TProfile(Form("tp_SqRes3TPC%s",NameEtaPID[i]),Form("Resolution^{2} for v_{3} %s ",NameEtaPID[i]),nBinCent,0,nBinCent);
      
    }

    Int_t p_name = 0;

    /*
    for(Int_t par = 0; par < 3; par++){
      for(Int_t sign = 0; sign < 2; sign++){
        h2_nSigmaM2_all[sign][par] = new TH2D(Form("h1_nSigmaM2_%s_ch%i", particles[par],sign),Form("n#sigma(%s)(m^2);n#sigma(%s)", partLateX[p_name], partLateX[p_name]), 1000, 0., 5. ,2001, -10., 10. );
      }
      p_name++;
    }
    p_name=0;

    for(Int_t par = 0; par < 3; par++){
      for(Int_t sign = 0; sign < 2; sign++){
        for(Int_t pti=0; pti<(int)ptBinRange.size()-1; pti++){
          h2_nSigma[sign][par][pti] = new TH2D(Form("h2_nSigmaM2AnddEdx_%s_ch%i_pt%i", particles[par],sign,pti),Form("m^{2} vs n#sigma(%s) %.1f<p_T<%.1f GeV/c;n#sigma(%s);m^{2},(GeV/c^{2})^{2}", partLateX[p_name], ptBinRange[pti], ptBinRange[pti+1] ,partLateX[sign]), 2001, -10, 10, 2001, -10, 10 );
          h2_nSigma2[sign][par][pti] = new TH2D(Form("h2_nSigmaM2AnddEdx2_%s_ch%i_pt%i", particles[par],sign,pti),Form("m^{2} vs n#sigma(%s) %.1f<p_T<%.1f GeV/c;n#sigma(%s);m^{2},(GeV/c^{2})^{2}", partLateX[p_name], ptBinRange[pti], ptBinRange[pti+1] ,partLateX[sign]), 2001, -10, 10, 2001, -10, 10 );
        }
        p_name++;
      }
    }

    p_name = 0;
    */

    for(Int_t par = 0; par < 3; par++){
      for(Int_t sign = 0; sign < 2; sign++){

        for(Int_t i = 0; i < nEtaGap; i++) {

          /*
          tp2_meanPt_CombPID[i][sign][par] = new TProfile2D(Form("tp_meanPt%s%s%sCombPID",particles[par],particlesSign[sign],NameEtaPID[i]),Form("Mean p_{t} for bins v_{2} %s %s ; bin; p_{t} [GeV/c]",partLateX[p_name],NameEtaPID[i]),100, 0.15, 5.15,nBinCent,0,nBinCent);

          tp2_v2_CombPID[i][sign][par] = new TProfile2D(Form("tp_v2ewTPC%s%s%sCombPID",particles[par],particlesSign[sign],NameEtaPID[i]),Form("v_{2} %s of p_{t} and cent by TPC %s;p_{t} [GeV/c];cent",partLateX[p_name],NameEtaPID[i]),100, 0.15, 5.15,nBinCent,0,nBinCent);
          tp2_v3_CombPID[i][sign][par] = new TProfile2D(Form("tp_v3ewTPC%s%s%sCombPID",particles[par],particlesSign[sign],NameEtaPID[i]),Form("v_{3} %s of p_{t} and cent by TPC %s;p_{t} [GeV/c];cent",partLateX[p_name],NameEtaPID[i]),100, 0.15, 5.15,nBinCent,0,nBinCent);

          tp2_v2_CombPID_cent[i][sign][par] = new TProfile(Form("tp_v2ewTPC%s%s%sCentCombPID",particles[par],particlesSign[sign],NameEtaPID[i]),Form("v_{2} %s of p_{t} and cent by TPC %s;cent bin",partLateX[p_name],NameEtaPID[i]),nBinCent,0,nBinCent);
          tp2_v3_CombPID_cent[i][sign][par] = new TProfile(Form("tp_v3ewTPC%s%s%sCentCombPID",particles[par],particlesSign[sign],NameEtaPID[i]),Form("v_{3} %s of p_{t} and cent by TPC %s;cent bin",partLateX[p_name],NameEtaPID[i]),nBinCent,0,nBinCent);
          */
          tp2_meanPt_nSigmaPID[i][sign][par] = new TProfile2D(Form("tp_meanPt%s%s%snSigmaPID",particles[par],particlesSign[sign],NameEtaPID[i]),Form("Mean p_{t} for bins v_{2} %s %s ; bin; p_{t} [GeV/c]",partLateX[p_name],NameEtaPID[i]),100, 0.15, 5.15,nBinCent,0,nBinCent);

          tp2_v2_nSigmaPID[i][sign][par] = new TProfile2D(Form("tp_v2ewTPC%s%s%snSigmaPID",particles[par],particlesSign[sign],NameEtaPID[i]),Form("v_{2} %s of p_{t} and cent by TPC %s;p_{t} [GeV/c];cent",partLateX[p_name],NameEtaPID[i]),100, 0.15, 5.15,nBinCent,0,nBinCent);
          tp2_v3_nSigmaPID[i][sign][par] = new TProfile2D(Form("tp_v3ewTPC%s%s%snSigmaPID",particles[par],particlesSign[sign],NameEtaPID[i]),Form("v_{3} %s of p_{t} and cent by TPC %s;p_{t} [GeV/c];cent",partLateX[p_name],NameEtaPID[i]),100, 0.15, 5.15,nBinCent,0,nBinCent);

          tp2_v2_nSigmaPID_cent[i][sign][par] = new TProfile(Form("tp_v2ewTPC%s%s%sCentnSigmaPID",particles[par],particlesSign[sign],NameEtaPID[i]),Form("v_{2} %s of p_{t} and cent by TPC %s;cent bin",partLateX[p_name],NameEtaPID[i]),nBinCent,0,nBinCent);
          tp2_v3_nSigmaPID_cent[i][sign][par] = new TProfile(Form("tp_v3ewTPC%s%s%sCentnSigmaPID",particles[par],particlesSign[sign],NameEtaPID[i]),Form("v_{3} %s of p_{t} and cent by TPC %s;cent bin",partLateX[p_name],NameEtaPID[i]),nBinCent,0,nBinCent);


          tp2_meanPt_TPC_TOF[i][sign][par] = new TProfile2D(Form("tp_meanPt%s%s%sTPCandTOF",particles[par],particlesSign[sign],NameEtaPID[i]),Form("Mean p_{t} for bins v_{2} %s %s ; bin; p_{t} [GeV/c]",partLateX[p_name],NameEtaPID[i]),100, 0.15, 5.15,nBinCent,0,nBinCent);

          tp_v2_TPC_TOF_cent[i][sign][par] = new TProfile(Form("tp_v2ewTPC%s%s%sCentTPCandTOF",particles[par],particlesSign[sign],NameEtaPID[i]),Form("v_{2} %s of cent by TPC %s;cent bin",partLateX[p_name],NameEtaPID[i]),nBinCent,0,nBinCent);
          tp_v3_TPC_TOF_cent[i][sign][par] = new TProfile(Form("tp_v3ewTPC%s%s%sCentTPCandTOF",particles[par],particlesSign[sign],NameEtaPID[i]),Form("v_{3} %s of cent by TPC %s;cent bin",partLateX[p_name],NameEtaPID[i]),nBinCent,0,nBinCent);

          tp2_v2_TPC_TOF[i][sign][par] = new TProfile2D(Form("tp_v2ewTPC%s%s%sTPCandTOF",particles[par],particlesSign[sign],NameEtaPID[i]),Form("v_{2} %s of p_{t} and cent by TPC %s;p_{t} [GeV/c];cent",partLateX[p_name],NameEtaPID[i]),100, 0.15, 5.15,nBinCent,0,nBinCent);
          tp2_v3_TPC_TOF[i][sign][par] = new TProfile2D(Form("tp_v3ewTPC%s%s%sTPCandTOF",particles[par],particlesSign[sign],NameEtaPID[i]),Form("v_{3} %s of p_{t} and cent by TPC %s;p_{t} [GeV/c];cent",partLateX[p_name],NameEtaPID[i]),100, 0.15, 5.15,nBinCent,0,nBinCent);

          tp2_v2_TPC_TOF_NoWeight[i][sign][par] = new TProfile2D(Form("tp_v2ewTPC%s%s%sTPCandTOFnoWeight",particles[par],particlesSign[sign],NameEtaPID[i]),Form("v_{2} %s of p_{t} and cent by TPC %s;p_{t} [GeV/c];cent",partLateX[p_name],NameEtaPID[i]),100, 0.15, 5.15,nBinCent,0,nBinCent);
          tp2_v3_TPC_TOF_NoWeight[i][sign][par] = new TProfile2D(Form("tp_v3ewTPC%s%s%sTPCandTOFnoWeight",particles[par],particlesSign[sign],NameEtaPID[i]),Form("v_{3} %s of p_{t} and cent by TPC %s;p_{t} [GeV/c];cent",partLateX[p_name],NameEtaPID[i]),100, 0.15, 5.15,nBinCent,0,nBinCent);

        }// for(Int_t i = 0; i < n; i++){}
        p_name++;
      }// for(Int_t sign = 0; sign < 2; sign++){}
    }// for(Int_t par = 0; par < 3; par++){}

    FileRec = new TFile(Form("%s/OUT/%iGeV/raw_%iGeV_PID.root",path,energy,energy),"READ");
    FileRec -> cd();
    for(Int_t dir = 0; dir < 2; dir++) {
      for(Int_t i = 0; i < nEtaGap; i++) {
        for(Int_t z = 0; z < nBinVtxZ_PID; z++){
          
          tp_read_profile = (TProfile2D*) FileRec -> Get( Form("tp2_Qx2%s%sVtxZ%i",direction[dir],NameEtaPID[i],z) );
          
          for(Int_t c = 0; c < nBinCent; c++){
            for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
              binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
              binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
              if (binCont == 0. || binEntr == 0.) continue;
              mp_Qx2[dir][i][z][c][(Double_t)run] = binCont;
            }
          }
          delete tp_read_profile;    

          tp_read_profile = (TProfile2D*) FileRec -> Get( Form("tp2_Qy2%s%sVtxZ%i",direction[dir],NameEtaPID[i],z) );
          for(Int_t c = 0; c < nBinCent; c++){
            for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
              binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
              binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
              if (binCont == 0. || binEntr == 0.) continue;
              mp_Qy2[dir][i][z][c][(Double_t)run] = binCont;
            }
          }
          delete tp_read_profile;  

          tp_read_profile = (TProfile2D*) FileRec -> Get( Form("tp2_Qx3%s%sVtxZ%i",direction[dir],NameEtaPID[i],z) );
          for(Int_t c = 0; c < nBinCent; c++){
            for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
              binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
              binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
              if (binCont == 0. || binEntr == 0.) continue;
              mp_Qx3[dir][i][z][c][(Double_t)run] = binCont;
            }
          }
          delete tp_read_profile;  

          tp_read_profile = (TProfile2D*) FileRec -> Get( Form("tp2_Qy3%s%sVtxZ%i",direction[dir],NameEtaPID[i],z) );
          for(Int_t c = 0; c < nBinCent; c++){
            for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
              binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
              binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
              if (binCont == 0. || binEntr == 0.) continue;
              mp_Qy3[dir][i][z][c][(Double_t)run] = binCont;
            }
          }
          delete tp_read_profile;
        }
      }
    }
    FileRec->Close();

    FileFlow = new TFile(Form("%s/OUT/%iGeV/rec_%iGeV_PID.root", path,energy,energy),"READ");
    FileFlow -> cd();

    for(Int_t l = 0; l < 2; l++) {
      for(Int_t i = 0; i < nEtaGap; i++) {
        for(Int_t j = 0; j < 4; j++) {
          for(Int_t z = 0; z < nBinVtxZ_PID; z++){

            tp_read_profile =  (TProfile2D*) FileFlow -> Get( Form("tp2_%isinPsi2%s%sVtxZ%i",j+1,direction[l],NameEtaPID[i],z) );
            
            for(Int_t c = 0; c < nBinCent; c++){
              for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
                binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
                binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
                if (binCont == 0. || binEntr == 0.) continue;
                mp_sinPsi2[l][i][j][z][c][(Double_t)run] = binCont;
              }
            }
            delete tp_read_profile; 

            tp_read_profile =  (TProfile2D*) FileFlow -> Get( Form("tp2_%icosPsi2%s%sVtxZ%i",j+1,direction[l],NameEtaPID[i],z) );
            for(Int_t c = 0; c < nBinCent; c++){
              for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
                binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
                binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
                if (binCont == 0. || binEntr == 0.) continue;
                mp_cosPsi2[l][i][j][z][c][(Double_t)run] = binCont;
              }
            }
            delete tp_read_profile;

            tp_read_profile =  (TProfile2D*) FileFlow -> Get( Form("tp2_%isinPsi3%s%sVtxZ%i",j+1,direction[l],NameEtaPID[i],z) );
            for(Int_t c = 0; c < nBinCent; c++){
              for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
                binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
                binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
                if (binCont == 0. || binEntr == 0.) continue;
                mp_sinPsi3[l][i][j][z][c][(Double_t)run] = binCont;
              }
            }
            delete tp_read_profile;

            tp_read_profile =  (TProfile2D*) FileFlow -> Get( Form("tp2_%icosPsi3%s%sVtxZ%i",j+1,direction[l],NameEtaPID[i],z) );
            for(Int_t c = 0; c < nBinCent; c++){
              for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
                binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
                binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
                if (binCont == 0. || binEntr == 0.) continue;
                mp_cosPsi3[l][i][j][z][c][(Double_t)run] = binCont;
              }
            }
            delete tp_read_profile;
          }
        }
      }
    }
    FileFlow->Close();

    FileFitM2 = new TFile(Form("%s/OUT/%iGeV/FitFunM2_%iGeVRun10.root",path, energy,energy),"READ");
    FileFitM2->cd();
    /*
    for(Int_t ch=0; ch<2; ch++){
      for(Int_t pti=0; pti<(int)ptBinRange.size()-1; pti++){
        //std::cout<<(int)ptBinRange.size()<<"\t\t"<<Form("tf1_%s_m2_charge%i_pt%i",particles[par],ch, pti)<<"\n";
        tf2_fitFun = (TF2*) FileFitM2 -> Get( Form("Three2x2DGaus%ich%i",pti,ch) );
        
        for(Int_t par=0; par<3; par++){
          meanNSigma_xy[ch][par][pti] = tf2_fitFun->GetParameter(1+5*par);
          sigmaNSigma_xy[ch][par][pti] = tf2_fitFun->GetParameter(2+5*par);
          meanM2_xy[ch][par][pti] = tf2_fitFun->GetParameter(3+5*par);
          sigmaM2_xy[ch][par][pti] = tf2_fitFun->GetParameter(4+5*par);
        }
        delete tf2_fitFun;
      }
    }
    */

    for(Int_t ch=0; ch<2; ch++){
      for(Int_t par=0; par<3; par++){
        if(par==0 || par==1){
          tf1_FitMeanM2[ch][par] = new TF1(Form("f_meanM2_%s_ch%i",particles[par],ch),"[0] + [1]*x + [2]*x*x", 0.2, 3.2 );
          tf1_fitFun = (TF1*) FileFitM2 -> Get( Form("f1_meanM2_%s_ch%i",particles[par],ch) );
          tf1_FitMeanM2[ch][par]->SetParameter(0,tf1_fitFun->GetParameter(0));
          tf1_FitMeanM2[ch][par]->SetParameter(1,tf1_fitFun->GetParameter(1));
          tf1_FitMeanM2[ch][par]->SetParameter(2,tf1_fitFun->GetParameter(2));
          delete tf1_fitFun;

        }
        if(par==2){
          tf1_FitMeanM2[ch][par] = new TF1(Form("f_meanM2_%s_ch%i",particles[par],ch),"[0] + [1]*x + [2]*x*x + [3]*x*x*x", 0.2, 3.2 );
          tf1_fitFun = (TF1*) FileFitM2 -> Get( Form("f1_meanM2_%s_ch%i",particles[par],ch) );
          tf1_FitMeanM2[ch][par]->SetParameter(0,tf1_fitFun->GetParameter(0));
          tf1_FitMeanM2[ch][par]->SetParameter(1,tf1_fitFun->GetParameter(1));
          tf1_FitMeanM2[ch][par]->SetParameter(2,tf1_fitFun->GetParameter(2));
          tf1_FitMeanM2[ch][par]->SetParameter(3,tf1_fitFun->GetParameter(3));
          delete tf1_fitFun;
        }
        
        tf1_FitSigmaM2[ch][par] = new TF1(Form("f_sigmaM2_%s_ch%i",particles[par],ch),"[0] + [1]*x + [2]*x*x", 0.2, 3.2 );
        tf1_fitFun = (TF1*) FileFitM2 -> Get( Form("f1_sigmaM2_%s_ch%i",particles[par],ch) );
        tf1_FitSigmaM2[ch][par]->SetParameter(0,tf1_fitFun->GetParameter(0));
        tf1_FitSigmaM2[ch][par]->SetParameter(1,tf1_fitFun->GetParameter(1));
        tf1_FitSigmaM2[ch][par]->SetParameter(2,tf1_fitFun->GetParameter(2));
        delete tf1_fitFun;

      }
    }

    FileFitM2->Close();

  }// if( mode_flow==true )
  

  if( !femtoReader->chain() ) {
    std::cout << "No chain has been found." << std::endl;
  }
  Long64_t eventsInTree = femtoReader->tree()->GetEntries();
  std::cout << "eventsInTree: "  << eventsInTree << std::endl;
  Long64_t events2read = femtoReader->chain()->GetEntries();

  std::cout << "Number of events to read: " << events2read << std::endl;

  // Loop over events
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

    if ( iEvent % 10000 == 0) {
      std::cout << "Working on event #[" << (iEvent+1)
                << "/" << events2read << "]" << std::endl;
    }

    if (iEvent == events2read-1) {
      std::cout << "Working on event #[" << (events2read)
                << "/" << events2read << "]" << std::endl;
    }

    Bool_t readEvent = femtoReader->readFemtoEvent(iEvent);
    if( !readEvent ) {
      std::cout << "Something went wrong, Master! Nothing to analyze..." << std::endl;
      break;
    }

    // Retrieve femtoDst
    StFemtoDst *dst = femtoReader->femtoDst();

    // Retrieve event information
    StFemtoEvent *event = dst->event();
    if( !event ) {
      std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
      break;
    }
 
    if( !isGoodEvent( event,  energy, mode_flow) ) continue;
    if( !TofMatchedCut(dst, 4) ) continue;
    
    RunID = event -> runId();
    
    if( nBinCent == 9 ){
      cent = event -> cent9();
    }
    if( nBinCent == 16 ){
      cent = event -> cent16();
    }
    
    Int_t binVtxZ = GetBinVtxZ(event,energy); 
    if( binVtxZ == -1 ) continue;

    for(Int_t dir = 0; dir < 2; dir++) {
      for(Int_t s = 0; s < nEtaGap; s++) { 
        mult[dir][s] = 0;
        Q2vecTPC[dir][s].Set(0.,0.);
        Q3vecTPC[dir][s].Set(0.,0.);
      }
    }

    // Track analysis
    Int_t nTracks = dst->numberOfTracks();
    //if( nTracks >= 3000) continue;
      
    // Track loop
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

      // Retrieve i-th femto track
      StFemtoTrack *femtoTrack = dst->track(iTrk);
      
      if( !isGoodTrackEP(event, femtoTrack, energy) ) continue; 

      //h_nHitsdEdx->Fill(femtoTrack->nHitsDedx());

      Int_t EtaDir = GetEtaDirection(femtoTrack); // 0 - east; 1 - west ; -1 - bad
      Int_t EtaBin = GetBinEta(femtoTrack);
      if( EtaDir == -1 || EtaBin == -1) continue;
      
      for(Int_t i = 0; i <= EtaBin; i++ ){
        Q2vecTPC[EtaDir][i] = Q2vecTPC[EtaDir][i] + femtoTrack -> pt() * CalculateQvector(2,femtoTrack);
        Q3vecTPC[EtaDir][i] = Q3vecTPC[EtaDir][i] + femtoTrack -> pt() * CalculateQvector(3,femtoTrack);
        mult[EtaDir][i] = mult[EtaDir][i] + 1;
      }
      //h_PtAndP->Fill( femtoTrack -> pt() - femtoTrack -> p() );
    } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

    MeanQ2.Set(0.,0.);
    MeanQ3.Set(0.,0.);

    for(Int_t dir = 0; dir < 2; dir++) {
      for(Int_t i = 0; i < nEtaGap; i++) {
        if(mode_rec == true || mode_flow == true){
          MeanQ2.Set( mp_Qx2[dir][i][binVtxZ][cent][(Double_t)RunID] , mp_Qy2[dir][i][binVtxZ][cent][(Double_t)RunID]);
          MeanQ3.Set( mp_Qx3[dir][i][binVtxZ][cent][(Double_t)RunID] , mp_Qy3[dir][i][binVtxZ][cent][(Double_t)RunID]);
        }
        if( mult[dir][i] != 0 ){
          Q2vecTPC[dir][i] = Q2vecTPC[dir][i] / mult[dir][i] - MeanQ2;
          Q3vecTPC[dir][i] = Q3vecTPC[dir][i] / mult[dir][i] - MeanQ3;
        }
      }
    }

    // Flattening stage
    if(mode_flow == true){
      for(Int_t dir = 0; dir < 2; dir++) {
        for(Int_t i = 0; i < nEtaGap; i++) {
          if( Q2vecTPC[dir][i].Mod() != 0. && Q3vecTPC[dir][i].Mod() != 0.) {
            
            Psi2TPC[dir][i] = Q2vecTPC[dir][i].Phi() / 2.0;
            Psi3TPC[dir][i] = Q3vecTPC[dir][i].Phi() / 3.0;
            
            for(Int_t k = 0; k < 4; k++) {
              
              sinPsi2 = mp_sinPsi2[dir][i][k][binVtxZ][(int)cent][(Double_t)RunID];
              cosPsi2 = mp_cosPsi2[dir][i][k][binVtxZ][(int)cent][(Double_t)RunID];
              sinPsi3 = mp_sinPsi3[dir][i][k][binVtxZ][(int)cent][(Double_t)RunID];
              cosPsi3 = mp_cosPsi3[dir][i][k][binVtxZ][(int)cent][(Double_t)RunID];

              dPsi2 += -2.0*( sinPsi2 * TMath::Cos( (Double_t)(k+1)*2.0*Psi2TPC[dir][i] ) )/( 2.0*(Double_t)(k+1) ) 
                       +2.0*( cosPsi2 * TMath::Sin( (Double_t)(k+1)*2.0*Psi2TPC[dir][i] ) )/( 2.0*(Double_t)(k+1) );

              dPsi3 += -2.0*( sinPsi3 * TMath::Cos( (Double_t)(k+1)*3.0*Psi3TPC[dir][i] ) )/( 3.0*(Double_t)(k+1) ) 
                       +2.0*( cosPsi3 * TMath::Sin( (Double_t)(k+1)*3.0*Psi3TPC[dir][i] ) )/( 3.0*(Double_t)(k+1) );
            } //for(Int_t k = 0; k < 4; k++)
            Psi2TPC[dir][i] += dPsi2;
            Psi3TPC[dir][i] += dPsi3;
            dPsi2 = 0.;
            dPsi3 = 0.;
          } 
        }// for(Int_t i = 0; i < nEtaGap; i++){}
      }// for(Int_t dir = 0; dir < 3; dir++){} 
    
      for(Int_t i = 0; i < nEtaGap; i++) {
        if(Q2vecTPC[0][i].Mod() != 0. && Q3vecTPC[0][i].Mod() != 0. && Q2vecTPC[1][i].Mod() != 0. && Q3vecTPC[1][i].Mod() != 0.) {
          tp_SqRes2[i] -> Fill(cent, TMath::Cos( 2*(Psi2TPC[1][i] - Psi2TPC[0][i]) ) );
          tp_SqRes3[i] -> Fill(cent, TMath::Cos( 3*(Psi3TPC[1][i] - Psi3TPC[0][i]) ) );
        }
      }

      // Track loop
      for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

        // Retrieve i-th femto track
        StFemtoTrack *femtoTrack = dst->track(iTrk);
        
        if( !isGoodTrackFlow(event, femtoTrack, energy) ) continue; 
        
        Int_t EtaDir = GetEtaDirection(femtoTrack); // 0 - east; 1 - west ; -1 - bad
        Int_t EtaBin = GetBinEta(femtoTrack);
        Int_t charge = GetCharge(femtoTrack);
        Int_t ptBin = GetBinPtRange(femtoTrack);

        /*
        ///////Проверочные истограммы
        h2_nSigmaM2_all[charge][0]->Fill(femtoTrack -> pt(), (femtoTrack->massSqr() - tf1_FitMeanM2[charge][0]->Eval( femtoTrack->pt() )) / tf1_FitSigmaM2[charge][0]->Eval( femtoTrack->pt() )  );
        h2_nSigmaM2_all[charge][1]->Fill(femtoTrack -> pt(), (femtoTrack->massSqr() - tf1_FitMeanM2[charge][1]->Eval( femtoTrack->pt() )) / tf1_FitSigmaM2[charge][1]->Eval( femtoTrack->pt() )  );
        h2_nSigmaM2_all[charge][2]->Fill(femtoTrack -> pt(), (femtoTrack->massSqr() - tf1_FitMeanM2[charge][2]->Eval( femtoTrack->pt() )) / tf1_FitSigmaM2[charge][2]->Eval( femtoTrack->pt() )  );
        
        if(ptBin != -1){
          h2_nSigma[charge][0][ptBin]->Fill(femtoTrack->nSigmaPion(),  (femtoTrack->massSqr() - tf1_FitMeanM2[charge][0]->Eval( femtoTrack->pt() )) / tf1_FitSigmaM2[charge][0]->Eval( femtoTrack->pt() ) );
          h2_nSigma[charge][1][ptBin]->Fill(femtoTrack->nSigmaKaon(),  (femtoTrack->massSqr() - tf1_FitMeanM2[charge][1]->Eval( femtoTrack->pt() )) / tf1_FitSigmaM2[charge][1]->Eval( femtoTrack->pt() ) );
          h2_nSigma[charge][2][ptBin]->Fill(femtoTrack->nSigmaProton(),  (femtoTrack->massSqr() - tf1_FitMeanM2[charge][2]->Eval( femtoTrack->pt() )) / tf1_FitSigmaM2[charge][2]->Eval( femtoTrack->pt() ) );
        }
        */

        //h2_m2VsPt_all->Fill(femtoTrack->charge() * femtoTrack -> pt(), femtoTrack->massSqr());
        //h2_dEdxVsPt_all->Fill( femtoTrack->charge() * femtoTrack -> pt(), femtoTrack->dEdx() * 1e6 );

        //if( EtaDir == -1 || EtaBin == -1 || charge == -1 ||  ptBin>13 || ptBin == -1) continue;
        if( EtaDir == -1 || EtaBin == -1 || charge == -1  || ptBin == -1) continue;

        /*
        Double_t NewX = CalcNewXY(femtoTrack, sigmaNSigma_xy[charge][0][ptBin], sigmaM2_xy[charge][0][ptBin], 
                                                                              meanNSigma_xy[charge][1][ptBin], meanM2_xy[charge][1][ptBin], 
                                                                              meanNSigma_xy[charge][0][ptBin], meanM2_xy[charge][0][ptBin]).X();
        Double_t NewY = CalcNewXY(femtoTrack, sigmaNSigma_xy[charge][0][ptBin], sigmaM2_xy[charge][0][ptBin], 
                                                                              meanNSigma_xy[charge][1][ptBin], meanM2_xy[charge][1][ptBin], 
                                                                              meanNSigma_xy[charge][0][ptBin], meanM2_xy[charge][0][ptBin]).Y();
        */
        //Int_t parCombPID = PID_CombPID90(femtoTrack,energy, NewX, NewY, tf2_FitNewXY[charge][ptBin] ,charge, ptBin);
        //Int_t parCombPID = PID_CombPIDfix95(femtoTrack,energy, NewX, NewY ,charge, ptBin, tf1_FitMeanM2, tf1_FitSigmaM2, 2.0);
        //if( parTPC == -1 && parTPC == -1 && parTPCandTOF == -1 && parCombPID == -1) continue;
        //if(  parTPCandTOF == -1 && parCombPID == -1) continue;

        Int_t parTPCandTOF = PID_TPC_TOF(femtoTrack,energy);
        Int_t parnSigma = PID_nSigma(femtoTrack,energy, charge, tf1_FitMeanM2, tf1_FitSigmaM2, 2.0, 2.5);

        if(  parTPCandTOF == -1 && parnSigma == -1) continue;
        
        weight = GetWeight(femtoTrack);

        Phi = femtoTrack -> phi();
        pt = femtoTrack -> pt();

        v2 = 0.;
        v3 = 0.;

        for(Int_t s = 0; s <= EtaBin; s++) { 

          v2 = TMath::Cos( 2.0*(Phi - Psi2TPC[TMath::Abs(EtaDir-1)][s]) );
          v3 = TMath::Cos( 3.0*(Phi - Psi3TPC[TMath::Abs(EtaDir-1)][s]) );
          
          /*
          if( parCombPID != -1 ){
            tp2_v2_CombPID[s][charge][parCombPID] -> Fill(pt, (Double_t)cent, v2, 1.0);
            tp2_v3_CombPID[s][charge][parCombPID] -> Fill(pt, (Double_t)cent, v3, 1.0);

            tp2_v2_CombPID_cent[s][charge][parCombPID] -> Fill((Double_t)cent, v2);
            tp2_v3_CombPID_cent[s][charge][parCombPID] -> Fill((Double_t)cent, v3);

            tp2_meanPt_CombPID[s][charge][parCombPID] -> Fill(pt,(Double_t)cent,pt);
          }
          */
    

          if( parnSigma != -1 ){
            tp2_v2_nSigmaPID[s][charge][parnSigma] -> Fill(pt, (Double_t)cent, v2, 1.0);
            tp2_v3_nSigmaPID[s][charge][parnSigma] -> Fill(pt, (Double_t)cent, v3, 1.0);

            tp2_v2_nSigmaPID_cent[s][charge][parnSigma] -> Fill((Double_t)cent, v2);
            tp2_v3_nSigmaPID_cent[s][charge][parnSigma] -> Fill((Double_t)cent, v3);

            tp2_meanPt_nSigmaPID[s][charge][parnSigma] -> Fill(pt,(Double_t)cent,pt);

        //h2_nSigma2[charge][0][ptBin]->Fill(femtoTrack->nSigmaPion(),  (femtoTrack->massSqr() - tf1_FitMeanM2[charge][0]->Eval( femtoTrack->pt() )) / tf1_FitSigmaM2[charge][0]->Eval( femtoTrack->pt() ) );
        //h2_nSigma2[charge][1][ptBin]->Fill(femtoTrack->nSigmaKaon(),  (femtoTrack->massSqr() - tf1_FitMeanM2[charge][1]->Eval( femtoTrack->pt() )) / tf1_FitSigmaM2[charge][1]->Eval( femtoTrack->pt() ) );
        //h2_nSigma2[charge][2][ptBin]->Fill(femtoTrack->nSigmaProton(),  (femtoTrack->massSqr() - tf1_FitMeanM2[charge][2]->Eval( femtoTrack->pt() )) / tf1_FitSigmaM2[charge][2]->Eval( femtoTrack->pt() ) );


          }
          

          if( parTPCandTOF != -1 ){

            tp_v2_TPC_TOF_cent[s][charge][parTPCandTOF] -> Fill((Double_t)cent, v2);
            tp_v3_TPC_TOF_cent[s][charge][parTPCandTOF] -> Fill((Double_t)cent, v3);

            tp2_v2_TPC_TOF[s][charge][parTPCandTOF] -> Fill(pt, (Double_t)cent, v2, pt);
            tp2_v3_TPC_TOF[s][charge][parTPCandTOF] -> Fill(pt, (Double_t)cent, v3, pt);

            tp2_v2_TPC_TOF_NoWeight[s][charge][parTPCandTOF] -> Fill(pt, (Double_t)cent, v2, 1.0);
            tp2_v3_TPC_TOF_NoWeight[s][charge][parTPCandTOF] -> Fill(pt, (Double_t)cent, v3, 1.0);

            tp2_meanPt_TPC_TOF[s][charge][parTPCandTOF] -> Fill(pt,(Double_t)cent,pt);
          }
        }// for(Int_t eta = 0; eta < n; eta++)
      } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)
    }// if(mode_flow == true){}

    if(check_QvectorAndPsi == true){
      for(Int_t dir = 0; dir < 2; dir ++) {
        for(Int_t i = 0; i < nEtaGap; i++) {
          if( Q2vecTPC[dir][i].Mod() != 0. && Q3vecTPC[dir][i].Mod() != 0.) {
            if(mode_raw == true || mode_rec == true){
              h_Qx2[dir][i][cent][binVtxZ] -> Fill( (Double_t)Q2vecTPC[dir][i].X() );
              h_Qy2[dir][i][cent][binVtxZ] -> Fill( (Double_t)Q2vecTPC[dir][i].Y() );
              h_Qx3[dir][i][cent][binVtxZ] -> Fill( (Double_t)Q3vecTPC[dir][i].X() );
              h_Qy3[dir][i][cent][binVtxZ] -> Fill( (Double_t)Q3vecTPC[dir][i].Y() );
              
              h_Psi2[dir][i][cent][binVtxZ] -> Fill( Q2vecTPC[dir][i].Phi() / 2.0 );
              h_Psi3[dir][i][cent][binVtxZ] -> Fill( Q3vecTPC[dir][i].Phi() / 3.0 );
            }
            if(mode_flow == true){
              h_Psi2[dir][i][cent][binVtxZ] -> Fill( Psi2TPC[dir][i] );
              h_Psi3[dir][i][cent][binVtxZ] -> Fill( Psi3TPC[dir][i] );
            } 
          }
        }
      }
    }

    if(mode_raw == true ){
      //h_BTofMatchedVsRefMult->Fill(event->refMult(),number_tof);
      for(Int_t dir = 0; dir < 2; dir ++) {
        for(Int_t i = 0; i < nEtaGap; i++) {
          if( Q2vecTPC[dir][i].Mod() != 0. && Q3vecTPC[dir][i].Mod() != 0.) {
            tp2_fill_Qx2[dir][i][binVtxZ] -> Fill( (Double_t)RunID, (Double_t)cent, (Double_t)Q2vecTPC[dir][i].X() );
            tp2_fill_Qy2[dir][i][binVtxZ] -> Fill( (Double_t)RunID, (Double_t)cent, (Double_t)Q2vecTPC[dir][i].Y() );
            tp2_fill_Qx3[dir][i][binVtxZ] -> Fill( (Double_t)RunID, (Double_t)cent, (Double_t)Q3vecTPC[dir][i].X() );
            tp2_fill_Qy3[dir][i][binVtxZ] -> Fill( (Double_t)RunID, (Double_t)cent, (Double_t)Q3vecTPC[dir][i].Y() );
          }
        }// for(Int_t i = 0; i < nEtaGap; i++){}  
      }// for(Int_t dir = 0; dir < 2; dir++){}
    }// if(mode_raw==true){}

    if(mode_rec == true ){
      for(Int_t i = 0; i < nEtaGap; i++) {
        for(Int_t j = 0; j < 4; j++) {
          if( Q2vecTPC[0][i].Mod() != 0. && Q3vecTPC[0][i].Mod() != 0.) {
            tp2_sinPsi2East[i][j][binVtxZ] -> Fill((Double_t)RunID, (Double_t)cent, TMath::Sin( (Double_t)(j+1)*Q2vecTPC[0][i].Phi() ) );
            tp2_cosPsi2East[i][j][binVtxZ] -> Fill((Double_t)RunID, (Double_t)cent, TMath::Cos( (Double_t)(j+1)*Q2vecTPC[0][i].Phi() ) );
            tp2_sinPsi3East[i][j][binVtxZ] -> Fill((Double_t)RunID, (Double_t)cent, TMath::Sin( (Double_t)(j+1)*Q3vecTPC[0][i].Phi() ) ); 
            tp2_cosPsi3East[i][j][binVtxZ] -> Fill((Double_t)RunID, (Double_t)cent, TMath::Cos( (Double_t)(j+1)*Q3vecTPC[0][i].Phi() ) );
          }
          if( Q2vecTPC[1][i].Mod() != 0. && Q3vecTPC[1][i].Mod() != 0.) {
            tp2_sinPsi2West[i][j][binVtxZ] -> Fill((Double_t)RunID, (Double_t)cent, TMath::Sin( (Double_t)(j+1)*Q2vecTPC[1][i].Phi() ) );
            tp2_cosPsi2West[i][j][binVtxZ] -> Fill((Double_t)RunID, (Double_t)cent, TMath::Cos( (Double_t)(j+1)*Q2vecTPC[1][i].Phi() ) );
            tp2_sinPsi3West[i][j][binVtxZ] -> Fill((Double_t)RunID, (Double_t)cent, TMath::Sin( (Double_t)(j+1)*Q3vecTPC[1][i].Phi() ) ); 
            tp2_cosPsi3West[i][j][binVtxZ] -> Fill((Double_t)RunID, (Double_t)cent, TMath::Cos( (Double_t)(j+1)*Q3vecTPC[1][i].Phi() ) );
          }
        }//for(Int_t j =0; j < 4; j++){}
      }// for(Int_t i = 0; i < nEtaGap; i++){}  
    }// if(mode_raw==true){}



  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

  outFile->Write();
  outFile->Close();

  femtoReader->Finish();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
	          << std::endl;

  timer1.Stop();
  timer1.Print();
}

  /*////////////////////////////////////////////////////////////////////////////////////////*/
 /*___________________________DESCRIPTION OF FUNCTIONS_____________________________________*/
/*////////////////////////////////////////////////////////////////////////////////////////*/

//********************CHECK EVENT ON GOOD********************//
Bool_t isGoodEvent(StFemtoEvent *const &event, const Int_t _energy, const Bool_t BadRunIdKeyFlowStage) {
  
  TVector3 pVtx = event->primaryVertex();
  // Reject vertices that are far from the central membrane along the beam
  if( TMath::Abs( pVtx.Z() ) > CutVtxZ.at(_energy) ) return false;
  //if( sqrt( pow( pVtx.X(), 2) + pow(pVtx.Y() + 0.8847, 2) ) > 1 ) check = false; 14.5
  if( sqrt( pow( pVtx.X(), 2) + pow(pVtx.Y() + CutDeltaVtxY.at(_energy), 2) ) > CutVtxR.at(_energy) ) return false;
  
  std::vector<Int_t> Bad = BadRuns.at(_energy);
  if(BadRunIdKeyFlowStage==true){
    if( std::find(Bad.begin(), Bad.end(), event -> runId()) != Bad.end() ) return false;
  }
  // if( event->numberOfTofMatched() <= 4) return false; //в данных 2010 года numberOfTofMatched не посчитан. Считать вручную
    /*  
    std::vector<Int_t> trig = Trigger.at(_energy);
      for(Int_t i = 0; i < trig.size(); i++){
        if( event->isTrigger(trig[i]) ) return true;
      }
      */
  return true;
}// isGoodEvent(){}

//********************CHECK TRACK ON GOOD********************//
Bool_t isGoodTrackEP(StFemtoEvent *const &event, StFemtoTrack *const &track, const Int_t _energy) {

  TVector3 pVtx = event->primaryVertex();

  if ( !track ) return false;
  // Must be a primary track
  if ( !track->isPrimary() ) return false;
  if ( ( track -> dEdx() ) == 0. ) return false;
  // Simple single-track cut
  if( track -> gMom().Mag() < 0.1) return false; 
  if( track -> gDCA(pVtx).Mag() > CutDCApidEP.at(_energy) ) return false;    
  if( TMath::Abs( track -> eta() ) > 1.0 ) return false; 
  if( track -> nHits() < CutnHits) return false;
  if( track -> pt() < CutPtotEPMin_PID) return false;
  if( track -> pt() > CutPtotEPMax_PID) return false; 
  if( ( (Double_t)track -> nHits() )/( (Double_t)track -> nHitsPoss() )  < CutnHitsRatio ) return false;

  return true; 
}// isGoodTrack(){}

//********************CHECK TRACK FLOW ON GOOD********************//
Bool_t isGoodTrackFlow(StFemtoEvent *const &event, StFemtoTrack *const &track, const Int_t _energy) {

  TVector3 pVtx = event->primaryVertex();

  if ( !track ) return false;
  // Must be a primary track
  if ( !track->isPrimary() ) return false;
  if ( ( track -> dEdx() ) == 0. ) return false;
  // Simple single-track cut
  if( track -> gMom().Mag() < 0.1) return false; 
  if( track -> gDCA(pVtx).Mag() > CutDCApidFlow.at(_energy) ) return false;
  if( TMath::Abs( track -> eta() ) > 1.0 ) return false; 
  if( track -> nHits() < CutnHits) return false;
  if( track -> pt() < CutPtotFlowMin_PID) return false;
  if( track -> pt() > CutPtotFlowMax_PID) return false; 
  if( ( (Double_t)track -> nHits() )/( (Double_t)track -> nHitsPoss() )  < CutnHitsRatio ) return false;

  return true; 
}// isGoodTrack(){}

//***********************CALCULATE Q-VECTOR**********************//
TVector2 CalculateQvector(Int_t harmonic, StFemtoTrack *const &track){
  
  TVector2 qv(0.,0.);
  qv.Set( TMath::Cos( harmonic * track -> phi() ) , TMath::Sin( harmonic * track -> phi() ) );
  return qv;

}// CalculateQvector(){}

Int_t GetEtaDirection(StFemtoTrack *const &track){

  if(track -> eta() < 0.) return 0;
  if(track -> eta() > 0.) return 1;
  return -1;
}

Int_t GetBinVtxZ(StFemtoEvent *const &event, const Int_t _energy){
  
  TVector3 pVtx = event->primaryVertex();

  if(pVtx.Z() == (-1.0 * CutVtxZ.at(_energy))) return 0;
  if(pVtx.Z() == CutVtxZ.at(_energy)) return (nBinVtxZ_PID-1);

  Int_t bin = -1;

  bin = (Int_t)( TMath::Abs(CutVtxZ.at(_energy) + pVtx.Z()) / (2.0*CutVtxZ.at(_energy) / nBinVtxZ_PID));
  
  //std::cout << bin <<"\t\t" << pVtx.Z()<<std::endl;

  if(bin > nBinVtxZ_PID)return -1;

  return bin;

}

Int_t GetBinEta(StFemtoTrack *const &track){

  for(Int_t i = nEtaGap-1 ; i >= 0; i--){
    if( TMath::Abs(track -> eta()) > EtaVecPID[i]){ 
      return i;
    } 
  }
  return -1;
}

Float_t GetRapidity(StFemtoTrack *const &track, Int_t particle){

  Float_t E = sqrt( pow( track->p(), 2) + pow( GetMass(particle) ,2) );
  return 0.5 * log( (E + track->pMom().Z() ) / ( E - track->pMom().Z() ) );
}

Float_t GetMass(Int_t particle){
  if(particle==0){
    return pion_mass;
  }
  if(particle==1){
    return kaon_mass;
  }
  if(particle==2){
    return proton_mass;
  }
  return 100.;
}

Double_t GetWeight(StFemtoTrack *const &track){
  
  Double_t w;
  if (track->pt() < 2.0){
    w = track->pt();
  }
  else{
    w = 2.0;
  }
  return w;
}


int PID_TPC_TOF(StFemtoTrack *const &track, const Int_t _energy){

  if ( track->isTofTrack() ){
    if( TMath::Abs( track->nSigmaPion() ) < 2.0 && track->massSqr() > SqMdown[0] && track->massSqr() < SqMup[0]){
      return 0;
    }
    if( TMath::Abs( track->nSigmaKaon() ) < 2.0 && track->massSqr()> SqMdown[1] && track->massSqr() < SqMup[1]){
      return 1;
    }
    if( TMath::Abs( track->nSigmaProton() ) < 2.0 && track->massSqr()> SqMdown[2] && track->massSqr() < SqMup[2]){
      return 2;
    }
  }
  return -1;
}// PID

int PID_CombPIDfix95(StFemtoTrack *const &track, const Int_t _energy, Double_t x, Double_t y, Int_t charge, Int_t pt, TF1 *funMean[][3], TF1 *funSigma[][3], Double_t nSigma){

  if ( track->isTofTrack() ){
    if(track->pt() < 1.2){
      return PID_CombPID_lowPt(track, _energy, charge, funMean, funSigma, nSigma);
    }
    if( CutPion27GeV_95[charge][pt][0] < x && CutPion27GeV_95[charge][pt][1] > x && CutPion27GeV_95[charge][pt][2] < y && CutPion27GeV_95[charge][pt][3] > y ){
      return 0;
    }

    if( CutKaon27GeV_95[charge][pt][0] < x && CutKaon27GeV_95[charge][pt][1] > x && CutKaon27GeV_95[charge][pt][2] < y && CutKaon27GeV_95[charge][pt][3] > y){
      return 1;
    }
  }
  return -1;
}// PID

int PID_nSigma(StFemtoTrack *const &track, const Int_t _energy, Int_t charge, TF1 *funMean[][3], TF1 *funSigma[][3], Double_t nSigma, Double_t nSigmaAntSimm){

  if(sqrt( pow( GetNSigmaM2(track->massSqr(), funMean[charge][0]->Eval( track->pt() ), funSigma[charge][0]->Eval( track->pt())), 2) + pow( 2.0 * track->nSigmaPion(), 2)) < nSigma &&
     sqrt( pow( GetNSigmaM2(track->massSqr(), funMean[charge][1]->Eval( track->pt() ), funSigma[charge][1]->Eval( track->pt())), 2) + pow( 2.0 * track->nSigmaKaon(), 2)) > nSigmaAntSimm &&
     sqrt( pow( GetNSigmaM2(track->massSqr(), funMean[charge][2]->Eval( track->pt() ), funSigma[charge][2]->Eval( track->pt())), 2) + pow( 2.0 * track->nSigmaProton(), 2)) > nSigmaAntSimm){
    return 0;
  }

  if(sqrt( pow( GetNSigmaM2(track->massSqr(), funMean[charge][1]->Eval( track->pt() ), funSigma[charge][1]->Eval( track->pt())), 2) + pow( 2.0 * track->nSigmaKaon(), 2)) < nSigma &&
     sqrt( pow( GetNSigmaM2(track->massSqr(), funMean[charge][0]->Eval( track->pt() ), funSigma[charge][0]->Eval( track->pt())), 2) + pow( 2.0 * track->nSigmaPion(), 2)) > nSigmaAntSimm &&
     sqrt( pow( GetNSigmaM2(track->massSqr(), funMean[charge][2]->Eval( track->pt() ), funSigma[charge][2]->Eval( track->pt())), 2) + pow( 2.0 * track->nSigmaProton(), 2)) > nSigmaAntSimm){
    return 1;
  }

  if(sqrt( pow( GetNSigmaM2(track->massSqr(), funMean[charge][2]->Eval( track->pt() ), funSigma[charge][2]->Eval( track->pt())), 2) + pow( 2.0 * track->nSigmaProton(), 2)) < nSigma &&
     sqrt( pow( GetNSigmaM2(track->massSqr(), funMean[charge][0]->Eval( track->pt() ), funSigma[charge][0]->Eval( track->pt())), 2) + pow( 2.0 * track->nSigmaKaon(), 2)) > nSigmaAntSimm &&
     sqrt( pow( GetNSigmaM2(track->massSqr(), funMean[charge][1]->Eval( track->pt() ), funSigma[charge][1]->Eval( track->pt())), 2) + pow( 2.0 * track->nSigmaPion(), 2)) > nSigmaAntSimm){
    return 2;
  }

  return -1;
}// PID

Double_t GetNSigmaM2(Double_t x_pt, Double_t mean, Double_t sigma){
  return (x_pt - mean)/sigma;
}

int PID_CombPID_lowPt(StFemtoTrack *const &track, const Int_t _energy, Int_t charge, TF1 *funMean[][3], TF1 *funSigma[][3], Double_t nSigma){

  if( TMath::Abs( (track->massSqr() - funMean[charge][0]->Eval( track->pt() )) / funSigma[charge][0]->Eval( track->pt() ) ) < nSigma &&
    TMath::Abs( (track->massSqr() - funMean[charge][1]->Eval( track->pt() )) / funSigma[charge][1]->Eval( track->pt() ) ) > 3.0 &&
    TMath::Abs( (track->massSqr() - funMean[charge][2]->Eval( track->pt() )) / funSigma[charge][2]->Eval( track->pt() ) ) > 3.0){
    return 0;
  }

  if( TMath::Abs( (track->massSqr() - funMean[charge][1]->Eval( track->pt() )) / funSigma[charge][1]->Eval( track->pt() ) ) < nSigma &&
    TMath::Abs( (track->massSqr() - funMean[charge][0]->Eval( track->pt() )) / funSigma[charge][0]->Eval( track->pt() ) ) > 3.0 &&
    TMath::Abs( (track->massSqr() - funMean[charge][2]->Eval( track->pt() )) / funSigma[charge][2]->Eval( track->pt() ) ) > 3.0){
    return 1;
  }

  if( TMath::Abs( (track->massSqr() - funMean[charge][2]->Eval( track->pt() )) / funSigma[charge][2]->Eval( track->pt() ) ) < nSigma &&
    TMath::Abs( (track->massSqr() - funMean[charge][1]->Eval( track->pt() )) / funSigma[charge][1]->Eval( track->pt() ) ) > 3.0 &&
    TMath::Abs( (track->massSqr() - funMean[charge][0]->Eval( track->pt() )) / funSigma[charge][0]->Eval( track->pt() ) ) > 3.0){
    return 2;
  }

  return -1;
}// PID

int GetCharge(StFemtoTrack *const &track){

  if((Int_t)track -> charge() > 0) return 0;
  if((Int_t)track -> charge() < 0) return 1;

  return -1;
}// PID

Int_t GetBinPtRange(StFemtoTrack *const &track){

  if ( track->isTofTrack() ){
    for(Int_t i = 0; i < (int)ptBinRange.size()-1; i++){
      if( track->pt() >= ptBinRange[i] && track->pt() < ptBinRange[i+1]){
        return i;
      }
    }
  }
  return -1;
}

Bool_t TofMatchedCut(StFemtoDst *const &dst, Int_t cutTofMatched){

    Int_t nTrack = dst->numberOfTracks();
    Int_t number_tof=0;

    for(Int_t iTrk=0; iTrk<nTrack; iTrk++) {
      StFemtoTrack *femtoTrack = dst->track(iTrk);
      if ( !femtoTrack ) continue;
      if ( femtoTrack->isTofTrack()){
        number_tof++;
      }
      if(number_tof > cutTofMatched) return true;
    }

    return false;

}

TVector2 CalcNewXY(StFemtoTrack *const &track, Double_t width_nsigma, Double_t width_m2, 
                                               Double_t mean_nsigma_Kaon, Double_t mean_m2_Kaon,
                                               Double_t mean_nsigma_Pion, Double_t mean_m2_Pion){

  TVector2 qv(0.,0.);

  Double_t f_scale = width_nsigma / width_m2;

  Double_t den = (mean_nsigma_Kaon - mean_nsigma_Pion)/f_scale;
  
  Double_t alpha = -1*TMath::ATan2(mean_m2_Kaon - mean_m2_Pion , den); 
  
  Double_t x1 = ( track->nSigmaPion() - mean_nsigma_Pion ) / f_scale;
  Double_t y1 = track->massSqr() - mean_m2_Pion;

  qv.Set(TMath::Cos(alpha)*x1 - TMath::Sin(alpha)*y1 , TMath::Sin(alpha)*x1 + TMath::Cos(alpha)*y1 );

  return qv;
}

