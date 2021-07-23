// Constant
#include "/mnt/pool/rhic/1/demanov/basov/hpc_scripts/macro/Constants.h"

Int_t RunIdRange = RunIdMax.at(energy) - RunIdMin.at(energy);

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
TH1D *h_MassSqr[2][30];

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

TProfile2D *tp2_v2_CombPID[nEtaGap][2][3][7];
TProfile2D *tp2_v3_CombPID[nEtaGap][2][3][7];

TProfile2D *tp2_v2_TPC_TOF_rapidity[nEtaGap][2][3][nBinCent];
TProfile2D *tp2_v3_TPC_TOF_rapidity[nEtaGap][2][3][nBinCent];

TProfile *tp_v2_TPC_TOF_cent[nEtaGap][2][3];
TProfile *tp_v3_TPC_TOF_cent[nEtaGap][2][3];


//********mean pt *******//
TProfile2D *tp2_meanPt_CombPID[nEtaGap][2][3][7];
TProfile2D *tp2_meanPt_TPC[nEtaGap][2][3];
TProfile2D *tp2_meanPt_TOF[nEtaGap][2][3];
TProfile2D *tp2_meanPt_TPC_TOF[nEtaGap][2][3];

TH1D *h_MultPID_CombPID[nEtaGap][2][3];
TH1D *h_MultPID_TPC[nEtaGap][2][3];
TH1D *h_MultPID_TOF[nEtaGap][2][3];
TH1D *h_MultPID_TPC_TOF[nEtaGap][2][3];

//чекнуть идентификацию
TH2D *h2_dEdxVsPt_all;
TH2D *h2_dEdxVsPt_all_TPCandTOF;
TH2D *h2_dEdxVsPt_all_combPID;
TH2D *h2_m2VsPt_all;
TH2D *h2_m2VsPt_all_TPCandTOF;
TH2D *h2_m2VsPt_all_combPID;

Double_t meanNSigma_xy[2][3][30] = {0.};
Double_t sigmaNSigma_xy[2][3][30] = {0.};
Double_t meanM2_xy[2][3][30] = {0.};
Double_t sigmaM2_xy[2][3][30] = {0.};

TH2D *h2_m2vsnSigmaPion_piKp[2][20];
TH2D *h2_m2vsnSigmaPion_piKp_new[2][20];

TH2D *h2_m2vsnSigmaPion_piKp_tof[2][20];
TH2D *h2_m2vsnSigmaPion_piKp_tof_new[2][20];

void checkQvectorAndPsi(bool mode1, bool mode2){

  // check histograms
  for(Int_t l = 0; l < 2; l++) {
    //loop by cent                                                                                                  
    for(Int_t c = 0; c < nBinCent; c++) {
      //loop by eta-gap 
      for(Int_t i = 0; i < nEtaGap; i++) {
        //loop by VtxZ bins
        for(Int_t z=0; z < nBinVtxZ_PID; z++){
          if(mode1 == true || mode2 == true ){
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

  if(check_QvectorAndPsi == true){

  }

void Mode_raw(){




}

  /*////////////////////////////////////////////////////////////////////////////////////////*/
 /*___________________________DESCRIPTION OF FUNCTIONS_____________________________________*/
/*////////////////////////////////////////////////////////////////////////////////////////*/

//********************CHECK EVENT ON GOOD********************//
