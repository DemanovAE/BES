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

// FemtoDst headers

#include "StFemtoDstReader.h"
#include "StFemtoDst.h"
#include "StFemtoEvent.h"
#include "StFemtoTrack.h"
#include "StFemtoV0.h"
#include "StFemtoXi.h"

// Constant
#include "Constants.h"
#include "functions.C"

// Load libraries (for ROOT_VERSTION_CODE >= 393215)
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
R__LOAD_LIBRARY(libStFemtoDst.so)
#endif

// inFile - is a name of name.FemtoDst.root file or a name
//          of a name.lis(t) files, that contains a list of
//          name1.FemtoDst.root, name2.FemtoDst.root, ... files

// Used function
Bool_t isGoodEvent(StFemtoEvent *const &event, const Int_t _energy, const Bool_t BadRunIdKeyFlowStage);
Bool_t isGoodTrackEP(StFemtoEvent *const &event, StFemtoTrack *const &track, const Int_t _energy);
Bool_t isGoodTrackFlow(StFemtoEvent *const &event, StFemtoTrack *const &track, const Int_t _energy);

//_________________
void FemtoDstAnalyzer_Hadrons(const Char_t *inFile = "st_physics_12150008_raw_4030001.femtoDst.root",
                          const Char_t *outFileName = "oTest.root",
                          const Char_t *mode = "QAmode",
                          const Int_t energy = 39) {

  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

  #if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
    gSystem->Load("libStFemtoDst.so");
  #endif

  StFemtoDstReader* femtoReader = new StFemtoDstReader(inFile);
  femtoReader->Init();

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

  Double_t Psi2TPC[2][nEtaGapHadrons];
  Double_t Psi3TPC[2][nEtaGapHadrons];

  Int_t mult[2][nEtaGapHadrons]; 
  TVector2 Q2vecTPC[2][nEtaGapHadrons]; 
  TVector2 Q3vecTPC[2][nEtaGapHadrons];
  TVector2 MeanQ2;
  TVector2 MeanQ3;

  // Mode works
  Bool_t mode_raw = false; 
  Bool_t mode_rec = false; 
  Bool_t mode_flow = false;
  Bool_t Check_histo = true; 
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
  TFile *FileRec, *FileFlow;
  outFile->cd();
  
  //histogram for Q-vectors and event planes with eta-gap and without error
  TH1D *h_Qx2[2][nEtaGapHadrons][nBinCent][nBinVtxZ_Hadrons];
  TH1D *h_Qy2[2][nEtaGapHadrons][nBinCent][nBinVtxZ_Hadrons];
  TH1D *h_Qx3[2][nEtaGapHadrons][nBinCent][nBinVtxZ_Hadrons];
  TH1D *h_Qy3[2][nEtaGapHadrons][nBinCent][nBinVtxZ_Hadrons];
  TH1D *h_Psi2[2][nEtaGapHadrons][nBinCent][nBinVtxZ_Hadrons];
  TH1D *h_Psi3[2][nEtaGapHadrons][nBinCent][nBinVtxZ_Hadrons];

  //***************TProfile for read in file **//
  TProfile2D *tp_read_profile;
  double binCont;
  double binEntr;

  //***** map for recentering **********************************//
  std::map<Double_t,Double_t> mp_Qx2[2][nEtaGapHadrons][nBinVtxZ_Hadrons][nBinCent];
  std::map<Double_t,Double_t> mp_Qy2[2][nEtaGapHadrons][nBinVtxZ_Hadrons][nBinCent];
  std::map<Double_t,Double_t> mp_Qx3[2][nEtaGapHadrons][nBinVtxZ_Hadrons][nBinCent];
  std::map<Double_t,Double_t> mp_Qy3[2][nEtaGapHadrons][nBinVtxZ_Hadrons][nBinCent];

  //***** map for flattening **********************************//
  std::map<Double_t,Double_t> mp_sinPsi2[2][nEtaGapHadrons][4][nBinVtxZ_Hadrons][nBinCent];
  std::map<Double_t,Double_t> mp_cosPsi2[2][nEtaGapHadrons][4][nBinVtxZ_Hadrons][nBinCent];
  std::map<Double_t,Double_t> mp_sinPsi3[2][nEtaGapHadrons][4][nBinVtxZ_Hadrons][nBinCent];
  std::map<Double_t,Double_t> mp_cosPsi3[2][nEtaGapHadrons][4][nBinVtxZ_Hadrons][nBinCent];

  //*************** FOR Q-VECTOR***************//
  TProfile2D *tp2_fill_Qx2[2][nEtaGapHadrons][nBinVtxZ_Hadrons]; 
  TProfile2D *tp2_fill_Qx3[2][nEtaGapHadrons][nBinVtxZ_Hadrons];
  TProfile2D *tp2_fill_Qy2[2][nEtaGapHadrons][nBinVtxZ_Hadrons];
  TProfile2D *tp2_fill_Qy3[2][nEtaGapHadrons][nBinVtxZ_Hadrons];
  
  //*************** FOR COS AND SIN***************//
  TProfile2D *tp2_sinPsi2East[nEtaGapHadrons][4][nBinVtxZ_Hadrons];
  TProfile2D *tp2_cosPsi2East[nEtaGapHadrons][4][nBinVtxZ_Hadrons];
  TProfile2D *tp2_sinPsi3East[nEtaGapHadrons][4][nBinVtxZ_Hadrons];
  TProfile2D *tp2_cosPsi3East[nEtaGapHadrons][4][nBinVtxZ_Hadrons];

  TProfile2D *tp2_sinPsi2West[nEtaGapHadrons][4][nBinVtxZ_Hadrons];
  TProfile2D *tp2_cosPsi2West[nEtaGapHadrons][4][nBinVtxZ_Hadrons];
  TProfile2D *tp2_sinPsi3West[nEtaGapHadrons][4][nBinVtxZ_Hadrons];
  TProfile2D *tp2_cosPsi3West[nEtaGapHadrons][4][nBinVtxZ_Hadrons];
  
  //*****Resolution^2 for Psi2 and Psi3*****// 
  TProfile *tp_SqRes2[nEtaGapHadrons];
  TProfile *tp_SqRes3[nEtaGapHadrons];

  //*****Flows*****//
  TProfile2D *tp2_v2hadrTPC[nEtaGapHadrons];
  TProfile2D *tp2_v3hadrTPC[nEtaGapHadrons];

  TProfile *tp_v2hadrTPC_cent[nEtaGapHadrons];
  TProfile *tp_v3hadrTPC_cent[nEtaGapHadrons];
  
  TProfile2D *tp2_v2hadrTPC_eta[nEtaGapHadrons][nBinCent];
  TProfile2D *tp2_v3hadrTPC_eta[nEtaGapHadrons][nBinCent];

  TH1D *hEta[nBinCent];
  TH1D *hEta2[nBinCent];

  //********mean pt *******//
  TProfile2D *tp2_meanPt_hadrons[nEtaGapHadrons];

	if(Check_histo==true){
	  // check histograms
	  for(Int_t l = 0; l < 2; l++) {
	    //loop by cent                                                                                                  
	    for(Int_t c = 0; c < nBinCent; c++) {
	      //loop by eta-gap 
	      for(Int_t i = 0; i < nEtaGapHadrons; i++) {
	        //loop by VtxZ bins
	        for(Int_t z=0; z < nBinVtxZ_Hadrons; z++){
	          //historam of Q-vectors for eta-gap 
	          h_Qx2[l][i][c][z] = new TH1D(Form("h_Qx2%s%scent%iVtxZ%i",direction[l],NameEtaHadrons[i],c,z),Form("Q_{x} for #psi_{2} %s %s cent %i VtxZ %i;Q_{x}; Entries",direction[l],NameEtaHadrons[i],c,z),300,-1.5,1.5);
	          h_Qy2[l][i][c][z] = new TH1D(Form("h_Qy2%s%scent%iVtxZ%i",direction[l],NameEtaHadrons[i],c,z),Form("Q_{y} for #psi_{2} %s %s cent %i VtxZ %i;Q_{y}; Entries",direction[l],NameEtaHadrons[i],c,z),300,-1.5,1.5);
	          h_Qx3[l][i][c][z] = new TH1D(Form("h_Qx3%s%scent%iVtxZ%i",direction[l],NameEtaHadrons[i],c,z),Form("Q_{x} for #psi_{3} %s %s cent %i VtxZ %i;Q_{x}; Entries",direction[l],NameEtaHadrons[i],c,z),300,-1.5,1.5);
	          h_Qy3[l][i][c][z] = new TH1D(Form("h_Qy3%s%scent%iVtxZ%i",direction[l],NameEtaHadrons[i],c,z),Form("Q_{y} for #psi_{3} %s %s cent %i VtxZ %i;Q_{y}; Entries",direction[l],NameEtaHadrons[i],c,z),300,-1.5,1.5);

	          //historam of event planes for east eta-gap 
	          h_Psi2[l][i][c][z] = new TH1D(Form("h_Psi2%s%scent%iVtxZ%i",direction[l],NameEtaHadrons[i],c,z),Form("#psi_{2} %s %s cent %i;#psi_{2}; VtxZ %i",direction[l],NameEtaHadrons[i],c,z),100,-0.05,3.2); 
	          h_Psi3[l][i][c][z] = new TH1D(Form("h_Psi3%s%scent%iVtxZ%i",direction[l],NameEtaHadrons[i],c,z),Form("#psi_{3} %s %s cent %i;#psi_{3}; VtxZ %i",direction[l],NameEtaHadrons[i],c,z),100,-0.05,2.15);
	        }// loop by VtxZ
	      }// loop by n  
	    }// loop by cent
	  }// loop by direstion 
	}
  
	if(mode_raw==true){
		// loop by direction 
		for(Int_t l = 0; l < 2; l++) {
		   //loop by systematics value
		  for(Int_t i = 0; i < nEtaGapHadrons; i++) {
		    //loop by VtxZ
		    for(Int_t z = 0; z < nBinVtxZ_Hadrons; z++){
		      //TProfile2D for recentering 
		      tp2_fill_Qx2[l][i][z] = new TProfile2D(Form("tp2_Qx2%s%sVtxZ%i",direction[l],NameEtaHadrons[i],z),Form("<Q_{x}> for #psi_{2} %s %s VtxZ %i ;RunID;cent",direction[l],NameEtaHadrons[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
		      tp2_fill_Qy2[l][i][z] = new TProfile2D(Form("tp2_Qy2%s%sVtxZ%i",direction[l],NameEtaHadrons[i],z),Form("<Q_{y}> for #psi_{2} %s %s VtxZ %i ;RunID;cent",direction[l],NameEtaHadrons[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
		      tp2_fill_Qx3[l][i][z] = new TProfile2D(Form("tp2_Qx3%s%sVtxZ%i",direction[l],NameEtaHadrons[i],z),Form("<Q_{x}> for #psi_{3} %s %s VtxZ %i ;RunID;cent",direction[l],NameEtaHadrons[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
		      tp2_fill_Qy3[l][i][z] = new TProfile2D(Form("tp2_Qy3%s%sVtxZ%i",direction[l],NameEtaHadrons[i],z),Form("<Q_{y}> for #psi_{3} %s %s VtxZ %i ;RunID;cent",direction[l],NameEtaHadrons[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
		    }// loop by VtxZ
		  }//loop by systematics value (Eta)
		}// loop by direction
	}

  	if( mode_rec==true ) { 

	  //loop by systematics value
	    for(Int_t i = 0; i < nEtaGapHadrons; i++) {
	      for(Int_t j = 0; j < 4; j++) {
	        for(Int_t z=0; z<nBinVtxZ_Hadrons; z++){
	          tp2_sinPsi2East[i][j][z] = new TProfile2D(Form("tp2_%isinPsi2%s%sVtxZ%i",j+1,direction[0],NameEtaHadrons[i],z),Form("<sin(%i*2#psi_{2})> %s %s VtxZ %i ",j+1,direction[0],NameEtaHadrons[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
	          tp2_cosPsi2East[i][j][z] = new TProfile2D(Form("tp2_%icosPsi2%s%sVtxZ%i",j+1,direction[0],NameEtaHadrons[i],z),Form("<cos(%i*2#psi_{2})> %s %s VtxZ %i ",j+1,direction[0],NameEtaHadrons[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
	          tp2_sinPsi3East[i][j][z] = new TProfile2D(Form("tp2_%isinPsi3%s%sVtxZ%i",j+1,direction[0],NameEtaHadrons[i],z),Form("<sin(%i*3#psi_{3})> %s %s VtxZ %i ",j+1,direction[0],NameEtaHadrons[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
	          tp2_cosPsi3East[i][j][z] = new TProfile2D(Form("tp2_%icosPsi3%s%sVtxZ%i",j+1,direction[0],NameEtaHadrons[i],z),Form("<cos(%i*3#psi_{3})> %s %s VtxZ %i ",j+1,direction[0],NameEtaHadrons[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);          
	        
	          tp2_sinPsi2West[i][j][z] = new TProfile2D(Form("tp2_%isinPsi2%s%sVtxZ%i",j+1,direction[1],NameEtaHadrons[i],z),Form("<sin(%i*2#psi_{2})> %s %s VtxZ %i ",j+1,direction[1],NameEtaHadrons[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
	          tp2_cosPsi2West[i][j][z] = new TProfile2D(Form("tp2_%icosPsi2%s%sVtxZ%i",j+1,direction[1],NameEtaHadrons[i],z),Form("<cos(%i*2#psi_{2})> %s %s VtxZ %i ",j+1,direction[1],NameEtaHadrons[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
	          tp2_sinPsi3West[i][j][z] = new TProfile2D(Form("tp2_%isinPsi3%s%sVtxZ%i",j+1,direction[1],NameEtaHadrons[i],z),Form("<sin(%i*3#psi_{3})> %s %s VtxZ %i ",j+1,direction[1],NameEtaHadrons[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);
	          tp2_cosPsi3West[i][j][z] = new TProfile2D(Form("tp2_%icosPsi3%s%sVtxZ%i",j+1,direction[1],NameEtaHadrons[i],z),Form("<cos(%i*3#psi_{3})> %s %s VtxZ %i ",j+1,direction[1],NameEtaHadrons[i],z),RunIdRange,RunIdMin.at(energy),RunIdMax.at(energy),nBinCent,0,nBinCent);          
	        }
	      }
	    }


      //FileRec = new TFile(Form("%s/OUT/%iGeV/raw_%iGeV_Hadrons_2binYesTof_dcaEP2_hist.root", path,energy, energy),"READ");
	    FileRec = new TFile(Form("%s/OUT/%iGeV/raw_%iGeV_Hadrons.root", path,energy, energy),"READ");
	    FileRec->cd();

	    for(Int_t dir = 0; dir < 2; dir++) {
	      for(Int_t i = 0; i < nEtaGapHadrons; i++) {
	        for(Int_t z=0; z<nBinVtxZ_Hadrons; z++){
	          
	          tp_read_profile = (TProfile2D*) FileRec -> Get( Form("tp2_Qx2%s%sVtxZ%i",direction[dir],NameEtaHadrons[i],z) );
	          for(Int_t c = 0; c < nBinCent; c++){
	            for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
	              binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
	              binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
	              if (binCont == 0. || binEntr == 0.) continue;
	              mp_Qx2[dir][i][z][c][(Double_t)run] = binCont;
	            }
	          }
	          delete tp_read_profile;    

	          tp_read_profile = (TProfile2D*) FileRec -> Get( Form("tp2_Qy2%s%sVtxZ%i",direction[dir],NameEtaHadrons[i],z) );
	          for(Int_t c = 0; c < nBinCent; c++){
	            for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
	              binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
	              binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
	              if (binCont == 0. || binEntr == 0.) continue;
	              mp_Qy2[dir][i][z][c][(Double_t)run] = binCont;
	            }
	          }
	          delete tp_read_profile;  

	          tp_read_profile = (TProfile2D*) FileRec -> Get( Form("tp2_Qx3%s%sVtxZ%i",direction[dir],NameEtaHadrons[i],z) );
	          for(Int_t c = 0; c < nBinCent; c++){
	            for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
	              binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
	              binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
	              if (binCont == 0. || binEntr == 0.) continue;
	              mp_Qx3[dir][i][z][c][(Double_t)run] = binCont;
	            }
	          }
	          delete tp_read_profile;  

	          tp_read_profile = (TProfile2D*) FileRec -> Get( Form("tp2_Qy3%s%sVtxZ%i",direction[dir],NameEtaHadrons[i],z) );
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

	    for(Int_t i = 0; i < nEtaGapHadrons; i++) {
	      
	      tp_SqRes2[i] = new TProfile(Form("tp_SqRes2TPC%s",NameEtaHadrons[i]),Form("Resolution^{2} for v_{2} %s ",NameEtaHadrons[i]),nBinCent,0,nBinCent);
	      tp_SqRes3[i] = new TProfile(Form("tp_SqRes3TPC%s",NameEtaHadrons[i]),Form("Resolution^{2} for v_{3} %s ",NameEtaHadrons[i]),nBinCent,0,nBinCent);
	      
	    }
	    
	    for(Int_t c=0; c<nBinCent; c++){
  			hEta[c] = new TH1D(Form("hEtaCentBin%i", c),Form("Track pseudorapidity for Cent bin %i;#eta;Entries", c), 110, -1.1, 1.1 );
  			hEta2[c] = new TH1D(Form("hEtaCentBin%ipt2", c),Form("Track pseudorapidity for Cent bin %i;#eta;Entries", c), 110, -1.1, 1.1 );
  		}

	    for(Int_t i = 0; i < nEtaGapHadrons; i++) {
	    
	      tp2_meanPt_hadrons[i] = new TProfile2D(Form("tp_meanPtHadrons%s",NameEtaHadrons[i]),Form("Mean p_{t} for bins v_{2} %s ; bin; p_{t} [GeV/c]",NameEtaHadrons[i]),200,0.2,10.2,nBinCent,0,nBinCent);

	      tp_v2hadrTPC_cent[i] = new TProfile(Form("tp_v2ewTPCHadrons%sCent",NameEtaHadrons[i]),Form("v_{2} of cent by TPC %s;cent bin",NameEtaHadrons[i]),nBinCent,0,nBinCent);
	      tp_v3hadrTPC_cent[i] = new TProfile(Form("tp_v3ewTPCHadrons%sCent",NameEtaHadrons[i]),Form("v_{3} of cent by TPC %s;cent bin",NameEtaHadrons[i]),nBinCent,0,nBinCent);

	      for(Int_t c = 0; c < nBinCent; c++){
	    	tp2_v2hadrTPC_eta[i][c] = new TProfile2D(Form("tp_v2ewTPCHadrons%sEtaCentBin%i",NameEtaHadrons[i],c),Form("v_{2} of Eta by TPC %s CentBin %i;#eta; p_{t} [GeV/c] ",NameEtaHadrons[i], c),100,-1.0,1.0, 100, 0.2, 5.2);
	      	tp2_v3hadrTPC_eta[i][c] = new TProfile2D(Form("tp_v3ewTPCHadrons%sEtaCentBin%i",NameEtaHadrons[i],c),Form("v_{3} of Eta by TPC %s CentBin %i;#eta; p_{t} [GeV/c] ",NameEtaHadrons[i], c),100,-1.0,1.0, 100, 0.2, 5.2);
	      }

	      tp2_v2hadrTPC[i] = new TProfile2D(Form("tp_v2ewTPCHadrons%s",NameEtaHadrons[i]),Form("v_{2} of p_{t} and cent by TPC %s;p_{t} [GeV/c];cent",NameEtaHadrons[i]),200,0.2,10.2,nBinCent,0,nBinCent);
	      tp2_v3hadrTPC[i] = new TProfile2D(Form("tp_v3ewTPCHadrons%s",NameEtaHadrons[i]),Form("v_{3} of p_{t} and cent by TPC %s;p_{t} [GeV/c];cent",NameEtaHadrons[i]),200,0.2,10.2,nBinCent,0,nBinCent);

	    }// for(Int_t i = 0; i < n; i++)


	    FileRec = new TFile(Form("%s/OUT/%iGeV/raw_%iGeV_Hadrons.root",path,energy, energy),"READ");
	    FileRec -> cd();
	    for(Int_t dir = 0; dir < 2; dir++) {
	      for(Int_t i = 0; i < nEtaGapHadrons; i++) {
	        for(Int_t z = 0; z < nBinVtxZ_Hadrons; z++){
	          
	          tp_read_profile = (TProfile2D*) FileRec -> Get( Form("tp2_Qx2%s%sVtxZ%i",direction[dir],NameEtaHadrons[i],z) );
	          
	          for(Int_t c = 0; c < nBinCent; c++){
	            for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
	              binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
	              binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
	              if (binCont == 0. || binEntr == 0.) continue;
	              mp_Qx2[dir][i][z][c][(Double_t)run] = binCont;
	            }
	          }
	          delete tp_read_profile;    

	          tp_read_profile = (TProfile2D*) FileRec -> Get( Form("tp2_Qy2%s%sVtxZ%i",direction[dir],NameEtaHadrons[i],z) );
	          for(Int_t c = 0; c < nBinCent; c++){
	            for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
	              binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
	              binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
	              if (binCont == 0. || binEntr == 0.) continue;
	              mp_Qy2[dir][i][z][c][(Double_t)run] = binCont;
	            }
	          }
	          delete tp_read_profile;  

	          tp_read_profile = (TProfile2D*) FileRec -> Get( Form("tp2_Qx3%s%sVtxZ%i",direction[dir],NameEtaHadrons[i],z) );
	          for(Int_t c = 0; c < nBinCent; c++){
	            for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
	              binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
	              binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
	              if (binCont == 0. || binEntr == 0.) continue;
	              mp_Qx3[dir][i][z][c][(Double_t)run] = binCont;
	            }
	          }
	          delete tp_read_profile;  

	          tp_read_profile = (TProfile2D*) FileRec -> Get( Form("tp2_Qy3%s%sVtxZ%i",direction[dir],NameEtaHadrons[i],z) );
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

	    FileFlow = new TFile(Form("%s/OUT/%iGeV/rec_%iGeV_Hadrons.root", path,energy, energy),"READ");
	    FileFlow -> cd();

	    for(Int_t l = 0; l < 2; l++) {
	      for(Int_t i = 0; i < nEtaGapHadrons; i++) {
	        for(Int_t j = 0; j < 4; j++) {
	          for(Int_t z = 0; z < nBinVtxZ_Hadrons; z++){

	            tp_read_profile =  (TProfile2D*) FileFlow -> Get( Form("tp2_%isinPsi2%s%sVtxZ%i",j+1,direction[l],NameEtaHadrons[i],z) );
	            
	            for(Int_t c = 0; c < nBinCent; c++){
	              for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
	                binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
	                binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
	                if (binCont == 0. || binEntr == 0.) continue;
	                mp_sinPsi2[l][i][j][z][c][(Double_t)run] = binCont;
	              }
	            }
	            delete tp_read_profile; 

	            tp_read_profile =  (TProfile2D*) FileFlow -> Get( Form("tp2_%icosPsi2%s%sVtxZ%i",j+1,direction[l],NameEtaHadrons[i],z) );
	            for(Int_t c = 0; c < nBinCent; c++){
	              for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
	                binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
	                binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
	                if (binCont == 0. || binEntr == 0.) continue;
	                mp_cosPsi2[l][i][j][z][c][(Double_t)run] = binCont;
	              }
	            }
	            delete tp_read_profile;

	            tp_read_profile =  (TProfile2D*) FileFlow -> Get( Form("tp2_%isinPsi3%s%sVtxZ%i",j+1,direction[l],NameEtaHadrons[i],z) );
	            for(Int_t c = 0; c < nBinCent; c++){
	              for(Int_t run = RunIdMin.at(energy); run < RunIdMax.at(energy); run++){
	                binCont = tp_read_profile->GetBinContent(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
	                binEntr = tp_read_profile->GetBinEntries(tp_read_profile->FindBin((Double_t)run, (Double_t)c));
	                if (binCont == 0. || binEntr == 0.) continue;
	                mp_sinPsi3[l][i][j][z][c][(Double_t)run] = binCont;
	              }
	            }
	            delete tp_read_profile;

	            tp_read_profile =  (TProfile2D*) FileFlow -> Get( Form("tp2_%icosPsi3%s%sVtxZ%i",j+1,direction[l],NameEtaHadrons[i],z) );
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

    // Cut for Tof Mached
    Int_t nTracksTOF = dst->numberOfTracks();
    Bool_t check_TofMatch=false;
    Int_t number_tof = 0;
    for(Int_t iTrk=0; iTrk<nTracksTOF; iTrk++) {
      StFemtoTrack *femtoTrack = dst->track(iTrk);
      if ( !femtoTrack ) continue;
      if ( femtoTrack->isTofTrack()){
        number_tof++;
      }
      if(number_tof > 4){
        check_TofMatch = true;
        iTrk = nTracksTOF;
      }
    }
    if(check_TofMatch == false)continue;


    if( !isGoodEvent( event,  energy, mode_flow) ) continue;

    RunID = event -> runId();
    
    if( nBinCent == 9 ){
      cent = event -> cent9();
    }
    if( nBinCent == 16 ){
      cent = event -> cent16();
    }
    
    Int_t signVtxZ = GetBinVtxZ(event,energy); // 0 - VtxZ < 0. ; 1 - VtxZ > 0.
    if( signVtxZ == -1 ) continue;

    for(Int_t dir = 0; dir < 2; dir++) {
      for(Int_t s = 0; s < nEtaGapHadrons; s++) { 
        mult[dir][s] = 0;
        Q2vecTPC[dir][s].Set(0.,0.);
        Q3vecTPC[dir][s].Set(0.,0.);
      }
    }

    // Track analysis
    Int_t nTracks = dst->numberOfTracks();

    // Track loop
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

      // Retrieve i-th femto track
      StFemtoTrack *femtoTrack = dst->track(iTrk);
      
      if( !isGoodTrackEP(event, femtoTrack, energy) ) continue; 
      
      Int_t EtaDir = GetEtaDirection(femtoTrack); // 0 - east; 1 - west ; -1 - bad
      Int_t EtaBin = GetBinEta(femtoTrack);
      if( EtaDir == -1 || EtaBin == -1) continue;
      
      for(Int_t i = 0; i <= EtaBin; i++ ){
        Q2vecTPC[EtaDir][i] = Q2vecTPC[EtaDir][i] + femtoTrack -> pt() * CalculateQvector(2,femtoTrack);
        Q3vecTPC[EtaDir][i] = Q3vecTPC[EtaDir][i] + femtoTrack -> pt() * CalculateQvector(3,femtoTrack);
        mult[EtaDir][i] = mult[EtaDir][i] + 1;
      }

    } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

    MeanQ2.Set(0.,0.);
    MeanQ3.Set(0.,0.);

    for(Int_t dir = 0; dir < 2; dir++) {
      for(Int_t i = 0; i < nEtaGapHadrons; i++) {
        if(mode_rec == true || mode_flow == true){
          MeanQ2.Set( mp_Qx2[dir][i][signVtxZ][cent][(Double_t)RunID] , mp_Qy2[dir][i][signVtxZ][cent][(Double_t)RunID]);
          MeanQ3.Set( mp_Qx3[dir][i][signVtxZ][cent][(Double_t)RunID] , mp_Qy3[dir][i][signVtxZ][cent][(Double_t)RunID]);
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
        for(Int_t i = 0; i < nEtaGapHadrons; i++) {
          if( Q2vecTPC[dir][i].Mod() != 0. && Q3vecTPC[dir][i].Mod() != 0.) {
            
            Psi2TPC[dir][i] = Q2vecTPC[dir][i].Phi() / 2.0;
            Psi3TPC[dir][i] = Q3vecTPC[dir][i].Phi() / 3.0;
            
            for(Int_t k = 0; k < 4; k++) {
              
              sinPsi2 = mp_sinPsi2[dir][i][k][signVtxZ][(int)cent][(Double_t)RunID];
              cosPsi2 = mp_cosPsi2[dir][i][k][signVtxZ][(int)cent][(Double_t)RunID];
              sinPsi3 = mp_sinPsi3[dir][i][k][signVtxZ][(int)cent][(Double_t)RunID];
              cosPsi3 = mp_cosPsi3[dir][i][k][signVtxZ][(int)cent][(Double_t)RunID];

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
        }// for(Int_t i = 0; i < nEtaGapHadrons; i++){}
      }// for(Int_t dir = 0; dir < 3; dir++){} 
    
      for(Int_t i = 0; i < nEtaGapHadrons; i++) {
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
        if( EtaDir == -1 || EtaBin == -1) continue;

        weight = GetWeight(femtoTrack);

        Phi = femtoTrack -> phi();
        pt = femtoTrack -> pt();

        v2 = 0.;
        v3 = 0.;
        if(pt < 3.2){
	        hEta[cent]->Fill( (Double_t)femtoTrack -> eta() );
        }
        if(pt < 2.0){
        	hEta2[cent]->Fill( (Double_t)femtoTrack -> eta() );
        }

        for(Int_t s = 0; s <= EtaBin; s++) { 

          v2 = TMath::Cos( 2.0*(Phi - Psi2TPC[TMath::Abs(EtaDir-1)][s]) );
          v3 = TMath::Cos( 3.0*(Phi - Psi3TPC[TMath::Abs(EtaDir-1)][s]) );

          tp_v2hadrTPC_cent[s] -> Fill((Double_t)cent, v2);
          tp_v3hadrTPC_cent[s] -> Fill((Double_t)cent, v3);

          tp2_v2hadrTPC_eta[s][cent] -> Fill((Double_t)femtoTrack -> eta(), pt, v2, weight);
          tp2_v3hadrTPC_eta[s][cent] -> Fill((Double_t)femtoTrack -> eta(), pt, v3, weight);

          tp2_v2hadrTPC[s] -> Fill(pt, (Double_t)cent, v2, weight);
          tp2_v3hadrTPC[s] -> Fill(pt, (Double_t)cent, v3, weight);
          tp2_meanPt_hadrons[s] -> Fill(pt,(Double_t)cent,pt);

        }// for(Int_t eta = 0; eta < n; eta++)
      } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)
    }// if(mode_flow == true){}

    if(Check_histo==true){
		for(Int_t dir = 0; dir < 2; dir ++) {
		  for(Int_t i = 0; i < nEtaGapHadrons; i++) {
		    if( Q2vecTPC[dir][i].Mod() != 0. && Q3vecTPC[dir][i].Mod() != 0.) {
		      if(mode_raw == true || mode_rec == true){
		        h_Qx2[dir][i][cent][signVtxZ] -> Fill( (Double_t)Q2vecTPC[dir][i].X() );
		        h_Qy2[dir][i][cent][signVtxZ] -> Fill( (Double_t)Q2vecTPC[dir][i].Y() );
		        h_Qx3[dir][i][cent][signVtxZ] -> Fill( (Double_t)Q3vecTPC[dir][i].X() );
		        h_Qy3[dir][i][cent][signVtxZ] -> Fill( (Double_t)Q3vecTPC[dir][i].Y() );
		        
		        h_Psi2[dir][i][cent][signVtxZ] -> Fill( Q2vecTPC[dir][i].Phi() / 2.0 );
		        h_Psi3[dir][i][cent][signVtxZ] -> Fill( Q3vecTPC[dir][i].Phi() / 3.0 );
		      }
		      if(mode_flow == true){
		        h_Psi2[dir][i][cent][signVtxZ] -> Fill( Psi2TPC[dir][i] );
		        h_Psi3[dir][i][cent][signVtxZ] -> Fill( Psi3TPC[dir][i] );
		      } 
		    }
		  }
		}
	}

    if(mode_raw == true ){
      for(Int_t dir = 0; dir < 2; dir ++) {
        for(Int_t i = 0; i < nEtaGapHadrons; i++) {
          if( Q2vecTPC[dir][i].Mod() != 0. && Q3vecTPC[dir][i].Mod() != 0.) {
            tp2_fill_Qx2[dir][i][signVtxZ] -> Fill( (Double_t)RunID, (Double_t)cent, (Double_t)Q2vecTPC[dir][i].X() );
            tp2_fill_Qy2[dir][i][signVtxZ] -> Fill( (Double_t)RunID, (Double_t)cent, (Double_t)Q2vecTPC[dir][i].Y() );
            tp2_fill_Qx3[dir][i][signVtxZ] -> Fill( (Double_t)RunID, (Double_t)cent, (Double_t)Q3vecTPC[dir][i].X() );
            tp2_fill_Qy3[dir][i][signVtxZ] -> Fill( (Double_t)RunID, (Double_t)cent, (Double_t)Q3vecTPC[dir][i].Y() );
          }
        }// for(Int_t i = 0; i < nEtaGapHadrons; i++){}  
      }// for(Int_t dir = 0; dir < 2; dir++){}
    }// if(mode_raw==true){}

    if(mode_rec == true ){
      for(Int_t i = 0; i < nEtaGapHadrons; i++) {
        for(Int_t j = 0; j < 4; j++) {
          if( Q2vecTPC[0][i].Mod() != 0. && Q3vecTPC[0][i].Mod() != 0.) {
            tp2_sinPsi2East[i][j][signVtxZ] -> Fill((Double_t)RunID, (Double_t)cent, TMath::Sin( (Double_t)(j+1)*Q2vecTPC[0][i].Phi() ) );
            tp2_cosPsi2East[i][j][signVtxZ] -> Fill((Double_t)RunID, (Double_t)cent, TMath::Cos( (Double_t)(j+1)*Q2vecTPC[0][i].Phi() ) );
            tp2_sinPsi3East[i][j][signVtxZ] -> Fill((Double_t)RunID, (Double_t)cent, TMath::Sin( (Double_t)(j+1)*Q3vecTPC[0][i].Phi() ) ); 
            tp2_cosPsi3East[i][j][signVtxZ] -> Fill((Double_t)RunID, (Double_t)cent, TMath::Cos( (Double_t)(j+1)*Q3vecTPC[0][i].Phi() ) );
          }
          if( Q2vecTPC[1][i].Mod() != 0. && Q3vecTPC[1][i].Mod() != 0.) {
            tp2_sinPsi2West[i][j][signVtxZ] -> Fill((Double_t)RunID, (Double_t)cent, TMath::Sin( (Double_t)(j+1)*Q2vecTPC[1][i].Phi() ) );
            tp2_cosPsi2West[i][j][signVtxZ] -> Fill((Double_t)RunID, (Double_t)cent, TMath::Cos( (Double_t)(j+1)*Q2vecTPC[1][i].Phi() ) );
            tp2_sinPsi3West[i][j][signVtxZ] -> Fill((Double_t)RunID, (Double_t)cent, TMath::Sin( (Double_t)(j+1)*Q3vecTPC[1][i].Phi() ) ); 
            tp2_cosPsi3West[i][j][signVtxZ] -> Fill((Double_t)RunID, (Double_t)cent, TMath::Cos( (Double_t)(j+1)*Q3vecTPC[1][i].Phi() ) );
          }
        }//for(Int_t j =0; j < 4; j++){}
      }// for(Int_t i = 0; i < nEtaGapHadrons; i++){}  
    }// if(mode_raw==true){}



  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

  outFile->Write();
  outFile->Close();

  femtoReader->Finish();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
	          << std::endl;
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
  /*
  std::vector<Int_t> Bad = BadRuns.at(_energy);
  if(BadRunIdKeyFlowStage==true){
    if( std::find(Bad.begin(), Bad.end(), event -> runId()) != Bad.end() ) return false;
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
  if( track -> gDCA(pVtx).Mag() > CutDCAHadronsEP.at(_energy) ) return false;    
  if( TMath::Abs( track -> eta() ) > 1.0 ) return false; 
  if( track -> nHits() < CutnHits) return false;
  if( track -> pt() < CutPtEPMin_Hadrons) return false;
  if( track -> pt() > CutPtEPMax_Hadrons) return false; 
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
  if( track -> gDCA(pVtx).Mag() > CutDCAHadronsFlow.at(_energy) ) return false;
  if( TMath::Abs( track -> eta() ) > 1.0 ) return false; 
  if( track -> nHits() < CutnHits) return false;
  if( track -> pt() < CutPtEPMin_Hadrons) return false;
  if( track -> pt() > CutPtFlowMax_Hadrons) return false; 
  if( track -> p() < CutPtEPMin_Hadrons) return false;
  if( track -> p() > CutPtFlowMax_Hadrons) return false; 
  if( ( (Double_t)track -> nHits() )/( (Double_t)track -> nHitsPoss() )  < CutnHitsRatio ) return false;

  return true; 
}// isGoodTrack(){}
