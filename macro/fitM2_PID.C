// C++ headers
#include <vector>
#include <iostream>
#include <fstream> 
#include <map>
#include <string>

// ROOT headers
#include "TGraph.h"
#include "TPad.h"
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

// Constant
#include "./Constants.h"

//_______________________________________________________________________________________________________________________________________________________________________
TF1 *fitFun1DGaus(TH1D *histo, Double_t cons, Double_t mean, Double_t sigma, Int_t numberFixPar, Double_t min, Double_t max, Int_t flag_fit);
TF1 *fitFunTwo1DGaus(TH1D *histo, TF1 *FirstGaus, TF1 *SecondGaus, Double_t min, Double_t max );
TF1 *fitFunThree1DGausNew(TH1D *histo,TF1 *FirstGaus, TF1 *SecondGaus, TF1 *ThirdGaus, Double_t min, Double_t max );
//_______________________________
void FitSquareOfMass(TFile *fileM2, TFile *fileOut);
TF1 *FitSquareOfMassInPtBin(TFile *fileM2, TFile *fileOut, Int_t ptBin, Int_t charge);


void fitM2_PID(const Char_t *inFileNameWithHisto = "iTest.root", const Char_t *inFileNameWithFitFun = "iTest2.root", const Char_t *outFileName = "oTest.root", const Int_t energy = 39, const Int_t step = 30){

  	TFile *inFileWithHisto;
  	TFile *inFileWithFitFun;
  	TFile *outFile;

  	if(step==1){
		
		inFileWithHisto = new TFile(Form("%s",inFileNameWithHisto), "READ");
  		outFile = new TFile(Form("%s",outFileName), "RECREATE");

  		FitSquareOfMass(inFileWithHisto, outFile);

  		outFile->Close();
  		inFileWithHisto->Close();

  	}  	

}


TF1 *fitFun1DGaus(TH1D *histo, Double_t cons, Double_t mean, Double_t sigma, Int_t numberFixPar, Double_t min, Double_t max, Int_t flag_fit){
	
	TF1 *gaus = new TF1("","gaus",min,max);
	gaus->SetParameter(0,cons);
	gaus->SetParameter(1,mean);
	gaus->SetParameter(2,sigma);
	if(flag_fit==0){
		histo->Fit(gaus,"RM");
	}
	return gaus;

}

TF1 *fitFunTwo1DGaus(TH1D *histo, TF1 *FirstGaus, TF1 *SecondGaus, Double_t min, Double_t max ){

	TF1 *gaus = new TF1("Two1DGaus","gaus(0)+gaus(3)",min,max);
	gaus->SetParameter(0,FirstGaus->GetParameter(0));
	gaus->SetParameter(1,FirstGaus->GetParameter(1));
	gaus->SetParameter(2,FirstGaus->GetParameter(2));
	gaus->SetParameter(3,SecondGaus->GetParameter(0));
	gaus->SetParameter(4,SecondGaus->GetParameter(1));
	gaus->SetParameter(5,SecondGaus->GetParameter(2));
	histo->Fit(gaus,"RM");
	return gaus;

}

TF1 *fitFunThree1DGausNew(TH1D *histo,TF1 *FirstGaus, TF1 *SecondGaus, TF1 *ThirdGaus, Double_t min, Double_t max ){

	TF1 *gaus = new TF1("","gaus(0)+gaus(3)+gaus(6)",min,max);
	gaus->SetParameter(0,FirstGaus->GetParameter(0));
	gaus->SetParameter(1,FirstGaus->GetParameter(1));
	gaus->SetParameter(2,FirstGaus->GetParameter(2));
	gaus->SetParameter(3,SecondGaus->GetParameter(0));
	gaus->SetParameter(4,SecondGaus->GetParameter(1));
	gaus->SetParameter(5,SecondGaus->GetParameter(2));
	gaus->SetParameter(6,ThirdGaus->GetParameter(0));
	gaus->SetParameter(7,ThirdGaus->GetParameter(1));
	gaus->SetParameter(8,ThirdGaus->GetParameter(2));
	histo->Fit(gaus,"RM");
	return gaus;

}

void FitSquareOfMass(TFile *fileM2, TFile *fileOut){

	gStyle->SetOptStat("n");
    gStyle->SetOptFit(1);

    Int_t n_pti = (int)ptBinRange.size();
    
    Double_t mean[3][2][100]={0.};
    Double_t sigma[3][2][100]={0.};
    Double_t x_pt[2][100]={0.};
    Double_t err[100]={0.};

    TF1 *FitGaus[2][100];

	for(Int_t ch=0; ch<2; ch++){
		for(Int_t pti=0; pti<n_pti-1; pti++){
			FitGaus[ch][pti] = FitSquareOfMassInPtBin(fileM2, fileOut, pti, ch);
			for(Int_t par=0; par<3; par++){
				mean[par][ch][pti]  = FitGaus[ch][pti]->GetParameter(1 + 3 * par);
				sigma[par][ch][pti] = FitGaus[ch][pti]->GetParameter(2 + 3 * par);
			}
			x_pt[ch][pti] = (ptBinRange[pti] + ptBinRange[pti+1]) / 2.0;
		}
	}
	
	TGraphErrors *graphMean[2][3];
    TGraphErrors *graphSigma[2][3];
	
	for(Int_t par=0; par<3; par++){
		for(Int_t ch=0; ch<2; ch++){
			graphMean[ch][par] = new TGraphErrors((int)ptBinRange.size()-1, x_pt[ch], mean[par][ch], err, err);
			graphMean[ch][par]->SetName(Form("gr_Mean%s_ch%i",particles[par],ch));
			graphMean[ch][par]->SetTitle(Form("Mean m^{2} (%s) vs p_{T} bin",partLateX[2*par+ch]));
			graphMean[ch][par]->GetXaxis()->SetTitle("#font[42]{p_{T}[GeV/c]}");
			graphMean[ch][par]->GetYaxis()->SetTitle(Form("#font[42]{#mu(m^{2})(%s)}",partLateX[2*par+ch]));

			graphMean[ch][par]->Write();

			graphSigma[ch][par] = new TGraphErrors((int)ptBinRange.size()-1, x_pt[ch], sigma[par][ch], err, err);
			graphSigma[ch][par]->SetName(Form("gr_Sigma%s_ch%i",particles[par],ch));
			graphSigma[ch][par]->SetTitle(Form("Sigma m^{2} (%s) vs p_{T} bin",partLateX[2*par+ch]));
			graphSigma[ch][par]->GetXaxis()->SetTitle("#font[42]{p_{T}[GeV/c]}");
			graphSigma[ch][par]->GetYaxis()->SetTitle(Form("#font[42]{#sigma(m^{2})(%s)}",partLateX[2*par+ch]));

			graphSigma[ch][par]->Write();

		}
	}

	TF1 *f_mean[2][3];
	TF1 *f_sigma[2][3];
	Double_t min[4]={0.2, 0.3, 0.4};
    Double_t max[4]={1.6, 1.5, 2.1};
    
    Double_t range_mean[2][3]={{-0.051, 0.0, 0.82},
    						   {0.04, 0.31, 0.892}};
    Double_t range_sigma[2][3]={{0.0, 0.0, 0.0},
    						    {0.241, 0.31, 0.31}};

   	for(Int_t par=0; par<3; par++){						    
	    for(Int_t ch=0; ch<2; ch++){

	    	if(par==0 || par==1){
		    	f_mean[ch][par] = new TF1(Form("f1_meanM2_%s_ch%i",particles[par],ch),"[0] + [1]*x + [2]*x*x", min[par] , max[par] );
		    	//graphMean[ch][par]->GetYaxis()->SetRangeUser( range_mean[0][par],range_mean[1][par]);
				graphMean[ch][par]->Fit(f_mean[ch][par],"RM");
	    	}

	    	if(par==2){
		    	f_mean[ch][par] = new TF1(Form("f1_meanM2_%s_ch%i",particles[par],ch),"[0] + [1]*x + [2]*x*x + [3]*x*x*x", min[par] , max[par] );
		    	//graphMean[ch][par]->GetYaxis()->SetRangeUser( range_mean[0][par],range_mean[1][par]);
				graphMean[ch][par]->Fit(f_mean[ch][par],"RM");
	    	}
			
			f_sigma[ch][par] = new TF1(Form("f1_sigmaM2_%s_ch%i",particles[par],ch),"[0] + [1]*x + [2]*x*x", min[par] , max[par] );
	    	//graphSigma[ch][par]->GetYaxis()->SetRangeUser( range_sigma[0][par],range_sigma[1][par]);
	    	graphSigma[ch][par]->Fit(f_sigma[ch][par],"RM");

	    	f_mean[ch][par]->Write();
	    	f_sigma[ch][par]->Write();

	    }
	}


}

TF1 *FitSquareOfMassInPtBin(TFile *fileM2, TFile *fileOut, Int_t ptBin, Int_t charge){

	gStyle->SetOptStat("n");
    gStyle->SetOptFit(1);


    TF1 *tf1_fitFun;
    TH1D *h1_m2;

    Double_t min[4]={-0.5, 0.15, 0.7, -0.5};
    Double_t max[4]={0.1, 0.4, 1.2, 0.5};
    Double_t massSqrt[3]={pion_mass*pion_mass, kaon_mass*kaon_mass, proton_mass*proton_mass};

    TH2D *h2_m2VsPt = (TH2D*)fileM2->Get(Form("h2_m2VsPt_all_ch%i", charge));
    TH1D *h1_m2_projY = (TH1D*)h2_m2VsPt->ProjectionX();

    fileOut->cd();

	h1_m2 = (TH1D*)h2_m2VsPt->ProjectionY(Form("h1_m2_charge%i_pt%i", charge,ptBin),h1_m2_projY->FindBin(ptBinRange[ptBin]),h1_m2_projY->FindBin(ptBinRange[ptBin+1]));
    h1_m2->SetTitle(Form("m^{2}, %.1f < p_{T} < %.1f [GeV/c]",ptBinRange[ptBin],ptBinRange[ptBin+1]));
    
    if(ptBinRange[ptBin] > 1.5) h1_m2->Rebin(2);

    TF1 *Three1DGaus;
    TF1 *gausParticle[3];

    gausParticle[0] = fitFun1DGaus(h1_m2, 10., massSqrt[0], 0.01, 10, -0.15, 0.1, 0);
    gausParticle[1] = fitFun1DGaus(h1_m2, 10., massSqrt[1], 0.01, 10,  0.15, 0.4, 0);
    gausParticle[2] = fitFun1DGaus(h1_m2, 10., massSqrt[2], 0.01, 10,  0.6,  1.2, 0);
    
	Three1DGaus = fitFunThree1DGausNew(h1_m2, gausParticle[0], gausParticle[1], gausParticle[2], -0.2, 1.4);
    Three1DGaus->SetName(Form("tf1_FitM2Three1DGaus_charge%i_pt%i",charge, ptBin)); 
	Three1DGaus->SetTitle(Form("Three 1D Gaus, %.1f < p_{T} < %.1f [GeV/c]",ptBinRange[ptBin],ptBinRange[ptBin+1]));

	h1_m2->Write();
	Three1DGaus->Write();

	return Three1DGaus;
			
}
