// C++ headers
#include <vector>
#include <iostream>
#include <fstream> 
#include <map>
#include <string>

// ROOT headers
#include <TStopwatch.h>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TCanvas.h"
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
#include "Constants.h"

//_______________________________________________________________________________________________________________________________________________________________________
TF1 *fitFun1DGaus(TH1D *histo, Double_t cons, Double_t mean, Double_t sigma, Int_t numberFixPar, Double_t min, Double_t max, Int_t flag_fit);
TF1 *fitFunTwo1DGaus(TH1D *histo, TF1 *FirstGaus, TF1 *SecondGaus, Double_t min, Double_t max );
TF1 *fitFunThree1DGausNew(TH1D *histo,TF1 *FirstGaus, TF1 *SecondGaus, TF1 *ThirdGaus, Double_t min, Double_t max );
//_______________________________
void FitSquareOfMass(TFile *fileM2, TFile *fileOut);
TF1 *FitSquareOfMassInPtBin(TFile *fileM2, TFile *fileOut, Int_t ptBin, Int_t charge);

Double_t readParFromFun(TFile *inFile, TString name, Double_t pt);
Double_t readParFromGraph(TFile *inFile, TString name, Double_t pt);
void FitTF1_dEdx_nSigma(TFile *fileIn, TFile *fileWithFitM2, Int_t PID_code, Int_t charge);
void FirstFitThree2x2DGaus(TFile *inFileHisto, TFile *inFileFitFun, Int_t pt);
void FitThree2x2DGaus(TFile *inFile, TFile *outFile, TH2D *histo, Int_t pt, Int_t charge);
void DrawAndSavePar(const Char_t *name, TFile *outFile);
void GetParFunction(TFile *outFile);


void FirstFitThree2x2DGausNewXY(TFile *inFileHisto, TFile *outFile, Int_t pt);
void FitThree2x2DGausNewXY(TFile *outFile, TH2D *histo, Int_t pt, Int_t charge);

void combPID_Gaus(const Char_t *inFileNameWithHisto = "iTest.root", const Char_t *inFileNameWithFitFun = "iTest1.root", const Char_t *outFileName = "oTest2.root", const Int_t energy = 39, Int_t step = 30, const Int_t pt=2){

	TStopwatch timer1;
  	timer1.Start();

  	TFile *inFileWithHisto;
  	TFile *inFileWithFitFun;
  	TFile *outFile;
  	TFile *inFileAfterFirstFit;

  	if(step==1){
		
		inFileWithHisto = new TFile(Form("%s",inFileNameWithHisto), "READ");
  		outFile = new TFile(Form("%s",inFileNameWithFitFun), "RECREATE");

  		FitSquareOfMass(inFileWithHisto, outFile);

  		outFile->Close();

  		outFile = new TFile(Form("%s",inFileNameWithFitFun), "UPDATE");

  		for(Int_t par=0; par<3; par++){
  			for(Int_t ch=0; ch<2; ch++){
  				FitTF1_dEdx_nSigma(inFileWithHisto, outFile, par, ch);
  			}
  		}

  		inFileWithHisto->Close();
  		outFile->Close();
  		//step=2;

  	}

  	if(step==2){

  		inFileWithHisto = new TFile(inFileNameWithHisto, "READ");
  		inFileWithFitFun = new TFile(inFileNameWithFitFun, "UPDATE");
		
		FirstFitThree2x2DGaus(inFileWithHisto, inFileWithFitFun, pt);	

  		inFileWithHisto->Close();
		inFileWithFitFun->Close();
  	}

  	if(step==3){

  		outFile = new TFile(outFileName, "UPDATE");

		DrawAndSavePar(inFileNameWithFitFun, outFile);	

  		outFile->Close();
  	}
  	/*
  	if(step==444){

  		outFile = new TFile(outFileName, "UPDATE");

		GetParFunction(outFile);	

  		outFile->Close();
  	}
	*/

  	if(step==4){

  		inFileWithHisto = new TFile(inFileNameWithHisto, "READ");
  		inFileWithFitFun = new TFile(outFileName, "RECREATE");
		
		FirstFitThree2x2DGausNewXY(inFileWithHisto, inFileWithFitFun, pt);	

  		inFileWithHisto->Close();
		inFileWithFitFun->Close();
  	}


  	timer1.Stop();
  	timer1.Print();

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
    Double_t max[4]={1.4, 1.45, 2.1};
    
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
				graphMean[ch][par]->Fit(f_mean[ch][par],"RM");
	    	}

	    	if(par==2){
		    	f_mean[ch][par] = new TF1(Form("f1_meanM2_%s_ch%i",particles[par],ch),"[0] + [1]*x + [2]*x*x + [3]*x*x*x", min[par] , max[par] );
		    	//graphMean[ch][par]->GetYaxis()->SetRangeUser( range_mean[0][par],range_mean[1][par]);
				graphMean[ch][par]->Fit(f_mean[ch][par],"RM");
				graphMean[ch][par]->Fit(f_mean[ch][par],"RM");
	    	}
			
			f_sigma[ch][par] = new TF1(Form("f1_sigmaM2_%s_ch%i",particles[par],ch),"[0] + [1]*x + [2]*x*x", min[par] , max[par] );
	    	//graphSigma[ch][par]->GetYaxis()->SetRangeUser( range_sigma[0][par],range_sigma[1][par]);
	    	graphSigma[ch][par]->Fit(f_sigma[ch][par],"RM");
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



void FitTF1_dEdx_nSigma(TFile *fileIn, TFile *fileWithFitM2, Int_t PID_code, Int_t charge){

	TH2D *histo;
	TH1D *histoProjY;
	TH1D *histoProjX;

	TF1 *f1_meanM2;
	TF1 *f1_sigmaM2;
	TF1 *f1_fit;

	Double_t factor=2.0; // 27GeV Run10, 62.4 GeV Run10
	//Double_t factor=1.0; // 7.7, 11.5, 14.5, 19.6, 27 Run18, 39


	TGraphErrors* gr_Mean;
	TGraphErrors* gr_Sigma;
	
	Double_t mean_nSigma[100];
	Double_t sigma_nSigma[100];
	Double_t pt_x[100]={0.0};
	Double_t err[100]={0.0};

	Double_t min[3]={-2.*factor,-3.*factor,-3.*factor};
	Double_t max[3]={ 2.*factor, 6.*factor, 6.*factor};

	Double_t meanM2=0.0;
	Double_t sigmaM2=0.0;

	fileWithFitM2->cd();
		
	f1_meanM2 = (TF1*)fileWithFitM2 -> Get(Form("f1_meanM2_%s_ch%i",particles[PID_code],charge));
	f1_sigmaM2 = (TF1*)fileWithFitM2 -> Get(Form("f1_sigmaM2_%s_ch%i",particles[PID_code],charge));
	
	for(Int_t pt=0; pt<(int)ptBinRange.size()-1; pt++){
		std::cout<<"good\n";

		meanM2 = f1_meanM2->Eval( (ptBinRange[pt]+ptBinRange[pt+1])/2.0 );
		sigmaM2 = f1_sigmaM2->Eval( (ptBinRange[pt]+ptBinRange[pt+1])/2.0 );
		
		histo = (TH2D*)fileIn -> Get(Form("h2_m2VsnSigma_piKp_%i_pt%i", charge, pt));
		histoProjY = (TH1D*)histo->ProjectionY();
		histoProjX = (TH1D*)histo->ProjectionX(Form("h2_m2VsnSigma_ptBin%i_charge%i_px_%s",pt,charge,particles[PID_code]), histoProjY->FindBin( meanM2 - 0.5*sigmaM2 ), histoProjY->FindBin( meanM2 +  0.5*sigmaM2 ));
		
		f1_fit = fitFun1DGaus(histoProjX, 10., (min[PID_code]+max[PID_code])/2., 0.5, 10, min[PID_code], max[PID_code], 0);
		histoProjX->Write();

		mean_nSigma[pt]=f1_fit->GetParameter(1);
		sigma_nSigma[pt]=f1_fit->GetParameter(2);
		pt_x[pt] = (ptBinRange[pt]+ptBinRange[pt+1])/2.0;

		delete f1_fit;
		delete histo;
		delete histoProjY;
		delete histoProjX;

	}

	gr_Mean = new TGraphErrors((int)ptBinRange.size()-1, pt_x, mean_nSigma, err, err );
	gr_Mean->SetName(Form("gr_%s_MeanNSigma_charge%i",particles[PID_code],charge)); 
	gr_Mean->SetTitle(Form("#mu(n#sigma), %s",partLateX[2*PID_code+charge]));
	gr_Mean->GetXaxis()->SetTitle("#font[42]{p_{T}[GeV/c]}");
	gr_Mean->GetYaxis()->SetTitle(Form("#font[42]{#mu(n#sigma)(%s)}",partLateX[2*PID_code+charge]));

	gr_Sigma = new TGraphErrors((int)ptBinRange.size()-1, pt_x, sigma_nSigma, err, err );
	gr_Sigma->SetName(Form("gr_%s_SigmaNSigma_charge%i",particles[PID_code],charge)); 
	gr_Sigma->SetTitle(Form("#sigma(n#sigma), %s",partLateX[2*PID_code+charge])); 
	gr_Sigma->GetXaxis()->SetTitle("#font[42]{p_{T}[GeV/c]}");
	gr_Sigma->GetYaxis()->SetTitle(Form("#font[42]{#sigma(n#sigma)(%s)}",partLateX[2*PID_code+charge]));

	fileWithFitM2->cd();

	gr_Mean->Write();
	gr_Sigma->Write();

}

void FirstFitThree2x2DGaus(TFile *inFileHisto, TFile *inFileFitFun, Int_t pt){
	
	TH2D *histoPos = (TH2D*)inFileHisto -> Get(Form("h2_m2VsnSigma_piKp_0_pt%i", pt));
	TH2D *histoNeg = (TH2D*)inFileHisto -> Get(Form("h2_m2VsnSigma_piKp_1_pt%i", pt));
	
	FitThree2x2DGaus(inFileFitFun,inFileFitFun,histoPos,pt,0);
	FitThree2x2DGaus(inFileFitFun,inFileFitFun,histoNeg,pt,1);
}

void FitThree2x2DGaus(TFile *inFile, TFile *outFile, TH2D *histo, Int_t pt, Int_t charge){

	gStyle->SetOptFit(kTRUE);

	if(ptBinRange[pt] > 1.5) histo->Rebin2D(2);
	

	Double_t factor=1.0;
	//Double_t factor=1.0;

	Double_t pt_point = (ptBinRange[pt]+ptBinRange[pt+1])/2;
	Double_t parThree2DGaus[15] = {0.};
	Double_t parThree2x2DGaus[30] = {0.};

	Double_t limitPar_meanM2[3]={0.02,0.04,0.1};
	Double_t limitPar_meanNSigma[3]={0.05*factor,0.2*factor,0.3*factor};

	Double_t max_nSigmadEdx = 6.0*factor;
	if(pt<2) max_nSigmadEdx = 10.0*factor;

	TF2 *Three2DGaus = new TF2(Form("Three2DGaus%ich%i",pt,charge),"xygaus(0)+xygaus(5)+xygaus(10)", -3.*factor,max_nSigmadEdx,-0.3,0.95);
	TF2 *Three2x2DGaus = new TF2(Form("Three2x2DGaus%ich%i",pt,charge),"xygaus(0)+xygaus(5)+xygaus(10)+xygaus(15)+xygaus(20)+xygaus(25)", -3*factor,max_nSigmadEdx,-0.4,1.4);

	inFile->cd();

	for(Int_t par=0; par<3; par++){

		Three2DGaus->SetParameter(0+5*par,100);
		Three2DGaus->SetParameter(1+5*par, readParFromGraph(inFile, Form("gr_%s_MeanNSigma_charge%i",particles[par],charge),pt_point));
		Three2DGaus->SetParameter(2+5*par, readParFromGraph(inFile, Form("gr_%s_SigmaNSigma_charge%i",particles[par],charge),pt_point));
		Three2DGaus->SetParameter(3+5*par, readParFromFun(inFile, Form("f1_meanM2_%s_ch%i",particles[par],charge), pt_point) );
		Three2DGaus->SetParameter(4+5*par, readParFromFun(inFile, Form("f1_sigmaM2_%s_ch%i",particles[par],charge), pt_point) );

		Three2DGaus->SetParLimits(1+5*par, readParFromGraph(inFile, Form("gr_%s_MeanNSigma_charge%i",particles[par],charge),pt_point)-limitPar_meanNSigma[par],
										   readParFromGraph(inFile, Form("gr_%s_MeanNSigma_charge%i",particles[par],charge),pt_point)+limitPar_meanNSigma[par] );
		Three2DGaus->SetParLimits(2+5*par, 0.4*readParFromGraph(inFile, Form("gr_%s_SigmaNSigma_charge%i",particles[par],charge),pt_point) ,
										   1.6*readParFromGraph(inFile, Form("gr_%s_SigmaNSigma_charge%i",particles[par],charge),pt_point));
		Three2DGaus->SetParLimits(3+5*par, readParFromFun(inFile, Form("f1_meanM2_%s_ch%i",particles[par],charge), pt_point)-limitPar_meanM2[par] ,
										   readParFromFun(inFile, Form("f1_meanM2_%s_ch%i",particles[par],charge), pt_point)+limitPar_meanM2[par]);
		Three2DGaus->SetParLimits(4+5*par, 0.4*readParFromFun(inFile, Form("f1_sigmaM2_%s_ch%i",particles[par],charge), pt_point) ,
										   1.6*readParFromFun(inFile, Form("f1_sigmaM2_%s_ch%i",particles[par],charge), pt_point) );

	}

	Three2DGaus->GetParameters(&parThree2DGaus[0]);
	std::cout<<"\n Par from 1d fit\n";
	for(Int_t p=0; p<15; p++){
		std::cout<<"\t"<<p<<"\t\t"<<parThree2DGaus[p]<<"\n";
	}
	std::cout<<"\n";

	histo->Fit(Three2DGaus,"RM");
	Three2DGaus->GetParameters(&parThree2DGaus[0]);
	Three2DGaus->SetParameters(parThree2DGaus);
	std::cout<<"\nStep 2\n";
	histo->Fit(Three2DGaus,"RM");
	std::cout<<"\nStep 3\n";
	histo->Fit(Three2DGaus,"RM");


	Three2DGaus->GetParameters(&parThree2x2DGaus[0]);
	Three2DGaus->GetParameters(&parThree2x2DGaus[15]);
	
	Three2x2DGaus->SetParameters(parThree2x2DGaus);

	for(Int_t i=0; i<6; i++){
		Three2x2DGaus->FixParameter(1+5*i,parThree2x2DGaus[1+5*i]);
		Three2x2DGaus->FixParameter(3+5*i,parThree2x2DGaus[3+5*i]);
	}

	for(Int_t i=0; i<3; i++){
		Three2x2DGaus->SetParLimits(0+5*i, 0.8*parThree2x2DGaus[0+5*i], 1.2*parThree2x2DGaus[0+5*i]);
		Three2x2DGaus->SetParLimits(2+5*i, 0.8*parThree2x2DGaus[2+5*i], 1.2*parThree2x2DGaus[2+5*i]);		
		Three2x2DGaus->SetParLimits(4+5*i, 0.8*parThree2x2DGaus[4+5*i], 1.2*parThree2x2DGaus[4+5*i]);
	}

	for(Int_t i=3; i<6; i++){
		Three2x2DGaus->SetParameter(0+5*i, parThree2x2DGaus[0+5*i]/20.);
		Three2x2DGaus->SetParameter(2+5*i, 1.4*parThree2x2DGaus[2+5*i]);
		Three2x2DGaus->SetParameter(4+5*i, 2.5*parThree2x2DGaus[4+5*i]);

		Three2x2DGaus->SetParLimits(0+5*i, 0.0, parThree2x2DGaus[0+5*i]/10);
		Three2x2DGaus->SetParLimits(2+5*i, 0.8*parThree2x2DGaus[2+5*i], 2.0*parThree2x2DGaus[2+5*i]);		
		Three2x2DGaus->SetParLimits(4+5*i, parThree2x2DGaus[4+5*i], 4.0*parThree2x2DGaus[4+5*i]);
	}

	std::cout<<"\nStep 1\n";
	histo->Fit(Three2x2DGaus,"RM");
	std::cout<<"\nStep 2\n";
	histo->Fit(Three2x2DGaus,"RM");

	TCanvas *can = new TCanvas(Form("canvas%i%i",pt,charge),"plot",1024,1024);
	histo->Draw("colz");
	Three2x2DGaus->Draw("same");
	//can->SaveAs(Form("./FitThree2x2DGaus_pt%i_ch%i_Run18.png",pt,charge));

	outFile->cd();

	histo->Write();
	Three2DGaus->Write();
	Three2x2DGaus->Write();


}

Double_t readParFromFun(TFile *inFile, TString name, Double_t pt){

	Double_t vol=0.;
	inFile->cd();
	TF1 *fun = (TF1*)inFile -> Get(name.Data());
	vol = fun->Eval(pt);

	return vol;
}

Double_t readParFromGraph(TFile *inFile, TString name, Double_t pt){
	
	Double_t vol=0.;
	inFile->cd();
	TGraphErrors *graph = (TGraphErrors*)inFile -> Get(name.Data());
	vol = graph->Eval(pt);

	return vol;
}

void DrawAndSavePar(const Char_t *name, TFile *outFile){

	TFile *inFileAfterFirstFit[100];
	TF2 *readFitInFile[100][2];
	TH2D *readHisto2DInFile[100][2];

	TGraphErrors *gr_MeanNS[3][2];
	TGraphErrors *gr_MeanM2[3][2];
	TGraphErrors *gr_SigmaNS[3][2];
	TGraphErrors *gr_SigmaM2[3][2];

	Double_t MeanNS[3][2][100];
	Double_t MeanM2[3][2][100];
	Double_t SigmaNS[3][2][100];
	Double_t SigmaM2[3][2][100];

	Double_t x_pt[100];
	Double_t x_pt_err[100]={0.0};
	Double_t y_err[100]={0.0};

	for(Int_t pti=0; pti<(int)ptBinRange.size()-1; pti++){
		inFileAfterFirstFit[pti] = new TFile(Form("%s/first_27GeV_ptBin_%i.root",name,pti),"READ");
		for(Int_t ch=0; ch<2; ch++){
			readFitInFile[pti][ch]= (TF2*)inFileAfterFirstFit[pti]-> Get(Form("Three2x2DGaus%ich%i",pti,ch));
			readHisto2DInFile[pti][ch]= (TH2D*)inFileAfterFirstFit[pti]-> Get(Form("h2_m2VsnSigma_piKp_%i_pt%i",ch,pti));
		}
	}

	for(Int_t pti=0; pti<(int)ptBinRange.size()-1; pti++){
		for(Int_t par=0; par<3; par++){
			for(Int_t ch=0; ch<2; ch++){

				x_pt[pti]=(ptBinRange[pti]+ptBinRange[pti+1])/2.;

				MeanNS[par][ch][pti] = readFitInFile[pti][ch]->GetParameter(1 + 5*par);
				SigmaNS[par][ch][pti] = readFitInFile[pti][ch]->GetParameter(2 + 5*par);

				MeanM2[par][ch][pti] = readFitInFile[pti][ch]->GetParameter(3 + 5*par);
				SigmaM2[par][ch][pti] = readFitInFile[pti][ch]->GetParameter(4 + 5*par);
			}
		}
	}

	for(Int_t par=0; par<3; par++){
		for(Int_t ch=0; ch<2; ch++){
			gr_MeanNS[par][ch] = new TGraphErrors(ptBinRange.size()-1, x_pt, MeanNS[par][ch], x_pt_err, y_err);
			gr_MeanNS[par][ch]->SetName(Form("gr_2DMeanNS_%s_ch%i",particles[par],ch));

			gr_MeanM2[par][ch] = new TGraphErrors(ptBinRange.size()-1, x_pt, MeanM2[par][ch], x_pt_err, y_err);
			gr_MeanM2[par][ch]->SetName(Form("gr_2DMeanM2_%s_ch%i",particles[par],ch));

			gr_SigmaNS[par][ch] = new TGraphErrors(ptBinRange.size()-1, x_pt, SigmaNS[par][ch], x_pt_err, y_err);
			gr_SigmaNS[par][ch]->SetName(Form("gr_2DSigmaNS_%s_ch%i",particles[par],ch));
			
			gr_SigmaM2[par][ch] = new TGraphErrors(ptBinRange.size()-1, x_pt, SigmaM2[par][ch], x_pt_err, y_err);
			gr_SigmaM2[par][ch]->SetName(Form("gr_2DSigmaM2_%s_ch%i",particles[par],ch));
		
		}
	}

	TF1 *funMean[3][2];
	TF1 *funSigma[3][2];

	for(Int_t par=0; par<3; par++){
		for(Int_t ch=0; ch<2; ch++){
			inFileAfterFirstFit[0]->cd();
			funMean[par][ch] = (TF1*)inFileAfterFirstFit[0]->Get(Form("f1_meanM2_%s_ch%i", particles[par],ch));
			funSigma[par][ch] = (TF1*)inFileAfterFirstFit[0]->Get(Form("f1_sigmaM2_%s_ch%i", particles[par],ch));
		}
	}

	outFile->cd();

	for(Int_t par=0; par<3; par++){
		for(Int_t ch=0; ch<2; ch++){
			funMean[par][ch]->Write();
			funSigma[par][ch]->Write();		
		}
	}

	for(Int_t pti=0; pti<(int)ptBinRange.size()-1; pti++){
		for(Int_t ch=0; ch<2; ch++){
			readFitInFile[pti][ch]->Write();
			readHisto2DInFile[pti][ch]->Write();
		}
	}

	for(Int_t ch=0; ch<2; ch++){
		for(Int_t par=0; par<3; par++){
			gr_MeanNS[par][ch]->Write();
			gr_SigmaNS[par][ch]->Write();
			gr_MeanM2[par][ch]->Write();
			gr_SigmaM2[par][ch]->Write();
		}
	}

	for(Int_t pti=0; pti<(int)ptBinRange.size()-1; pti++){
		inFileAfterFirstFit[pti]->Close();
	}

}


void GetParFunction(TFile *outFile){
	
	TCanvas *can[100];
	
	gStyle->SetOptFit(kTRUE);

	Int_t corr=0;

	TGraphErrors *graph[3][2][4];
	TString txt_charge[] = {"Pos","Neg"};
	TString txt_par[] = {"Pion","Kaon","Proton"};
	TString txt_volum[] = {"MeanNS","SigmaNS","MeanM2","SigmaM2"};
	TString txt_volum2[] = {"#mu(n#sigma_{#pi})","#sigma(n#sigma_{#pi})","#mu(m^2)","#sigma(m^2)"};
	TString txt_title1[] = {"Mean","Width"};
	TString txt_title2[] = {"n#sigma","m^{2}"};
	TString txt_bg[] = {"Peak","BG"};
	Int_t number_can=0;

	TF1 *f_meanNSigma[3][2];
	TF1 *f_meanM2[3][2];
	TF1 *f_wigthM2[3][2];
	TF1 *f_wigthNSigma[3][2];

	for(Int_t ch=0; ch<2; ch++){

		f_meanM2[0][ch]= new TF1(Form("f_meanM2_%s_ch%i",particles[0],ch),"[0] + [1]*x + [2]*x*x", 0.2 , 1.6);
		f_meanM2[1][ch]= new TF1(Form("f_meanM2_%s_ch%i",particles[1],ch),"[0] + [1]*x + [2]*x*x", 0.2 , 1.4);
		f_meanM2[2][ch]= new TF1(Form("f_meanM2_%s_ch%i",particles[2],ch),"[0] + [1]*x + [2]*x*x", 0.5 , 2.5);
		
		f_wigthM2[0][ch]= new TF1(Form("f_wigthM2_%s_ch%i",particles[0],ch),"[0] + [1]*x + [2]*x*x", 0.2 , 1.6);
		f_wigthM2[1][ch]= new TF1(Form("f_wigthM2_%s_ch%i",particles[1],ch),"[0] + [1]*x + [2]*x*x", 0.2 , 1.4);
		f_wigthM2[2][ch]= new TF1(Form("f_wigthM2_%s_ch%i",particles[2],ch),"[0] + [1]*x + [2]*x*x", 0.5 , 2.5);
		
		f_wigthNSigma[0][ch]= new TF1(Form("f_wigthNSigma_%s_ch%i",particles[0],ch),"[0] + [1]*x + [2]*x*x", 0.2 , 1.6);
		f_wigthNSigma[1][ch]= new TF1(Form("f_wigthNSigma_%s_ch%i",particles[1],ch),"[0] + [1]*x + [2]*x*x", 0.2 , 1.4);
		f_wigthNSigma[2][ch]= new TF1(Form("f_wigthNSigma_%s_ch%i",particles[2],ch),"[0] + [1]*x + [2]*x*x + [3]*x*x*x", 0.5 , 2.5);
		
		//f_meanNSigma[0][ch]= new TF1(Form("f_%s",txt_par[0].Data())),"[0]+[1]*x+[2]*x*x + [3]*x*x*x",0.2 , 1.6);
		f_meanNSigma[0][ch]= new TF1(Form("f_meanNSigma_%s_ch%i",particles[0],ch),"[0]*TMath::Power(1.0/x,[1]) + [2]*x + [3]",0.2 , 1.6);
		f_meanNSigma[1][ch]= new TF1(Form("f_meanNSigma_%s_ch%i",particles[1],ch),"[0]*TMath::Power(1.0/x,[1]) + [2]*x + [3]",0.2 , 1.4);
		f_meanNSigma[2][ch]= new TF1(Form("f_meanNSigma_%s_ch%i",particles[2],ch),"[0]*TMath::Power(1.0/x,[1]) + [2]*x + [3]",0.5 , 2.5);

	}

	outFile->cd();


	for(Int_t par=0; par<3; par++){
		for(Int_t ch=0; ch<2; ch++){
			for(Int_t vol=0; vol<4; vol++){

				graph[par][ch][vol] = (TGraphErrors*)outFile->Get(Form("gr_2D%s_%s_ch%i",txt_volum[vol].Data(),particles[par], ch));
          		//graph[par][ch][vol]->SetName(Form("gr_%s%s_%s_%s",txt_par[par].Data(),txt_charge[ch].Data(),txt_volum[vol].Data(),txt_bg.Data()));
          		//graph[par][ch][vol]->SetTitle(Form("%s %s %s",txt_title1[(int)(vol%2)].Data(),txt_title2[(int)(vol/2)].Data(),txt_par[par].Data()));

				can[number_can] = new TCanvas(Form("can_%i",number_can),"Can",640,640);
				can[number_can]->cd();
				graph[par][ch][vol]->Draw("ALP");
				graph[par][ch][vol]->GetXaxis()->SetTitle("p_{t}[GeV/c]");
    			graph[par][ch][vol]->GetYaxis()->SetTitle(Form("%s %s",txt_volum2[vol].Data(), txt_par[par].Data()));
    			graph[par][ch][vol]->GetYaxis()->SetTitleOffset(1.0);

				if(vol==0){
					//graph_star_meanSigma[par]->Draw("L*same");
					graph[par][ch][vol]->Fit(f_meanNSigma[par][ch],"RM");
					graph[par][ch][vol]->Fit(f_meanNSigma[par][ch],"RM");
				}
				if(vol==1){
					f_wigthNSigma[par][ch]->SetLineColor(kBlack);
					graph[par][ch][vol]->Fit(f_wigthNSigma[par][ch],"RM");
					graph[par][ch][vol]->Fit(f_wigthNSigma[par][ch],"RM");
					//graph[par][ch][vol]->Fit(f_wigthNSigma[par][ch],"RM","",0.4,2.6);
				}
				if(vol==2){
					f_meanM2[par][ch]->SetLineColor(kBlack);
					graph[par][ch][vol]->Fit(f_meanM2[par][ch],"RM");
					graph[par][ch][vol]->Fit(f_meanM2[par][ch],"RM");
				}
				if(vol==3){
					f_wigthM2[par][ch]->SetLineColor(kBlack);
					graph[par][ch][vol]->Fit(f_wigthM2[par][ch],"RM");
					graph[par][ch][vol]->Fit(f_wigthM2[par][ch],"RM");
				}
				//graph[par][vol][1]->Draw("L*same");
				//can[number_can]->SaveAs(Form("./new_code_pict/firstFitPar/%s%s_%s_Run18.png",txt_volum[vol].Data(),txt_par[par].Data(),txt_charge[ch].Data()));
				outFile->cd();
				can[number_can]->Write();
				number_can++;
			}
		}
	}

	for(Int_t ch=0; ch<2; ch++){
		for(Int_t par=0; par<3; par++){
			f_meanNSigma[par][ch]->SetRange(0.0,3.2);
			f_meanNSigma[par][ch]->Write();
			
			f_wigthNSigma[par][ch]->SetRange(0.0,3.2);
			f_wigthNSigma[par][ch]->Write();
			
			f_meanM2[par][ch]->SetRange(0.0,3.2);
			f_meanM2[par][ch]->Write();
			
			f_wigthM2[par][ch]->SetRange(0.0,3.2);
			f_wigthM2[par][ch]->Write();

		}
	}

}


void FirstFitThree2x2DGausNewXY(TFile *inFileHisto, TFile *outFile, Int_t pt){
	
	TH2D *histoPos = (TH2D*)inFileHisto -> Get(Form("h2_m2VsnSigma_piK_0_pt%i_new", pt));
	TH2D *histoNeg = (TH2D*)inFileHisto -> Get(Form("h2_m2VsnSigma_piK_1_pt%i_new", pt));
	
	//histoPos->RebinX(2);
	//histoNeg->RebinX(2);

	FitThree2x2DGausNewXY(outFile,histoPos,pt,0);
	FitThree2x2DGausNewXY(outFile,histoNeg,pt,1);
}

void FitThree2x2DGausNewXY(TFile *outFile, TH2D *histo, Int_t pt, Int_t charge){

	gStyle->SetOptFit(kTRUE);

	Double_t pt_point = (ptBinRange[pt]+ptBinRange[pt+1])/2;
	Double_t parThree2DGaus[15] = {0.};
	Double_t parThree2x2DGaus[30] = {0.};

	Double_t limitPar_meanNewX[3]={0.1,0.1,0.2};
	Double_t limitPar_meanNewY[3]={0.1,0.1,0.2};

	TH1D *histoProjX = (TH1D*)histo->ProjectionX(Form("px_ptBin%i_charge%i",pt,charge));

	TF1 *tf1_Gaus;
	TF1 *One1DGaus = new TF1(Form("One1DGaus%ich%i",pt,charge),"gaus", -1.0,0.1);
	TF1 *Three1DGaus = new TF1(Form("Three1DGaus%ich%i",pt,charge),"gaus(0)+gaus(3)+gaus(6)", -0.4,1.5);
	TF1 *Three2x1DGaus = new TF1(Form("Three2x1DGaus%ich%i",pt,charge),"gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)+gaus(15)", -0.4,1.5);
	TF2 *Three2DGaus = new TF2(Form("Three2DGaus%ich%iNewXY",pt,charge),"xygaus(0)+xygaus(5)+xygaus(10)", -0.3,1.2,-1.2,0.5);
	TF2 *Three2x2DGaus = new TF2(Form("Three2x2DGaus%ich%iNewXY",pt,charge),"xygaus(0)+xygaus(5)+xygaus(10)+xygaus(15)+xygaus(20)+xygaus(25)", -0.5,1.5,-1.3,0.6);

	Double_t meanNewX1d[3]={0.0,0.23,0.87};
	Double_t meanNewX1min[3]={-0.2,0.18,0.82};
	Double_t meanNewX1smax[3]={0.15,0.30,1.1};

	for(Int_t par=0; par<3; par++){
		tf1_Gaus = fitFun1DGaus(histoProjX , 10., meanNewX1d[par], 0.01, 10, meanNewX1min[par], meanNewX1smax[par], 0);
		Three1DGaus->SetParameter(0+3*par,tf1_Gaus->GetParameter(0));
		Three1DGaus->SetParameter(1+3*par,tf1_Gaus->GetParameter(1));
		Three1DGaus->SetParameter(2+3*par,tf1_Gaus->GetParameter(2));
		delete tf1_Gaus;
	}
	histoProjX->Fit(Three1DGaus,"RM");
	histoProjX->Fit(Three1DGaus,"RM");

	Double_t meanProtonNewX = Three1DGaus->GetParameter(7);
	Double_t sigmaProtonNewX = TMath::Abs(Three1DGaus->GetParameter(8));

	TH1D *histoProjY = (TH1D*)histo->ProjectionY(Form("py_ptBin%i_charge%i",pt,charge),histoProjX->FindBin(meanProtonNewX - sigmaProtonNewX),histoProjX->FindBin(meanProtonNewX + sigmaProtonNewX));
	histoProjY->Fit(One1DGaus,"RM");

	for(Int_t par=0; par<3; par++){

		Three2DGaus->SetParameter(0+5*par,Three1DGaus->GetParameter(0+3*par));
		Three2DGaus->SetParameter(1+5*par, Three1DGaus->GetParameter(1+3*par));
		Three2DGaus->SetParameter(2+5*par, TMath::Abs(Three1DGaus->GetParameter(2+3*par)));
		if(par==2){
			Three2DGaus->SetParameter(3+5*par, One1DGaus->GetParameter(1));
			Three2DGaus->SetParameter(4+5*par, TMath::Abs(One1DGaus->GetParameter(2)));
		}else{
			Three2DGaus->SetParameter(3+5*par, 0.0);
			Three2DGaus->SetParameter(4+5*par, 1.75*TMath::Abs(Three1DGaus->GetParameter(2+3*par)));
		}

		Three2DGaus->SetParLimits(1+5*par, Three1DGaus->GetParameter(1+3*par)-limitPar_meanNewX[par], Three1DGaus->GetParameter(1+3*par)+limitPar_meanNewX[par]);
		Three2DGaus->SetParLimits(2+5*par, 0.7*TMath::Abs( Three1DGaus->GetParameter(2+3*par)), 1.3*TMath::Abs( Three1DGaus->GetParameter(2+3*par)));
		if(par==2){
			Three2DGaus->SetParLimits(3+5*par, One1DGaus->GetParameter(1)-limitPar_meanNewY[par], One1DGaus->GetParameter(1)+limitPar_meanNewY[par]);
			Three2DGaus->SetParLimits(4+5*par, 0.5*TMath::Abs( One1DGaus->GetParameter(2)), 2.*TMath::Abs( One1DGaus->GetParameter(2)));
		}else{
			Three2DGaus->SetParLimits(3+5*par, 0.0-limitPar_meanNewY[par], 0.0+limitPar_meanNewY[par]);
			Three2DGaus->SetParLimits(4+5*par, 0.5*TMath::Abs( Three1DGaus->GetParameter(2+3*par)), 2*TMath::Abs( Three1DGaus->GetParameter(2+3*par)));
		}	

	}

	Three2DGaus->GetParameters(&parThree2DGaus[0]);
	std::cout<<"\n Par from 1d fit\n";
	for(Int_t p=0; p<15; p++){
		std::cout<<"\t"<<p<<"\t\t"<<parThree2DGaus[p]<<"\n";
	}
	std::cout<<"\n";

	histo->Fit(Three2DGaus,"RM");
	Three2DGaus->GetParameters(&parThree2DGaus[0]);
	Three2DGaus->SetParameters(parThree2DGaus);
	std::cout<<"\nStep 2\n";
	histo->Fit(Three2DGaus,"RM");


	Three2DGaus->GetParameters(&parThree2x2DGaus[0]);
	Three2DGaus->GetParameters(&parThree2x2DGaus[15]);
	
	Three2x2DGaus->SetParameters(parThree2x2DGaus);

	for(Int_t i=0; i<6; i++){
		Three2x2DGaus->FixParameter(1+5*i,parThree2x2DGaus[1+5*i]);
		Three2x2DGaus->FixParameter(3+5*i,parThree2x2DGaus[3+5*i]);
	}

	for(Int_t i=0; i<3; i++){
		Three2x2DGaus->SetParLimits(0+5*i, 0.8*parThree2x2DGaus[0+5*i], 1.2*parThree2x2DGaus[0+5*i]);
		Three2x2DGaus->SetParLimits(2+5*i, 0.8*parThree2x2DGaus[2+5*i], 1.2*parThree2x2DGaus[2+5*i]);		
		Three2x2DGaus->SetParLimits(4+5*i, 0.8*parThree2x2DGaus[4+5*i], 1.2*parThree2x2DGaus[4+5*i]);
	}

	for(Int_t i=3; i<6; i++){
		Three2x2DGaus->SetParameter(0+5*i, parThree2x2DGaus[0+5*i]/20.);
		Three2x2DGaus->SetParameter(2+5*i, 1.4*parThree2x2DGaus[2+5*i]);
		Three2x2DGaus->SetParameter(4+5*i, 2.5*parThree2x2DGaus[4+5*i]);

		Three2x2DGaus->SetParLimits(0+5*i, 0.0, parThree2x2DGaus[0+5*i]/10);
		Three2x2DGaus->SetParLimits(2+5*i, 0.8*parThree2x2DGaus[2+5*i], 2.0*parThree2x2DGaus[2+5*i]);		
		Three2x2DGaus->SetParLimits(4+5*i, parThree2x2DGaus[4+5*i], 4.0*parThree2x2DGaus[4+5*i]);
	}

	histo->Fit(Three2x2DGaus,"RM");
	histo->Fit(Three2x2DGaus,"RM");

	//TCanvas *can = new TCanvas(Form("canvas%i%i",pt,charge),"plot",1024,1024);
	//histo->Draw("colz");
	//Three2x2DGaus->Draw("same");
	//can->SaveAs(Form("./FitThree2x2DGausNewXY_pt%i_ch%i_Run18.png",pt,charge));

	outFile->cd();

	histoProjX->Write();
	histoProjY->Write();
	histo->Write();
	Three2DGaus->Write();
	Three2x2DGaus->Write();

}