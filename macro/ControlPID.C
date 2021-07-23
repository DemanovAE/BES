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
//#include "/mnt/pool/rhic/1/demanov/basov/hpc_scripts/macro/Constants.h"
//#include "./Constants.h"

const Float_t pion_mass = 0.13957061;
const Float_t kaon_mass = 0.493677;
const Float_t proton_mass = 0.9382720813;

const Int_t point = 21;
const Double_t m[point] = {0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.4,2.6,2.8,3.0,3.2};


const Char_t *particles[] = {"Pion","Kaon","Proton"};
const Char_t *partLateX[] = {"#pi^{+}","#pi^{-}","K^{+}","K^{-}","p","#bar{p}"};


//_______________________________________________________________________________________________________________________________________________________________________
void nSigmaDistributionMassSquare(TFile *fileM2, TFile *fileOut, Int_t PID_code);
void DrawParFitMassSquare(TFile *fileM2, TFile *fileOut, Int_t PID_code, Int_t draw_mode);
void DrawM2fit(TFile *fileFitFun, TFile *fileNSigma, TFile *fileOut, Int_t PID_code, Int_t pt_bin);

void FitTF1_dEdx_nSigma(TFile *fileInM2, TFile *fileOut, Int_t PID_code);
void Draw_FitTF1_dEdx_nSigma(TFile *fileInM2, TFile *fileOut, Int_t PID_code, Int_t draw_mode);

void FirstFitTree2x2DGaus(TFile *inFileHisto, TFile *inFileFitFun, TFile *outFile, Int_t pt);
void FitTree2x2DGaus(TFile *inFile, TFile *outFile, TH2D *histo, Int_t pt, Int_t charge);
Double_t readParFromFun(TFile *inFile, TString name, Double_t pt);
Double_t readParFromGraph(TFile *inFile, TString name, Double_t pt);
void WriteGausFit(TFile *inFile, TFile *outFile, Int_t pt, TString mode);
void GetParFit(TFile *outFile);
//void readFunctionPar(TFile *inFile, Int_t PID_code, Int_t charge, Double_t pt, Double_t MeanNSigma, Double_t SigmaNSigma, Double_t MeanM2, Double_t SigmaM2);
void FirstFitTree2x2DGausNewXY(TFile *inFileHisto, TFile *inFileFitFun, TFile *outFile, Int_t pt);
void FitTree2x2DGausNewXY(TFile *inFile, TFile *outFile, TH2D *histo, Int_t pt, Int_t charge);


TF1 *fitFun1DGaus(TH1D *histo, Double_t cons , Double_t mean , Double_t sigma , Double_t min , Double_t max , Int_t flag_fit);
TF1 *fitFunTwo1DGaus(TH1D *histo, Double_t cons1 , Double_t mean1 , Double_t sigma1 , Double_t cons2, Double_t mean2 , Double_t sigma2 , Double_t min , Double_t max );

//_______________________________

void ControlPID(const Char_t *inFileName = "iTest.root", const Char_t *outFileName = "oTest.root", const Int_t pt = 0, const Int_t energy = 39, const Int_t step = 30){

  	TFile *inFile;
  	TFile *inFileFitFun;
  	TFile *outFile;

  	if(step==1){
		
		inFile = new TFile(inFileName, "READ");
  		outFile = new TFile(outFileName, "RECREATE");

  		nSigmaDistributionMassSquare(inFile,outFile,0);
  		nSigmaDistributionMassSquare(inFile,outFile,1);
  		nSigmaDistributionMassSquare(inFile,outFile,2);

  		outFile->Close();

  		outFile = new TFile(outFileName, "UPDATE");

  		DrawParFitMassSquare(outFile,outFile,0,1);
  		DrawParFitMassSquare(outFile,outFile,1,1);
  		DrawParFitMassSquare(outFile,outFile,2,1);

 
  		outFile->Close();

  		outFile = new TFile(outFileName, "UPDATE");

  		FitTF1_dEdx_nSigma(inFile, outFile, 0);
  		FitTF1_dEdx_nSigma(inFile, outFile, 1);
  		FitTF1_dEdx_nSigma(inFile, outFile, 2);

  		inFile->Close();
  		outFile->Close();

  		outFile = new TFile(outFileName, "UPDATE");

  		Draw_FitTF1_dEdx_nSigma(outFile, outFile, 0, 1);
  		Draw_FitTF1_dEdx_nSigma(outFile, outFile, 1, 1);
  		Draw_FitTF1_dEdx_nSigma(outFile, outFile, 2, 1);

  		outFile->Close();

  	}

  	if(step==12){

  		inFile = new TFile(inFileName, "READ");
  		outFile = new TFile(outFileName, "RECREATE");
  		//inFileFitFun = new TFile("./FitFunM2_27GeVRun18.root", "READ");
  		inFileFitFun = new TFile("/mnt/pool/rhic/1/demanov/basov/hpc_scripts/PIDcomb/FitFunM2_27GeVRun18.root", "READ");
		
		FirstFitTree2x2DGaus(inFile, inFileFitFun, outFile, pt);
		outFile->Close();
		inFile->Close();
  	
  	}

  	if(step==13){

  		outFile = new TFile(outFileName, "UPDATE");

 		for(Int_t i=0; i<point-1; i++){
 			inFile = new TFile(Form("%s_ptBin_%i.root",inFileName,i),"READ");
 			WriteGausFit(inFile, outFile, i,"");
 			inFile->Close();
 		}
 
		outFile->Close();  	
  	
  	}

  	if(step==14){

 		inFile = new TFile(inFileName,"READ");

  		GetParFit(inFile);
 
		inFile->Close();  	
  	
  	}

  	if(step==15){

  		inFile = new TFile(inFileName, "READ");
  		outFile = new TFile(outFileName, "RECREATE");
  		//inFileFitFun = new TFile("./FitFunM2_27GeV.root", "READ");
  		inFileFitFun = new TFile("/mnt/pool/rhic/1/demanov/basov/hpc_scripts/PIDcomb/FitFunM2_27GeVRun18.root", "READ");
		
		FirstFitTree2x2DGausNewXY(inFile, inFileFitFun, outFile, pt);
		outFile->Close();
		inFile->Close();
  	
  	}

  	if(step==16){

  		outFile = new TFile(outFileName, "UPDATE");

 		for(Int_t i=0; i<point-1; i++){
 			inFile = new TFile(Form("%s_ptBin_%i.root",inFileName,i),"READ");
 			WriteGausFit(inFile, outFile, i,"NewXY");
 			inFile->Close();
 		}
 
		outFile->Close();  	
  	
  	}

  	if(step==2){
		inFile = new TFile(inFileName, "READ");

  		TFile *fileNSigma = new TFile("./fitGausRun18/FitFunM2_27GeVRun18.root", "READ");
  		for(Int_t pti=0; pti<20; pti++){
			DrawM2fit(inFile, fileNSigma, outFile, 0, pti);
			DrawM2fit(inFile, fileNSigma, outFile, 1, pti);
			DrawM2fit(inFile, fileNSigma, outFile, 2, pti);
  		}
  		fileNSigma->Close();
  		inFile->Close();
  		outFile->Close();
  	}

  	

}


TF1 *fitFun1DGaus(TH1D *histo, Double_t cons, Double_t mean, Double_t sigma, Double_t min, Double_t max, Int_t flag_fit){
	
	TF1 *gaus = new TF1("","gaus",min,max);
	gaus->SetParameter(0,cons);
	gaus->SetParameter(1,mean);
	gaus->SetParameter(2,sigma);
	if(flag_fit==0){
		histo->Fit(gaus,"RM");
	}
	return gaus;

}

TF1 *fitFunTwo1DGaus(TH1D *histo, Double_t cons1, Double_t mean1, Double_t sigma1, Double_t cons2, Double_t mean2, Double_t sigma2, Double_t min, Double_t max ){

	TF1 *gaus = new TF1("","gaus(0)+gaus(3)",min,max);
	gaus->SetParameter(0,cons1);
	gaus->SetParameter(1,mean1);
	gaus->SetParameter(2,sigma1);
	gaus->SetParameter(3,cons2);
	gaus->SetParameter(4,mean2);
	gaus->SetParameter(5,sigma2);
	histo->Fit(gaus,"RM");
	return gaus;

}

void nSigmaDistributionMassSquare(TFile *fileM2, TFile *fileOut, Int_t PID_code){

	std::vector<Double_t> ptBinRange;
	for(Int_t i=0; i<point; i++){
		ptBinRange.push_back(m[i]);
	}

	gStyle->SetOptStat("n");
    gStyle->SetOptFit(1);

    TF1 *tf1_fitFun[2][point];
    TH1D *h1_m2[2][point];

    Double_t mean[2][point]={0.};
    Double_t sigma[2][point]={0.};
    Double_t pt[point] = {0.};
    Double_t err[point] = {0.};

    Double_t min[4]={-0.5, 0.15, 0.7, -0.5};
    Double_t max[4]={0.1, 0.4, 1.2, 0.5};
    Double_t massSqrt[3]={pion_mass*pion_mass, kaon_mass*kaon_mass, proton_mass*proton_mass};

    TH2D *h2_m2VsPt[2]; 
    h2_m2VsPt[0] = (TH2D*)fileM2->Get(Form("h2_m2VsPt_all_ch0"));
    h2_m2VsPt[1] = (TH2D*)fileM2->Get(Form("h2_m2VsPt_all_ch1"));

    TH1D *h1_m2_projY = (TH1D*)h2_m2VsPt[0]->ProjectionX();

    fileOut->cd();

    if(PID_code==0){
	    h2_m2VsPt[0]->Write();
    	h2_m2VsPt[1]->Write();
    }

    for(Int_t ch=0; ch<2; ch++){
        for(Int_t pti=0; pti<(int)ptBinRange.size()-1; pti++){
	        h1_m2[ch][pti] = (TH1D*)h2_m2VsPt[ch]->ProjectionY(Form("h1_m2_%s_charge%i_pt%i", particles[PID_code],ch,pti),h1_m2_projY->FindBin(ptBinRange[pti]),h1_m2_projY->FindBin(ptBinRange[pti+1]));
    		h1_m2[ch][pti]->SetTitle(Form("m^{2}, %.1f < p_{T} < %.1f [GeV/c]",ptBinRange[pti],ptBinRange[pti+1]));
    		if(pti>2)h1_m2[ch][pti]->Rebin(2);
    	}
    }

    Int_t flag=0;
    TF1 *tf1_help[2];

    fileOut->cd();
	for(Int_t ch=0; ch<2; ch++){
	    for(Int_t pti=0; pti<(int)ptBinRange.size()-1; pti++){

	    	if( (PID_code==0 || PID_code==1) && ptBinRange[pti]>=1.5 ){
	    		flag=1;
	    	}else{
	    		flag=0;
	    	}
	    	
	    	if(flag==0){
	    		tf1_fitFun[ch][pti] = fitFun1DGaus(h1_m2[ch][pti], 10., massSqrt[PID_code], 0.01, min[PID_code], max[PID_code],0);
	    	}
	    	if(flag==1){
	    		tf1_help[0] = fitFun1DGaus(h1_m2[ch][pti], 10., massSqrt[0], 0.01, -0.4, 0.1, 0);
	    		tf1_help[1] = fitFun1DGaus(h1_m2[ch][pti], 10., massSqrt[1], 0.01,  0.18, 0.3, 0);
	    		
	    		tf1_fitFun[ch][pti] = fitFunTwo1DGaus(h1_m2[ch][pti], tf1_help[0]->GetParameter(0), tf1_help[0]->GetParameter(1), tf1_help[0]->GetParameter(2),
	    															  tf1_help[1]->GetParameter(0), tf1_help[1]->GetParameter(1), tf1_help[1]->GetParameter(2), min[3], max[3]);

	    		tf1_fitFun[ch][pti] = fitFun1DGaus(h1_m2[ch][pti], tf1_fitFun[ch][pti]->GetParameter(3*PID_code + 0), tf1_fitFun[ch][pti]->GetParameter(3*PID_code + 1), tf1_fitFun[ch][pti]->GetParameter(3*PID_code + 2), min[PID_code], max[PID_code],1);
			}

			tf1_fitFun[ch][pti]->SetName(Form("tf1_%s_m2_charge%i_pt%i",particles[PID_code],ch, pti)); 
			tf1_fitFun[ch][pti]->SetTitle(Form("1DGaus, %s, %.1f < p_{T} < %.1f [GeV/c]",partLateX[2*PID_code+ch],ptBinRange[pti],ptBinRange[pti+1])); 
			
			mean[ch][pti] = tf1_fitFun[ch][pti]->GetParameter(1);
			sigma[ch][pti] = tf1_fitFun[ch][pti]->GetParameter(2);
			
			if(ch==0){
				pt[pti] = (ptBinRange[pti] + ptBinRange[pti+1]) / 2.0;
			}
	    }
	}

    TGraphErrors *graphMean[2];
    TGraphErrors *graphSigma[2];
	for(Int_t ch=0; ch<2; ch++){
		graphMean[ch] = new TGraphErrors((int)ptBinRange.size()-1, pt, mean[ch], err, err);
		graphMean[ch]->SetName(Form("gr_Mean%s_ch%i",particles[PID_code],ch));
		graphMean[ch]->SetTitle(Form("Mean m^{2} (%s) vs p_{T} bin",partLateX[2*PID_code+ch]));
		graphMean[ch]->GetXaxis()->SetTitle("#font[42]{p_{T}[GeV/c]}");
		graphMean[ch]->GetYaxis()->SetTitle(Form("#font[42]{#mu(m^{2})(%s)}",partLateX[2*PID_code+ch]));

		graphMean[ch]->Write();

		graphSigma[ch] = new TGraphErrors((int)ptBinRange.size()-1, pt, sigma[ch], err, err);
		graphSigma[ch]->SetName(Form("gr_Sigma%s_ch%i",particles[PID_code],ch));
		graphSigma[ch]->SetTitle(Form("Sigma m^{2} (%s) vs p_{T} bin",partLateX[2*PID_code+ch]));
		graphSigma[ch]->GetXaxis()->SetTitle("#font[42]{p_{T}[GeV/c]}");
		graphSigma[ch]->GetYaxis()->SetTitle(Form("#font[42]{#sigma(m^{2})(%s)}",partLateX[2*PID_code+ch]));

		graphSigma[ch]->Write();

	}

	for(Int_t ch=0; ch<2; ch++){
	    for(Int_t pti=0; pti<(int)ptBinRange.size()-1; pti++){
			h1_m2[ch][pti]->Write(); 
	    }
	}
	for(Int_t ch=0; ch<2; ch++){
	    for(Int_t pti=0; pti<(int)ptBinRange.size()-1; pti++){
			tf1_fitFun[ch][pti]->Write(); 
	    }
	}

}

void DrawParFitMassSquare(TFile *fileM2, TFile *fileOut, Int_t PID_code, Int_t draw_mode){

	std::vector<Double_t> ptBinRange;
	for(Int_t i=0; i<point; i++){
		ptBinRange.push_back(m[i]);
	}

	gStyle->SetOptStat("n");
    gStyle->SetOptFit(1);

    TGraphErrors *gr_Mean[2]; 
    TGraphErrors *gr_Sigma[2];
    
    for(Int_t ch=0;ch<2;ch++){
		gr_Mean[ch] = (TGraphErrors*)fileM2->Get(Form("gr_Mean%s_ch%i",particles[PID_code],ch));
    	gr_Sigma[ch] = (TGraphErrors*)fileM2->Get(Form("gr_Sigma%s_ch%i",particles[PID_code],ch));
    }

	for(Int_t ch=0; ch<2; ch++){

		gr_Mean[ch]->GetYaxis()->SetTitleOffset(1.5);
		gr_Mean[ch]->SetMarkerStyle(24);
		gr_Mean[ch]->SetMarkerSize(2);
		gr_Mean[ch]->SetMarkerColor(kBlue);
		gr_Mean[ch]->SetLineWidth(2);
		gr_Mean[ch]->SetLineColor(kBlue);

		gr_Sigma[ch]->GetYaxis()->SetTitleOffset(1.5);
		gr_Sigma[ch]->SetMarkerStyle(24);
		gr_Sigma[ch]->SetMarkerSize(2);
		gr_Sigma[ch]->SetMarkerColor(kBlue);
		gr_Sigma[ch]->SetLineWidth(2);
		gr_Sigma[ch]->SetLineColor(kBlue);

	}

	TF1 *f_mean[2];
	TF1 *f_sigma[2];
	TF1 *f_Draw_mean[2];
	TF1 *f_Draw_sigma[2];
	Double_t min[4]={0.2, 0.2, 0.5};
    Double_t max[4]={2.0, 2.0, 2.5};
    
    Double_t range_mean[2][3]={{-0.051, 0.0, 0.82},
    						   {0.04, 0.31, 0.892}};
    Double_t range_sigma[2][3]={{0.0, 0.0, 0.0},
    						    {0.241, 0.31, 0.31}};

    for(Int_t ch=0; ch<2; ch++){

    	if(PID_code==0 || PID_code==1){
			f_Draw_mean[ch] = new TF1(Form("f_Draw_%s_ch%i",particles[PID_code],ch),"[0] + [1]*x + [2]*x*x", 0.2 , 3.2 );
	    	f_mean[ch] = new TF1(Form("f1_meanM2_%s_ch%i",particles[PID_code],ch),"[0] + [1]*x + [2]*x*x", min[PID_code] , max[PID_code] );
	    	gr_Mean[ch]->GetYaxis()->SetRangeUser( range_mean[0][PID_code],range_mean[1][PID_code]);
			gr_Mean[ch]->Fit(f_mean[ch],"RM");
			f_Draw_mean[ch]->SetParameter(0,f_mean[ch]->GetParameter(0));
			f_Draw_mean[ch]->SetParameter(1,f_mean[ch]->GetParameter(1));
			f_Draw_mean[ch]->SetParameter(2,f_mean[ch]->GetParameter(2));
    	}

    	if(PID_code==2){
			f_Draw_mean[ch] = new TF1(Form("f_Draw_%s_ch%i",particles[PID_code],ch),"[0] + [1]*x + [2]*x*x + [3]*x*x*x", 0.2 , 3.2 );
	    	f_mean[ch] = new TF1(Form("f1_meanM2_%s_ch%i",particles[PID_code],ch),"[0] + [1]*x + [2]*x*x + [3]*x*x*x", min[PID_code] , max[PID_code] );
	    	gr_Mean[ch]->GetYaxis()->SetRangeUser( range_mean[0][PID_code],range_mean[1][PID_code]);
			gr_Mean[ch]->Fit(f_mean[ch],"RM");
			f_Draw_mean[ch]->SetParameter(0,f_mean[ch]->GetParameter(0));
			f_Draw_mean[ch]->SetParameter(1,f_mean[ch]->GetParameter(1));
			f_Draw_mean[ch]->SetParameter(2,f_mean[ch]->GetParameter(2));
			f_Draw_mean[ch]->SetParameter(3,f_mean[ch]->GetParameter(3));
    	}
		
		
		f_Draw_sigma[ch] = new TF1(Form("f_Draw_%s_ch%i",particles[PID_code],ch),"[0] + [1]*x + [2]*x*x", 0.2 , 3.2 );
		f_sigma[ch] = new TF1(Form("f1_sigmaM2_%s_ch%i",particles[PID_code],ch),"[0] + [1]*x + [2]*x*x", min[PID_code] , max[PID_code] );
    	gr_Sigma[ch]->GetYaxis()->SetRangeUser( range_sigma[0][PID_code],range_sigma[1][PID_code]);
    	gr_Sigma[ch]->Fit(f_sigma[ch],"RM");
		f_Draw_sigma[ch]->SetParameter(0,f_sigma[ch]->GetParameter(0));
		f_Draw_sigma[ch]->SetParameter(1,f_sigma[ch]->GetParameter(1));
		f_Draw_sigma[ch]->SetParameter(2,f_sigma[ch]->GetParameter(2));

    }

    fileOut->cd();
    for(Int_t ch=0; ch<2;ch++){
    	f_mean[ch]->Write();
    	f_sigma[ch]->Write();
    }

    TCanvas *canvas;
    TPad *pad[4];

    if(draw_mode==0){
	    canvas = new TCanvas(Form("canvas_%s",particles[PID_code]),Form("plot_a"),1200,1200);
	    canvas->GetFrame()->SetFillColor(21);
	    canvas->GetFrame()->SetBorderSize(115);
	    canvas->cd();

	    pad[0] = new TPad(Form("pad_%i_0",PID_code), Form("p1"), 0.0, 0.5, 0.5, 1.0, 0, 0, 0);
	    pad[1] = new TPad(Form("pad_%i_1",PID_code), Form("p2"), 0.5, 0.5, 1.0, 1.0, 0, 0, 0);
	    pad[2] = new TPad(Form("pad_%i_2",PID_code), Form("p3"), 0.0, 0.0, 0.5, 0.5, 0, 0, 0);
	    pad[3] = new TPad(Form("pad_%i_3",PID_code), Form("p4"), 0.5, 0.0, 1.0, 0.5, 0, 0, 0);

	    for(Int_t i=0; i<4; i++){
	        pad[i]->Draw();
	    }

		TLegend *legend0 = new TLegend(0.1,0.1,0.5,0.25);
	    //legend0->SetHeader("dgdjgjdj");
	    legend0->SetHeader(Form("#font[42]{Range of fit: %.1f<p_{T}<%.1f}",min[PID_code],max[PID_code]));
	    legend0->AddEntry( f_Draw_mean[0], "#font[42]{fit}","l");

	    pad[0]->cd();
	    gr_Mean[0]->Draw("APL");
	   	f_Draw_mean[0]->Draw("same");
	    legend0->Draw("same");


	    //h1_fit[0]->SetTitle(Form("%s Pos, %.1f<p_{T}<%.1f", particles[PID_code],ptBinRange[pt_bin], ptBinRange[pt_bin+1] ));
	    
	    pad[1]->cd();
	  	gr_Sigma[0]->Draw("APL");
		f_Draw_sigma[0]->Draw("same");

	    
	    pad[2]->cd();
	   	gr_Mean[1]->Draw("APL");
	   	f_Draw_mean[1]->Draw("same");
	    //h1_fit[0]->SetTitle(Form("%s Pos, %.1f<p_{T}<%.1f", particles[PID_code],ptBinRange[pt_bin], ptBinRange[pt_bin+1] ));
	    
	    pad[3]->cd();
	  	gr_Sigma[1]->Draw("APL");
	  	f_Draw_sigma[1]->Draw("same");
	    
	    canvas->SaveAs(Form("./pict_proton/ParFitM2_%s_Run18.png",particles[PID_code]));
	}


}

void DrawM2fit(TFile *fileFitFun, TFile *fileNSigma, TFile *fileOut, Int_t PID_code, Int_t pt_bin){

	std::vector<Double_t> ptBinRange;
	for(Int_t i=0; i<point; i++){
		ptBinRange.push_back(m[i]);
	}

    TCanvas *canvas = new TCanvas(Form("canvas_%s_%i",particles[PID_code],pt_bin),Form("plot_%i",pt_bin),1200,1200);
    canvas->GetFrame()->SetFillColor(21);
    canvas->GetFrame()->SetBorderSize(115);
    canvas->cd();

    gStyle->SetOptStat("n");
    gStyle->SetOptFit(1);

    TPad *pad[4];

    pad[0] = new TPad(Form("pad_%i_0_%i",PID_code,pt_bin), Form("p1_%i",pt_bin), 0.0, 0.5, 0.5, 1.0, 0, 0, 0);
    pad[1] = new TPad(Form("pad_%i_1_%i",PID_code,pt_bin), Form("p2_%i",pt_bin), 0.5, 0.5, 1.0, 1.0, 0, 0, 0);
    pad[2] = new TPad(Form("pad_%i_2_%i",PID_code,pt_bin), Form("p3_%i",pt_bin), 0.0, 0.0, 0.5, 0.5, 0, 0, 0);
    pad[3] = new TPad(Form("pad_%i_3_%i",PID_code,pt_bin), Form("p4_%i",pt_bin), 0.5, 0.0, 1.0, 0.5, 0, 0, 0);

    for(Int_t i=0; i<4; i++){
        pad[i]->Draw();
    }

    TH1D *h1_fit[2];
    h1_fit[0] = (TH1D*)fileFitFun->Get(Form("h1_m2_%s_charge%i_pt%i", particles[PID_code],0,pt_bin));
    h1_fit[1] = (TH1D*)fileFitFun->Get(Form("h1_m2_%s_charge%i_pt%i", particles[PID_code],1,pt_bin));
    
    TH1D *h1_new[2];
    h1_new[0] = (TH1D*)fileNSigma->Get(Form("h2_m2VsnSigma_ptBin%i_charge%i_px_%s",pt_bin, 0, particles[PID_code]));
    //h1_new[0] = (TH1D*)fileNSigma->Get(Form("h1_m2_%s_charge%i_pt%i",particles[PID_code],0,pt_bin));
    h1_new[1] = (TH1D*)fileNSigma->Get(Form("h2_m2VsnSigma_ptBin%i_charge%i_px_%s",pt_bin, 1, particles[PID_code]));
    //h1_new[1] = (TH1D*)fileNSigma->Get(Form("h1_m2_%s_charge%i_pt%i",particles[PID_code],1,pt_bin));

    h1_new[0]->Rebin(2);
    h1_new[1]->Rebin(2);

    Double_t min=0.;
    Double_t max=0.;
    Int_t flag=0;
    
    if(ptBinRange[pt_bin]<1.3){
        if(PID_code==0){
        	min=-4.0;
        	max= 4.0;	
        }
        if(PID_code==1){
        	min= -3.0;
        	max= 3.0;	
        }
        if(PID_code==2){
        	min= -4.0;
        	max= 4.0;	
        }
    }
    if(ptBinRange[pt_bin]>=1.3){
        if(PID_code==0){
        	min=-4.0;
        	max= 2.0;	
        }
        if(PID_code==1){
        	min= -1.5;
        	max= 3.0;	
        }
        if(PID_code==2){
        	min= -2.5;
        	max= 4.0;	
        }
    }
    if(ptBinRange[pt_bin]>1.55){
        if(PID_code==0){
        	min=-3.0;
        	max= 1.0;	
        }
        if(PID_code==1){
        	min= -1.0;
        	max= 3.0;	
        }
        if(PID_code==2){
        	min= -1.5;
        	max= 3.0;	
        }
        
    }

    std::cout<<"good\n";
    TF1 *gaus1D[2];

    for(Int_t ch=0; ch<2; ch++){
   		gaus1D[ch] = new TF1(Form("gaus1D%spt_ch%i_pt%i", particles[PID_code], ch, pt_bin),"gaus", min, max);
    	gaus1D[ch]->SetLineWidth(2);
    	gaus1D[ch]->SetLineStyle(1);
    	gaus1D[ch]->SetLineColor(kRed);
    }

    pad[0]->cd();
    pad[0]->SetLogy();
    h1_fit[0]->SetLabelSize(0.05,"X");
    h1_fit[0]->SetLabelSize(0.05,"Y");
    h1_fit[0]->SetTitleSize(0.05,"X");
    h1_fit[0]->SetTitleSize(0.05,"Y");
    h1_fit[0]->SetTitleOffset(0.8,"X");
    h1_fit[0]->SetTitleOffset(0.8,"Y");
    //h1_fit[0]->SetStats(0);
    h1_fit[0]->SetTitle(Form("%s Pos, %.1f<p_{T}<%.1f", particles[PID_code],ptBinRange[pt_bin], ptBinRange[pt_bin+1] ));
    h1_fit[0]->Draw();
    
    std::cout<<"good\n";

    pad[1]->cd();
    h1_new[0]->SetLabelSize(0.05,"X");
    h1_new[0]->SetLabelSize(0.05,"Y");
    h1_new[0]->SetTitleSize(0.05,"X");
    h1_new[0]->SetTitleSize(0.05,"Y");
    h1_new[0]->SetTitleOffset(0.8,"X");
    h1_new[0]->SetTitleOffset(0.8,"Y");
    //h1_new[0]->SetStats(0);
    h1_new[0]->Draw();
    h1_new[0]->Fit(gaus1D[0],"RM");

    
    pad[2]->cd();
    pad[2]->SetLogy();
    h1_fit[1]->SetLabelSize(0.05,"X");
    h1_fit[1]->SetLabelSize(0.05,"Y");
    h1_fit[1]->SetTitleSize(0.05,"X");
    h1_fit[1]->SetTitleSize(0.05,"Y");
    h1_fit[1]->SetTitleOffset(0.8,"X");
    h1_fit[1]->SetTitleOffset(0.8,"Y");
    h1_fit[1]->SetTitle(Form("%s Neg, %.1f<p_{T}<%.1f", particles[PID_code],ptBinRange[pt_bin], ptBinRange[pt_bin+1] ));
    //h1_fit[1]->SetStats(0);
    h1_fit[1]->Draw();

    pad[3]->cd();
    if(pt_bin==0 && PID_code==2){
        pad[3]->SetLogy();
    }
    h1_new[1]->SetLabelSize(0.05,"X");
    h1_new[1]->SetLabelSize(0.05,"Y");
    h1_new[1]->SetTitleSize(0.05,"X");
    h1_new[1]->SetTitleSize(0.05,"Y");
    h1_new[1]->SetTitleOffset(0.8,"X");
    h1_new[1]->SetTitleOffset(0.8,"Y");
    //h1_new[1]->SetStats(0);
    h1_new[1]->Draw();
    h1_new[1]->Fit(gaus1D[1],"RM");
    
    canvas->SaveAs(Form("./pict_proton/nSigmaDistr_%s_ptBin%i_Run18.png",particles[PID_code],pt_bin));
    //canvas->SaveAs(Form("./pict_proton/FIT_ptBin%i.pdf",pt_bin));
    
}


void FitTF1_dEdx_nSigma(TFile *fileInM2, TFile *fileOut, Int_t PID_code){

	std::vector<Double_t> ptBinRange;
	for(Int_t i=0; i<point; i++){
		ptBinRange.push_back(m[i]);
	}

	TH2D *histo;
	TH1D *histoProjY;
	TH1D *histoProjX;

	TF1 *f1_meanM2;
	TF1 *f1_sigmaM2;
	TF1 *f1_fit;

	Double_t factor=2.0;

	TGraphErrors* gr_Mean[2];
	TGraphErrors* gr_Sigma[2];
	
	Double_t mean_nSigma[2][point-1];
	Double_t sigma_nSigma[2][point-1];
	Double_t pt_x[2][point-1]={0.0};
	Double_t err[point-1]={0.0};

	Double_t min[3]={-2.*factor,-3.*factor,-3.*factor};
	Double_t max[3]={ 2.*factor, 6.*factor, 6.*factor};

	Double_t meanM2=0.0;
	Double_t sigmaM2=0.0;

	fileOut->cd();

	for(Int_t ch=0; ch<2; ch++){
		
		f1_meanM2 = (TF1*)fileOut -> Get(Form("f1_meanM2_%s_ch%i",particles[PID_code],ch));
		f1_sigmaM2 = (TF1*)fileOut -> Get(Form("f1_sigmaM2_%s_ch%i",particles[PID_code],ch));
		
		for(Int_t pt=0; pt<(int)ptBinRange.size()-1; pt++){
			std::cout<<"good\n";

			meanM2 = f1_meanM2->Eval( (ptBinRange[pt]+ptBinRange[pt+1])/2.0 );
			sigmaM2 = 0.5 * f1_sigmaM2->Eval( (ptBinRange[pt]+ptBinRange[pt+1])/2.0 );
			
			histo = (TH2D*)fileInM2 -> Get(Form("h2_m2VsnSigma_piKp_%i_pt%i", ch, pt));
			histoProjY = (TH1D*)histo->ProjectionY();
			histoProjX = (TH1D*)histo->ProjectionX(Form("h2_m2VsnSigma_ptBin%i_charge%i_px_%s",pt,ch,particles[PID_code]), histoProjY->FindBin( meanM2 - sigmaM2 ), histoProjY->FindBin( meanM2 + sigmaM2 ));
			
			f1_fit = fitFun1DGaus(histoProjX, 10., (min[PID_code]+max[PID_code])/2., 0.5, min[PID_code], max[PID_code], 0);
			histoProjX->Write();

			mean_nSigma[ch][pt]=f1_fit->GetParameter(1);
			sigma_nSigma[ch][pt]=f1_fit->GetParameter(2);
			pt_x[ch][pt] = (ptBinRange[pt]+ptBinRange[pt+1])/2.0;

			delete f1_fit;
			delete histo;
			delete histoProjY;
			delete histoProjX;

		}

		gr_Mean[ch] = new TGraphErrors(point-1, pt_x[ch], mean_nSigma[ch], err, err );
		gr_Mean[ch]->SetName(Form("gr_%s_MeanNSigma_charge%i",particles[PID_code],ch)); 
		gr_Mean[ch]->SetTitle(Form("#mu(n#sigma), %s",partLateX[2*PID_code+ch]));
		gr_Mean[ch]->GetXaxis()->SetTitle("#font[42]{p_{T}[GeV/c]}");
		gr_Mean[ch]->GetYaxis()->SetTitle(Form("#font[42]{#mu(n#sigma)(%s)}",partLateX[2*PID_code+ch]));

		gr_Sigma[ch] = new TGraphErrors(point-1, pt_x[ch], sigma_nSigma[ch], err, err );
		gr_Sigma[ch]->SetName(Form("gr_%s_SigmaNSigma_charge%i",particles[PID_code],ch)); 
		gr_Sigma[ch]->SetTitle(Form("#sigma(n#sigma), %s",partLateX[2*PID_code+ch])); 
		gr_Sigma[ch]->GetXaxis()->SetTitle("#font[42]{p_{T}[GeV/c]}");
		gr_Sigma[ch]->GetYaxis()->SetTitle(Form("#font[42]{#sigma(n#sigma)(%s)}",partLateX[2*PID_code+ch]));


	}

	fileOut->cd();

	for(Int_t ch=0; ch<2; ch++){
		gr_Mean[ch]->Write();
		gr_Sigma[ch]->Write();
	}


}

void Draw_FitTF1_dEdx_nSigma(TFile *fileInM2, TFile *fileOut, Int_t PID_code, Int_t draw_mode){

	std::vector<Double_t> ptBinRange;
	for(Int_t i=0; i<point; i++){
		ptBinRange.push_back(m[i]);
	}

	gStyle->SetOptStat("n");
    gStyle->SetOptFit(1);

    TGraphErrors *gr_Mean[2]; 
    TGraphErrors *gr_Sigma[2];
    
    for(Int_t ch=0;ch<2;ch++){
		gr_Mean[ch] = (TGraphErrors*)fileInM2->Get(Form("gr_%s_MeanNSigma_charge%i",particles[PID_code],ch));
    	gr_Sigma[ch] = (TGraphErrors*)fileInM2->Get(Form("gr_%s_SigmaNSigma_charge%i",particles[PID_code],ch));
    }

	for(Int_t ch=0; ch<2; ch++){

		gr_Mean[ch]->GetYaxis()->SetTitleOffset(1.5);
		gr_Mean[ch]->SetMarkerStyle(24);
		gr_Mean[ch]->SetMarkerSize(2);
		gr_Mean[ch]->SetMarkerColor(kBlue);
		gr_Mean[ch]->SetLineWidth(2);
		gr_Mean[ch]->SetLineColor(kBlue);

		gr_Sigma[ch]->GetYaxis()->SetTitleOffset(1.5);
		gr_Sigma[ch]->SetMarkerStyle(24);
		gr_Sigma[ch]->SetMarkerSize(2);
		gr_Sigma[ch]->SetMarkerColor(kBlue);
		gr_Sigma[ch]->SetLineWidth(2);
		gr_Sigma[ch]->SetLineColor(kBlue);

	}

	TF1 *f_mean[2];
	TF1 *f_sigma[2];
	TF1 *f_Draw_mean[2];
	TF1 *f_Draw_sigma[2];
	Double_t factor=2.0;
	Double_t min[4]={0.4, 0.4, 0.5};
    Double_t max[4]={1.5, 1.5, 2.5};
    
    Double_t range_mean[2][3]={{-0.15, -2.5, -3.0},
    						   {0.13, 6.0, 7.0}};
    Double_t range_sigma[2][3]={{0.45, 0.48, 0.47},
    						    {0.65, 1.20, 1.3}};

    for(Int_t ch=0; ch<2; ch++){

		f_Draw_mean[ch] = new TF1(Form("f_Draw_%s_ch%i",particles[PID_code],ch),"[0]*TMath::Power(1.0/x,[1]) + [2]*x + [3]", 0.2 , 3.2 );
		f_mean[ch] = new TF1(Form("f1_meanNSigma_%s_ch%i",particles[PID_code],ch),"[0]*TMath::Power(1.0/x,[1]) + [2]*x + [3]", min[PID_code] , max[PID_code] );
		gr_Mean[ch]->GetYaxis()->SetRangeUser( range_mean[0][PID_code],range_mean[1][PID_code]);
		gr_Mean[ch]->Fit(f_mean[ch],"RM");
		f_Draw_mean[ch]->SetParameter(0,f_mean[ch]->GetParameter(0));
		f_Draw_mean[ch]->SetParameter(1,f_mean[ch]->GetParameter(1));
		f_Draw_mean[ch]->SetParameter(2,f_mean[ch]->GetParameter(2));
		f_Draw_mean[ch]->SetParameter(3,f_mean[ch]->GetParameter(3));

		f_Draw_sigma[ch] = new TF1(Form("f_Draw_%s_ch%i",particles[PID_code],ch),"[0]*TMath::Power(1.0/x,[1]) + [2]*x + [3]", 0.2 , 3.2 );
		f_sigma[ch] = new TF1(Form("f1_sigmaNSigma_%s_ch%i",particles[PID_code],ch),"[0]*TMath::Power(1.0/x,[1]) + [2]*x + [3]", min[PID_code] , max[PID_code] );
    	gr_Sigma[ch]->GetYaxis()->SetRangeUser( range_sigma[0][PID_code],range_sigma[1][PID_code]);
    	gr_Sigma[ch]->Fit(f_sigma[ch],"RM");
		f_Draw_sigma[ch]->SetParameter(0,f_sigma[ch]->GetParameter(0));
		f_Draw_sigma[ch]->SetParameter(1,f_sigma[ch]->GetParameter(1));
		f_Draw_sigma[ch]->SetParameter(2,f_sigma[ch]->GetParameter(2));
		f_Draw_sigma[ch]->SetParameter(3,f_sigma[ch]->GetParameter(3));
		//f_Draw_sigma[ch]->SetParameter(4,f_sigma[ch]->GetParameter(4));

    }

    fileOut->cd();
    for(Int_t ch=0; ch<2;ch++){
    	f_mean[ch]->Write();
    	f_sigma[ch]->Write();
    }

    TCanvas *canvas;
    TPad *pad[4];

    if(draw_mode==0){
	    canvas = new TCanvas(Form("canvas_%s",particles[PID_code]),Form("plot_a"),1200,1200);
	    canvas->GetFrame()->SetFillColor(21);
	    canvas->GetFrame()->SetBorderSize(115);
	    canvas->cd();

	    pad[0] = new TPad(Form("pad_%i_0",PID_code), Form("p1"), 0.0, 0.5, 0.5, 1.0, 0, 0, 0);
	    pad[1] = new TPad(Form("pad_%i_1",PID_code), Form("p2"), 0.5, 0.5, 1.0, 1.0, 0, 0, 0);
	    pad[2] = new TPad(Form("pad_%i_2",PID_code), Form("p3"), 0.0, 0.0, 0.5, 0.5, 0, 0, 0);
	    pad[3] = new TPad(Form("pad_%i_3",PID_code), Form("p4"), 0.5, 0.0, 1.0, 0.5, 0, 0, 0);

	    for(Int_t i=0; i<4; i++){
	        pad[i]->Draw();
	    }

		TLegend *legend0 = new TLegend(0.1,0.1,0.7,0.25);
	    //legend0->SetHeader("dgdjgjdj");
	    legend0->SetHeader(Form("#font[42]{Range of fit: %.1f<p_{T}<%.1f}",min[PID_code],max[PID_code]));
	    legend0->AddEntry( f_Draw_mean[0], "#font[42]{fit (p0*Power(1.0/x, p1)+p2*x+p3)}","l");

	    pad[0]->cd();
	    gr_Mean[0]->Draw("APL");
	   	f_Draw_mean[0]->Draw("same");
	    legend0->Draw("same");


	    //h1_fit[0]->SetTitle(Form("%s Pos, %.1f<p_{T}<%.1f", particles[PID_code],ptBinRange[pt_bin], ptBinRange[pt_bin+1] ));
	    
	    TLegend *legend1 = new TLegend(0.5,0.4,0.9,0.65);
	    //legend1->SetHeader("dgdjgjdj");
	    legend1->SetHeader(Form("#font[42]{Range of fit: %.1f<p_{T}<%.1f}",min[PID_code],max[PID_code]));
	    legend1->AddEntry( f_Draw_mean[0], "#font[42]{fit (p0*Power(1.0/x, p1)+p2*x+p3)}","l");

	    pad[1]->cd();
	  	gr_Sigma[0]->Draw("APL");
		f_Draw_sigma[0]->Draw("same");
	    legend1->Draw("same");


	    
	    pad[2]->cd();
	   	gr_Mean[1]->Draw("APL");
	   	f_Draw_mean[1]->Draw("same");
	    //h1_fit[0]->SetTitle(Form("%s Pos, %.1f<p_{T}<%.1f", particles[PID_code],ptBinRange[pt_bin], ptBinRange[pt_bin+1] ));
	    
	    pad[3]->cd();
	  	gr_Sigma[1]->Draw("APL");
	  	f_Draw_sigma[1]->Draw("same");
	    
	    canvas->SaveAs(Form("./pict_proton/ParFitNSIgma_%s_Run18.png",particles[PID_code]));
	}



}

void FirstFitTree2x2DGaus(TFile *inFileHisto, TFile *inFileFitFun, TFile *outFile, Int_t pt){
	
	TH2D *histoPos = (TH2D*)inFileHisto -> Get(Form("h2_m2VsnSigma_piKp_0_pt%i", pt));
	TH2D *histoNeg = (TH2D*)inFileHisto -> Get(Form("h2_m2VsnSigma_piKp_1_pt%i", pt));
	
	//histoPos->RebinX(2);
	//histoNeg->RebinX(2);

	FitTree2x2DGaus(inFileFitFun,outFile,histoPos,pt,0);
	FitTree2x2DGaus(inFileFitFun,outFile,histoNeg,pt,1);
}

void FitTree2x2DGaus(TFile *inFile, TFile *outFile, TH2D *histo, Int_t pt, Int_t charge){

	gStyle->SetOptFit(kTRUE);

	if(pt>2){
		histo->Rebin2D(2);
	}

	Double_t factor=1.8;

	Double_t pt_point = (m[pt]+m[pt+1])/2;
	Double_t parThree2DGaus[15] = {0.};
	Double_t parThree2x2DGaus[30] = {0.};

	Double_t limitPar_meanM2[3]={0.02,0.04,0.1};
	Double_t limitPar_meanNSigma[3]={0.05*factor,0.2*factor,0.3*factor};

	TF2 *Three2DGaus = new TF2(Form("Three2DGaus%ich%i",pt,charge),"xygaus(0)+xygaus(5)+xygaus(10)", -3.*factor,6.*factor,-0.3,0.95);
	TF2 *Three2x2DGaus = new TF2(Form("Three2x2DGaus%ich%i",pt,charge),"xygaus(0)+xygaus(5)+xygaus(10)+xygaus(15)+xygaus(20)+xygaus(25)", -3*factor,6*factor,-0.4,1.4);

	inFile->cd();

	for(Int_t par=0; par<3; par++){

		Three2DGaus->SetParameter(0+5*par,100);
		Three2DGaus->SetParameter(1+5*par, readParFromFun(inFile, Form("f1_meanNSigma_%s_ch%i",particles[par],charge),pt_point));
		Three2DGaus->SetParameter(2+5*par, readParFromFun(inFile, Form("f1_sigmaNSigma_%s_ch%i",particles[par],charge),pt_point));
		Three2DGaus->SetParameter(3+5*par, readParFromFun(inFile, Form("f1_meanM2_%s_ch%i",particles[par],charge), pt_point) );
		Three2DGaus->SetParameter(4+5*par, readParFromFun(inFile, Form("f1_sigmaM2_%s_ch%i",particles[par],charge), pt_point) );

		Three2DGaus->SetParLimits(1+5*par, readParFromFun(inFile, Form("f1_meanNSigma_%s_ch%i",particles[par],charge),pt_point)-limitPar_meanNSigma[par],
										   readParFromFun(inFile, Form("f1_meanNSigma_%s_ch%i",particles[par],charge),pt_point)+limitPar_meanNSigma[par] );
		Three2DGaus->SetParLimits(2+5*par, 0.4*readParFromFun(inFile, Form("f1_sigmaNSigma_%s_ch%i",particles[par],charge),pt_point) ,
										   1.6*readParFromFun(inFile, Form("f1_sigmaNSigma_%s_ch%i",particles[par],charge),pt_point));
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

void FirstFitTree2x2DGausNewXY(TFile *inFileHisto, TFile *inFileFitFun, TFile *outFile, Int_t pt){
	
	TH2D *histoPos = (TH2D*)inFileHisto -> Get(Form("h2_m2VsnSigma_piKp_0_pt%i_new", pt));
	TH2D *histoNeg = (TH2D*)inFileHisto -> Get(Form("h2_m2VsnSigma_piKp_1_pt%i_new", pt));
	
	//histoPos->RebinX(2);
	//histoNeg->RebinX(2);

	FitTree2x2DGausNewXY(inFileFitFun,outFile,histoPos,pt,0);
	FitTree2x2DGausNewXY(inFileFitFun,outFile,histoNeg,pt,1);
}

void FitTree2x2DGausNewXY(TFile *inFile, TFile *outFile, TH2D *histo, Int_t pt, Int_t charge){

	gStyle->SetOptFit(kTRUE);

	Double_t pt_point = (m[pt]+m[pt+1])/2;
	Double_t parThree2DGaus[15] = {0.};
	Double_t parThree2x2DGaus[30] = {0.};

	Double_t limitPar_meanNewX[3]={0.1,0.1,0.2};
	Double_t limitPar_meanNewY[3]={0.1,0.1,0.2};

	TH1D *histoProjX = (TH1D*)histo->ProjectionX(Form("px_ptBin%i_charge%i",pt,charge));

	TF1 *tf1_Gaus;
	TF1 *One1DGaus = new TF1(Form("One1DGaus%ich%i",pt,charge),"gaus", -1.0,0.1);
	TF1 *Three1DGaus = new TF1(Form("Three1DGaus%ich%i",pt,charge),"gaus(0)+gaus(3)+gaus(6)", -0.4,1.5);
	TF2 *Three2DGaus = new TF2(Form("Three2DGaus%ich%iNewXY",pt,charge),"xygaus(0)+xygaus(5)+xygaus(10)", -0.3,1.2,-1.2,0.5);
	TF2 *Three2x2DGaus = new TF2(Form("Three2x2DGaus%ich%iNewXY",pt,charge),"xygaus(0)+xygaus(5)+xygaus(10)+xygaus(15)+xygaus(20)+xygaus(25)", -0.5,1.5,-1.3,0.6);

	Double_t meanNewX1d[3]={0.0,0.23,0.87};
	Double_t meanNewX1min[3]={-0.2,0.18,0.82};
	Double_t meanNewX1smax[3]={0.15,0.30,1.1};

	for(Int_t par=0; par<3; par++){
		tf1_Gaus = fitFun1DGaus(histoProjX , 10., meanNewX1d[par], 0.01, meanNewX1min[par], meanNewX1smax[par], 0);
		Three1DGaus->SetParameter(0+3*par,tf1_Gaus->GetParameter(0));
		Three1DGaus->SetParameter(1+3*par,tf1_Gaus->GetParameter(1));
		Three1DGaus->SetParameter(2+3*par,tf1_Gaus->GetParameter(2));
		delete tf1_Gaus;
	}
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

	TCanvas *can = new TCanvas(Form("canvas%i%i",pt,charge),"plot",1024,1024);
	histo->Draw("colz");
	Three2x2DGaus->Draw("same");
	//can->SaveAs(Form("./FitThree2x2DGausNewXY_pt%i_ch%i_Run18.png",pt,charge));

	outFile->cd();

	histoProjX->Write();
	histoProjY->Write();
	histo->Write();
	Three2DGaus->Write();
	Three2x2DGaus->Write();


}

void WriteGausFit(TFile *inFile, TFile *outFile, Int_t pt, TString mode){

	TF1 *f1_15;

	for(Int_t ch=0; ch<2; ch++){

		inFile->cd();
		std::cout<<Form("Three2DGaus%ich%i%s",pt,ch,mode.Data())<<"\n";
		f1_15 = (TF1*)inFile -> Get(Form("Three2DGaus%ich%i%s",pt,ch,mode.Data()));
		outFile->cd();
		f1_15->Write();
	}

	for(Int_t ch=0; ch<2; ch++){

		inFile->cd();
		f1_15 = (TF1*)inFile -> Get(Form("Three2x2DGaus%ich%i%s",pt,ch,mode.Data()));
		outFile->cd();
		f1_15->Write();
	}

}


void GetParFit(TFile *outFile){
	
	std::vector<Double_t> par;
	
	TF2 *fun[point][2];
	TFile *file[point];
	TCanvas *can[point];
	
	std::vector<Double_t> ptBinRange;
	for(Int_t i=0; i<point; i++){
		ptBinRange.push_back(m[i]);
	}
	
	gStyle->SetOptFit(kTRUE);

	for(Int_t i=0; i<point; i++){
		can[i] = new TCanvas(Form("canvas%i",i),"plot",1024,1024);
	}

	for(Int_t i=0; i<ptBinRange.size()-1;i++){
		file[i] = new TFile(Form("./fitGaus/27GeV_new_ptBin_%i.root",i), "READ");
		fun[i][0] = (TF2*)file[i]->Get(Form("Three2x2DGaus%ich0",i));
		fun[i][1] = (TF2*)file[i]->Get(Form("Three2x2DGaus%ich1",i));
	}

	Double_t err[point]={0.0};
	Double_t pt[point]={0.};
	Double_t meanSigma[2][3][2][point];
	Double_t widthSigma[2][3][2][point];
	Double_t meanM2[2][3][2][point];
	Double_t widthM2[2][3][2][point];
	Int_t corr=0;

	
	for(Int_t ch=0; ch<2; ch++){
		for(Int_t par=0; par<3; par++){
			for(Int_t p=0; p<ptBinRange.size()-1; p++){
				
				meanSigma[ch][par][0][p]=fun[p][ch]->GetParameter(1+par*5);
				widthSigma[ch][par][0][p]=fun[p][ch]->GetParameter(2+par*5);
				
				meanM2[ch][par][0][p]=fun[p][ch]->GetParameter(3+par*5);
				widthM2[ch][par][0][p]=abs(fun[p][ch]->GetParameter(4+par*5));

				pt[p]= (ptBinRange[p]+ptBinRange[p+1])/2.0;
				
			}
		}
	}
	

	TGraphErrors *graph[3][4][2];
	TString txt_charge[] = {"Pos","Neg"};
	TString txt_par[] = {"Pion","Kaon","Proton"};
	TString txt_volum[] = {"MeanSigma","WidthSigma","MeanM2","WidthM2"};
	TString txt_volum2[] = {"#mu(n#sigma_{#pi})","#sigma(n#sigma_{#pi})","#mu(m^2)","#sigma(m^2)"};
	TString txt_title1[] = {"Mean","Width"};
	TString txt_title2[] = {"n#sigma","m^{2}"};
	TString txt_bg[] = {"Peak","BG"};
	Int_t number_can=0;

	TF1 *f_meanNSigma[3];
	TF1 *f_meanM2[3];
	TF1 *f_wigthM2[3];
	TF1 *f_wigthNSigma[3];

	f_meanM2[0] = new TF1("f_meanM2[0]","[0] + [1]*x + [2]*x*x", 0.2 , 1.8);
	f_meanM2[1] = new TF1("f_meanM2[1]","[0] + [1]*x + [2]*x*x", 0.2 , 1.8);
	f_meanM2[2] = new TF1("f_meanM2[2]","[0] + [1]*x + [2]*x*x", 0.2 , 2.5);
	
	f_wigthM2[0] = new TF1("f_wigthM2[0]","[0] + [1]*x + [2]*x*x", 0.2 , 1.8);
	f_wigthM2[1] = new TF1("f_wigthM2[1]","[0] + [1]*x + [2]*x*x", 0.2 , 1.8);
	f_wigthM2[2] = new TF1("f_wigthM2[2]","[0] + [1]*x + [2]*x*x", 0.5 , 2.5);
	f_wigthNSigma[0] = new TF1("f_wigthNSigma[0]","[0] + [1]*x + [2]*x*x", 0.2 , 1.8);
	f_wigthNSigma[1] = new TF1("f_wigthNSigma[1]","[0] + [1]*x + [2]*x*x", 0.2 , 1.8);
	f_wigthNSigma[2] = new TF1("f_wigthNSigma[2]","[0] + [1]*x + [2]*x*x + [3]*x*x*x", 0.5 , 2.5);
	
	//f_meanNSigma[0]= new TF1(Form("f_%s",txt_par[0].Data()),"[0]+[1]*x+[2]*x*x + [3]*x*x*x",0.2 , 1.6);
	f_meanNSigma[0]= new TF1(Form("f_%s",txt_par[0].Data()),"[0]*TMath::Power(1.0/x,[1]) + [2]*x + [3]",0.2 , 1.6);
	f_meanNSigma[1]= new TF1(Form("f_%s",txt_par[1].Data()),"[0]*TMath::Power(1.0/x,[1]) + [2]*x + [3]",0.2 , 1.6);
	f_meanNSigma[2]= new TF1(Form("f_%s",txt_par[2].Data()),"[0]*TMath::Power(1.0/x,[1]) + [2]*x + [3]",0.5 , 2.5);


	for(Int_t i=0; i<3; i++){
		f_meanNSigma[i]->SetLineColor(kBlack);
	}

	for(Int_t ch=0; ch<1; ch++){
		for(Int_t par=0; par<3; par++){
			for(Int_t vol=0; vol<4; vol++){
				for(Int_t cor=0; cor<2; cor++){
					
					if(vol==0)graph[par][vol][cor] = new TGraphErrors(point-2, pt, meanSigma[ch][par][cor], err, err);
					if(vol==1)graph[par][vol][cor] = new TGraphErrors(point-2, pt, widthSigma[ch][par][cor], err, err);
					if(vol==2)graph[par][vol][cor] = new TGraphErrors(point-2, pt, meanM2[ch][par][cor], err, err);
					if(vol==3)graph[par][vol][cor] = new TGraphErrors(point-2, pt, widthM2[ch][par][cor], err, err);
	          		
	          		graph[par][vol][cor]->SetName(Form("gr_%s%s_%s_%s",txt_par[par].Data(),txt_charge[ch].Data(),txt_volum[vol].Data(),txt_bg[cor].Data()));
	          		graph[par][vol][cor]->SetTitle(Form("%s %s %s",txt_title1[(int)(vol%2)].Data(),txt_title2[(int)(vol/2)].Data(),txt_par[par].Data()));
	          		if(cor==0){
	          			graph[par][vol][cor]->SetMarkerColor(kRed);
						graph[par][vol][cor]->SetLineColor(kRed);
					}else{
						graph[par][vol][cor]->SetMarkerColor(kBlue);
						graph[par][vol][cor]->SetLineColor(kBlue);
					}
					graph[par][vol][cor]->SetMarkerStyle(24);
					graph[par][vol][cor]->SetMarkerSize(2);
					graph[par][vol][cor]->SetLineWidth(2);
					outFile->cd();
					graph[par][vol][cor]->Write();
				}
				can[number_can]->cd();
				graph[par][vol][0]->Draw("ALP");
				graph[par][vol][0]->GetXaxis()->SetTitle("p_{t}[GeV/c]");
    			graph[par][vol][0]->GetYaxis()->SetTitle(Form("%s %s",txt_volum2[vol].Data(), txt_par[par].Data()));
    			graph[par][vol][0]->GetYaxis()->SetTitleOffset(1.0);

				if(vol==0){
					//graph_star_meanSigma[par]->Draw("L*same");
					graph[par][vol][0]->Fit(f_meanNSigma[par],"RM");
				}
				if(vol==1){
					f_wigthNSigma[par]->SetLineColor(kBlack);
					graph[par][vol][0]->Fit(f_wigthNSigma[par],"RM");
					//graph[par][vol][0]->Fit(f_wigthNSigma[par],"RM","",0.4,2.6);
				}
				if(vol==2){
					f_meanM2[par]->SetLineColor(kBlack);
					graph[par][vol][0]->Fit(f_meanM2[par],"RM");
				}
				if(vol==3){
					f_wigthM2[par]->SetLineColor(kBlack);
					graph[par][vol][0]->Fit(f_wigthM2[par],"RM");
				}
				//graph[par][vol][1]->Draw("L*same");
				can[number_can]->SaveAs(Form("./new_code_pict/firstFitPar/%s%s_%s_Run18.png",txt_volum[vol].Data(),txt_par[par].Data(),txt_charge[ch].Data()));
				outFile->cd();
				can[number_can]->Write();
				number_can++;
			}
		}
	}

}
