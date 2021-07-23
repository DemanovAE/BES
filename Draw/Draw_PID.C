#include <string>
#include <iostream>
#include <vector>
#include <TGraph.h>
#include <Math/SpecFuncMathMore.h>
#include <TMath.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include "TFile.h"
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>
typedef vector<Double_t> VEC;
//#include "./Constants.h"
#include "./Analysis_Function.h"

const Char_t *outPath = "/home/demanov/New_work/MyAnalysis";
const std::map<Int_t, const Char_t *> Energy_scan_text={{7, "7.7"},{11, "11.5"},{14, "14.5"},{19, "19.6"},{27, "27"},{39, "39"},{62, "62.4"}};
TString particle[]={"#pi^{+}","#pi^{-}","K^{+}","K^{-}","p","#bar{p}","#pi^{+}#pi^{-}","K^{+}K^{-}","p#bar{p}","Hadrons","","PID"};
const Int_t cent[10]={80,70,60,50,40,30,20,10,5,0};
TString bukva[] = {"(a)","(b)","(c)","(d)","(i)", "(f)", "(g)","(h)", "(i)", "(j)","(k)","(l)","(m)","(n)","(o)","(p)", "(q)"};

void makeplotstyle();
int GetNumberParticle(TString PID, TString charge);
VEC Rebin2(Int_t harmonic, Int_t Energy, TString PID, TString charge, Int_t CentBinMin, Int_t CentBinMax);
VEC Rebin2K(Int_t harmonic, Int_t Energy, TString PID, TString charge, Int_t CentBinMin, Int_t CentBinMax);
VEC Rebin2H(Int_t harmonic, Int_t Energy, TString PID, TString charge, Int_t CentBinMin, Int_t CentBinMax);
void DrawEtaPID(TFile *file, Int_t centMin, Int_t centMax, Double_t ptMin, Double_t ptMax, Int_t Energy, TString PID, TString charge, TString EtaGap);
//потоки как функция энергии
void FlowVsEnergyBESupdate(Double_t ptMin, Double_t ptMax, TString EtaGap, TString PID, TString charge, Double_t min_y_axis,  Double_t max_y_axis,Double_t min_x_axis, Double_t max_x_axis,Double_t MinYFlow, Double_t MaxYratio);
void DrawFlowVsRapidity(Int_t Energy);
void FlowVsRapidity(Int_t Energy, Int_t harmonic, Int_t CentBinMin, Int_t CentBinMax, Double_t ptMin, Double_t ptMax, TString EtaGap, TString PID, TString charge, Double_t min_y_axis,Double_t max_y_axis);
void FlowVsRapiditySymmetry(Int_t Energy, Int_t harmonic, Double_t ptMin, Double_t ptMax, TString EtaGap, TString PID, TString charge,  Double_t min_y_axis, Double_t max_y_axis);
void FlowVsRapidityForBES(Int_t harmonic, Double_t ptMin, Double_t ptMax, TString EtaGap, TString PID, TString charge,  Double_t min_y_axis, Double_t max_y_axis);
//void FlowVsRapidityEnergy(Int_t CentBinMin, Int_t CentBinMax, Double_t ptMin, Double_t ptMax, TString EtaGap, TString PID, TString charge, Double_t max_y_axis, Double_t max_x_axis,Double_t MinYratio, Double_t MaxYratio);
//потоки как функция pt в зависимости от энергий 
void FlowVsPtForBES(Int_t harmonic, Int_t CentBinMin, Int_t CentBinMax, Double_t ptMin, Double_t ptMax, TString EtaGap, TString PID, TString charge, Double_t max_y_axis, Double_t max_y_axis2, Double_t max_x_axis, Double_t MinYratio, Double_t MaxYratio,Double_t MinYratio2, Double_t MaxYratio2);
void DrawFlowVsPtForBES(TString mod);
//Разные eta-gap с ратио
void FlowDifferentEtaGap(Int_t Energy, Int_t harmonic, Int_t CentBinMin, Int_t CentBinMax, TString PID, TString charge, Double_t max_y_axis, Double_t max_x_axis, Double_t MinYratio, Double_t MaxYratio);
void DrawFlowDifferentEtaGap(TString mod, Int_t Energy);
//разные методы PID
void FlowDifferentMethodPID(Int_t Energy, Int_t harmonic, Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, TString PID, TString charge, Double_t max_y_axis, Double_t max_x_axis, Double_t MinYratio, Double_t MaxYratio);
void DrawFlowDifferentMethodPID(TString mod, Int_t Energy);
//разница частиц и античастиц в зависимости от энергии
void ParAntVsEnergyMultyPad( Int_t CentBinMin, Int_t CentBinMax, Double_t ptMin, Double_t ptMax, TString EtaGap,  Double_t min_x_axis, Double_t max_x_axis, Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio, Double_t min_y_axis_dif, Double_t max_y_axis_dif);
void FlowVsPtForBESmultyPad( Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, TString charge, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio, Double_t min_y_axis_dif, Double_t max_y_axis_dif);
void FlowVsPtForBESmultyPad_Int( Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, TString charge, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio, Double_t min_y_axis_dif, Double_t max_y_axis_dif);
void ParAntFlowVsPtForBES( Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, TString PID, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio, Double_t min_y_axis_dif, Double_t max_y_axis_dif);
void ParAntFlowVsCentForBES( Int_t ptMin, Int_t ptMax, TString EtaGap, TString PID, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio, Double_t min_y_axis_dif, Double_t max_y_axis_dif);
void ParAntFlowVsPtForBESupdate( Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, Double_t min_x_axis, Double_t max_x_axis);
void DrawFLowVsEnergyBES();

void PIDFlowVsPtVsCentForBES_Int( TString EtaGap, TString PID, TString charge, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis);
void PIDFlowVsPtVsCentForBES_noInt( TString EtaGap, TString PID, TString charge, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis);
void PIDFlowVsPtVsCentForBES_Int_cent060( TString EtaGap, TString PID, TString charge, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis);

void PrintPoint(Int_t Energy, Int_t harmonic, Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, TString PID, TString charge);
void funPrintPoint(Int_t Energy);

void DrawSravOldNew27(Int_t Energy, Int_t harmonic, TString EtaGap, TString PID, Double_t min_x_axis, Double_t max_x_axis, Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_dif, Double_t max_y_axis_dif, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio);

void Resolution();


int Analysis_BES(Int_t Energy){
	

	makeplotstyle();
	TFile *fstyle = new TFile("style.root");
	TStyle *tsty = (TStyle *)fstyle->Get("PlottingInStyle");
	tsty->cd();

	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	
	/*
	PIDFlowVsPtVsCentForBES_Int("Eta01", "Pion","Pos",-0.05, 2.86, -0.3, 5.67);
	PIDFlowVsPtVsCentForBES_Int("Eta01", "Pion","Neg",-0.05, 2.86, -0.3, 5.67);
	
	PIDFlowVsPtVsCentForBES_Int("Eta01", "Kaon","Pos",-0.05, 2.86, -0.3, 5.67);
	PIDFlowVsPtVsCentForBES_Int("Eta01", "Kaon","Neg",-0.05, 2.86, -0.3, 5.67);
	
	PIDFlowVsPtVsCentForBES_Int("Eta01", "Proton","Pos",-0.05, 2.86, -0.3, 5.67);
	PIDFlowVsPtVsCentForBES_Int("Eta01", "Proton","Neg",-0.05, 2.86, -0.3, 5.67);
	*/
	/*
	PIDFlowVsPtVsCentForBES_noInt("Eta01", "Pion","Pos",-0.05, 2.86, -0.02, 0.27);
	PIDFlowVsPtVsCentForBES_noInt("Eta01", "Pion","Neg",-0.05, 2.86, -0.02, 0.27);
	
	PIDFlowVsPtVsCentForBES_noInt("Eta01", "Kaon","Pos",-0.05, 2.86, -0.02, 0.27);
	PIDFlowVsPtVsCentForBES_noInt("Eta01", "Kaon","Neg",-0.05, 2.86, -0.02, 0.27);
	
	PIDFlowVsPtVsCentForBES_noInt("Eta01", "Proton","Pos",-0.05, 2.86, -0.02, 0.27);
	PIDFlowVsPtVsCentForBES_noInt("Eta01", "Proton","Neg",-0.05, 2.86, -0.02, 0.27);
	*/
	//funPrintPoint(Energy);

	//FlowVsEnergyBESupdate( 0.2, 2.0, "Eta15","Hadrons", "", -0.005, 0.115, 4, 75, 0.5, 1.15);
	//FlowVsEnergyBESupdate( 0.2, 3.2, "Eta15","Hadrons", "", -0.005, 0.115, 4, 75, 0.5, 1.15);
	/*
	DrawSravOldNew27(27, 2, "Eta01", "Pion", -0.06, 2.76, -0.02, 0.24, -0.03, 0.034, 0.82, 1.18);
	DrawSravOldNew27(27, 2, "Eta01", "Kaon", -0.06, 2.76, -0.02, 0.24, -0.03, 0.034, 0.82, 1.18);
	DrawSravOldNew27(27, 2, "Eta01", "Proton", -0.06, 2.76, -0.02, 0.32, -0.03, 0.034, 0.97, 2.18);

	DrawSravOldNew27(27, 3, "Eta01", "Pion", -0.06, 2.76, -0.002, 0.095, -0.03, 0.034, 0.82, 1.18);
	DrawSravOldNew27(27, 3, "Eta01", "Kaon", -0.06, 2.76, -0.002, 0.134, -0.03, 0.034, 0.82, 1.18);
	DrawSravOldNew27(27, 3, "Eta01", "Proton", -0.06, 2.76, -0.002, 0.115, -0.03, 0.034, 0.97, 2.18);
	*/
    //ParAntVsEnergyMultyPad( Int_t CentBinMin, Int_t CentBinMax, Double_t ptMin, Double_t ptMax, TString EtaGap,  Double_t min_x_axis, Double_t max_x_axis, Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio, Double_t min_y_axis_dif, Double_t max_y_axis_dif);
	//ParAntVsEnergyMultyPad(2,8,0.2,2.0,"Eta01", 6.8, 68, -0.005, 0.035, 0.82, 1.18, -0.03, 0.034);good
	//ParAntVsEnergyMultyPad(2,8,0.2,2.0,"Eta01", 6.8, 68, -0.006, 0.014, 0.82, 1.18, -0.03, 0.034);
	
	
 	//FlowVsPtForBESmultyPad( Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, TString charge, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio, Double_t min_y_axis_dif, Double_t max_y_axis_dif);
	//FlowVsPtForBESmultyPad(2, 8, "Eta01", "Pos",-0.04, 2.96, -0.02, 0.22, 0.82, 1.18, -0.03, 0.034);
	//FlowVsPtForBESmultyPad(2, 8, "Eta01", "Neg",-0.04, 2.96, -0.02, 0.22, 0.82, 1.18, -0.03, 0.034);

   	//FlowVsPtForBESmultyPad_Int(2, 8, "Eta01", "Pos",-0.04, 2.96, -0.3, 4.2, 0.82, 1.18, -0.03, 0.034);
	//FlowVsPtForBESmultyPad_Int(2, 8, "Eta01", "Neg",-0.04, 2.96, -0.3, 4.2, 0.82, 1.18, -0.03, 0.034);

	///////////////////////ParAntFlowVsPtForBESupdate(2, 8, "Eta01",-0.04, 2.96);
	ParAntFlowVsPtForBESupdate(2, 8, "Eta01", 0.04, 2.96);
	//ParAntFlowVsPtForBESupdate(2, 8, "Eta01", 0.04, 2.96);
  // ParAntFlowVsPtForBESupdate(7, 8, "Eta01",-0.04, 2.96);
  	//ParAntFlowVsPtForBESupdate(2, 6, "Eta01",-0.04, 2.96);
    //ParAntFlowVsPtForBES(2, 8, "Eta01", "Pion",-0.04, 2.96, -0.02, 0.22, 0.82, 1.18, -0.03, 0.034);

	/*
    ParAntFlowVsPtForBES(2, 8, "Eta01", "Kaon",-0.05, 2.26, -0.02, 0.22, 0.81, 1.19, -0.013, 0.013);
    ParAntFlowVsPtForBES(0, 8, "Eta01", "Kaon",-0.05, 2.26, -0.02, 0.22, 0.81, 1.19, -0.013, 0.013);
    ParAntFlowVsPtForBES(0, 3, "Eta01", "Kaon",-0.05, 2.26, -0.02, 0.22, 0.81, 1.19, -0.024, 0.024);
    ParAntFlowVsPtForBES(4, 6, "Eta01", "Kaon",-0.05, 2.26, -0.02, 0.22, 0.81, 1.19, -0.013, 0.013);
    ParAntFlowVsPtForBES(7, 8, "Eta01", "Kaon",-0.05, 2.26, -0.02, 0.22, 0.81, 1.19, -0.013, 0.013);
	*/

    //ParAntFlowVsPtForBES(2, 6, "Eta01", "Kaon",-0.05, 2.96, -0.02, 0.22, 0.72, 1.48, -0.03, 0.034);

    /*
    ParAntFlowVsPtForBES(2, 8, "Eta01", "Pion",-0.05, 2.96, -0.02, 0.22, 0.82, 1.18, -0.03, 0.034);
    ParAntFlowVsPtForBES(2, 8, "Eta01", "Kaon",-0.05, 2.96, -0.02, 0.22, 0.72, 1.48, -0.03, 0.034);
    ParAntFlowVsPtForBES(2, 8, "Eta01", "Proton",-0.05, 2.96, -0.02, 0.22, 0.89, 2.07, -0.01, 0.084);

    ParAntFlowVsPtForBES(2, 6, "Eta01", "Pion",-0.05, 2.96, -0.02, 0.22, 0.82, 1.18, -0.03, 0.034);
    ParAntFlowVsPtForBES(2, 6, "Eta01", "Kaon",-0.05, 2.96, -0.02, 0.22, 0.72, 1.48, -0.03, 0.034);
    ParAntFlowVsPtForBES(2, 6, "Eta01", "Proton",-0.05, 2.96, -0.02, 0.22, 0.89, 2.07, -0.01, 0.084);
 	*/
    //FlowVsRapidityForBES(, Double_t ptMin, Double_t ptMax, TString EtaGap, TString PID, TString charge,  Double_t min_y_axis, Double_t max_y_axis){
 	//FlowVsRapidityForBES(2,0.2,2.0,"Eta15","Hadrons", "",0.015,0.065);
 	//FlowVsRapidityForBES(2,0.2,3.2,"Eta15","Hadrons", "",0.015,0.065);
 	//FlowVsRapidityForBES(3,0.2,2.0,"Eta15","Hadrons", "",0.001,0.023);
 	//FlowVsRapidityForBES(3,0.2,3.2,"Eta15","Hadrons", "",0.001,0.023);
	//DrawFlowVsRapidity(Energy);
	//DrawFlowDifferentEtaGap("Hadrons", Energy);
	//DrawFlowDifferentEtaGap("PID", Energy);
	DrawFlowDifferentMethodPID("PID", Energy);
	//DrawFlowVsPtForBES("Hadrons");

	//Resolution();

	return 0;

}


void Resolution(){

	TFile *file_1 = new TFile(Form("%s/file_pid/Flow_11GeV_PID_10binYesTofdca11_sys.root",outPath),"READ");
	TFile *file_out = new TFile("OUT.root","RECREATE");

	Draw_Object *res2;
	Draw_Object *res3;

	Analysis_Function *Analysis = new Analysis_Function;
	
	Analysis -> SetParametrs(27, file_1, file_out, 2, "", "", "Eta01", "");
	Analysis -> ResolutionEP_EtaSub();

	Analysis -> SetParametrs(27, file_1, file_out, 3, "", "", "Eta01", "");
	Analysis -> ResolutionEP_EtaSub();



}


void PIDFlowVsPtVsCentForBES_noInt( TString EtaGap, TString PID, TString charge, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis){

	VEC rebin_vec;
	VEC rebin_vec_integral={0.2,2.0,5.2};
	if(strncmp(PID,"Proton",6)==0){
		rebin_vec_integral={0.5,2.0,5.2};	
	}
	TString prefix = "TPCandTOF";
	
	Int_t Arr_energy[6]={11,14,19,27,39,62};
	const Int_t color[]={6, 1, 4, 2, 1, 2, 46};
	const Int_t style[]={8, 22, 23, 24, 25, 34, 8, 29};
	const Double_t size[]={1.5,1.5,1.5,1.5,1.5,1.5,1.5};
	
	Int_t centBinMin[6] = {7,6,5,4,2};
	Int_t centBinMax[6] = {8,6,5,4,3};

	Float_t k[2]={1.0,2.0};

	TFile *file_1[6];
	for(Int_t i=0; i<6; i++){
		file_1[i] = new TFile(Form("%s/file_pid/Flow_%iGeV_PID_10binYesTofdca11_sys.root",outPath, Arr_energy[i]),"READ");
		//file_1[i] = new TFile(Form("%s/file_hadrons/Sys/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2_sys.root",outPath, Arr_energy[i]),"READ");
	}
	TFile *file_out = new TFile("OUT.root","RECREATE");

	Draw_Picture_new *canvasFlow = new Draw_Picture_new;
	Draw_Picture_new *canvasDifference = new Draw_Picture_new;
	Draw_Picture_new *canvasRatio = new Draw_Picture_new;
	
	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<std::vector<Draw_Object *>> flow;
	std::vector<Draw_Object *> object;
	
	Int_t b=0;

	for(Int_t c=0; c<5; c++){
		
		canvasFlow->SetTextInfo(Form("#font[42]{#scale[1.0]{%i-%i%%}}",cent[centBinMax[c]+1], cent[centBinMin[c]] ), "Pad", 0.35, 0.85*max_y_axis, 0.13, 0.0);
		
		for(Int_t harm=0; harm<2; harm++){
			for(Int_t j=0; j<6; j++){

				rebin_vec = Rebin2H( harm+2, Arr_energy[j], PID, charge, centBinMin[c], centBinMax[c]);
				
				object.push_back(new Draw_Object);
				Analysis -> SetParametrs(Arr_energy[j], file_1[j], file_out, harm+2, PID, charge, EtaGap, prefix);
				Analysis -> FlowVsPt_EtaSub(object.back(), centBinMin[c], centBinMax[c], rebin_vec);
	    		//Analysis -> FlowVsPtDivideIntegralFLow_EtaSub_PID(object.back(), centBinMin[c], centBinMax[c], rebin_vec, rebin_vec_integral, 0);
	    		object.back() -> SetParametrsGraph(Form("data_v%i",harm+2),Form("#font[42]{#scale[1.3]{%s GeV}}", Energy_scan_text.at(Arr_energy[j])), color[j], style[j], size[j]);
				if(harm==1){
					for(Int_t n=0; n < object.back()->GetSizeVector(); n++){
						object.back()->СhangePointGraph(n,
														object.back()->GetPointXGraph(n), 
														object.back()->GetPointYGraph(n) * k[1],
														object.back()->GetPointXErrorGraph(n),
														object.back()->GetPointYErrorGraph(n)* k[1],
														object.back()->GetPointYSysErrorGraph(n) * k[1]);
					}
				}

			}

			flow.push_back(object);
			object.clear();

			if(c==4){
				if(harm==0) canvasFlow->SetTextInfo(Form("#font[42]{(v_{%i})}",harm+2), Form("V%i",harm+2), 0.06, 0.7*max_y_axis, 0.13, 0.0);
				if(harm==1) canvasFlow->SetTextInfo(Form("#font[42]{(%.1f #times v_{%i})}",k[1],harm+2), Form("V%i",harm+2), 0.06, 0.7*max_y_axis, 0.13, 0.0);
			}

			canvasFlow->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.85*max_y_axis, 0.13, 0.0);
			b++;
		}
	std::cout<<"good\n";

	}

	canvasFlow->SetTextInfo("#font[42]{v_{n}}","Y",0.25,0.6,0.5,0.0);
	canvasFlow->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %s}}", particle[GetNumberParticle(PID,charge)].Data()),"Info",1.3,0.85*max_y_axis,0.13,0.0);
	//canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{ STAR Preliminary}}"),"Prel",.4,0.85*max_y_axis,0.13,0.0);
	canvasFlow->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasFlow->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.13);
	canvasFlow->SetDrawObject(flow[0],"OnePad");
	canvasFlow->SetLegend(0.12,0.5,0.98,0.82,3);
	TCanvas *result = canvasFlow->CanvasNxM(1000,1280, 5, 2, 0.4, 0.575, flow,1,0,0);

	result->SaveAs(Form("%s/picture/3_PID/Sys/int/SysBES_vn_%s%s_%s_pt_energy_0_60.pdf",outPath, PID.Data(), charge.Data(), EtaGap.Data()));
	result->SaveAs(Form("%s/picture/3_PID/Sys/int/SysBES_vn_%s%s_%s_pt_energy_0_60.png",outPath, PID.Data(), charge.Data(), EtaGap.Data()));
	delete result;

	for(Int_t i=0; i<6; i++){
		file_1[i]->Close();
		delete file_1[i];
	}

}


void PIDFlowVsPtVsCentForBES_Int( TString EtaGap, TString PID, TString charge, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis){

	VEC rebin_vec;
	VEC rebin_vec_integral={0.2,2.0,5.2};
	if(strncmp(PID,"Proton",6)==0){
		rebin_vec_integral={0.5,2.0,5.2};	
	}
	TString prefix = "TPCandTOF";
	
	Int_t Arr_energy[6]={11,14,19,27,39,62};
	const Int_t color[]={6, 1, 4, 2, 1, 2, 46};
	const Int_t style[]={8, 22, 23, 24, 25, 34, 8, 29};
	const Double_t size[]={1.5,1.5,1.5,1.5,1.5,1.5,1.5};
	
	Int_t centBinMin[6] = {7,6,5,4,2};
	Int_t centBinMax[6] = {8,6,5,4,3};

	Float_t k[2]={1.0,1.0};

	TFile *file_1[6];
	for(Int_t i=0; i<6; i++){
		file_1[i] = new TFile(Form("%s/file_pid/Flow_%iGeV_PID_10binYesTofdca11_sys.root",outPath, Arr_energy[i]),"READ");
		//file_1[i] = new TFile(Form("%s/file_hadrons/Sys/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2_sys.root",outPath, Arr_energy[i]),"READ");
	}
	TFile *file_out = new TFile("OUT.root","RECREATE");

	Draw_Picture_new *canvasFlow = new Draw_Picture_new;
	Draw_Picture_new *canvasDifference = new Draw_Picture_new;
	Draw_Picture_new *canvasRatio = new Draw_Picture_new;
	
	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<std::vector<Draw_Object *>> flow;
	std::vector<Draw_Object *> object;
	
	Int_t b=0;

	for(Int_t c=0; c<5; c++){
		
		canvasFlow->SetTextInfo(Form("#font[42]{#scale[1.0]{%i-%i%%}}",cent[centBinMax[c]+1], cent[centBinMin[c]] ), "Pad", 0.35, 0.85*max_y_axis, 0.13, 0.0);
		
		for(Int_t harm=0; harm<2; harm++){
			for(Int_t j=0; j<6; j++){

				rebin_vec = Rebin2H( harm+2, Arr_energy[j], PID, charge, centBinMin[c], centBinMax[c]);
				object.push_back(new Draw_Object);
				Analysis -> SetParametrs(Arr_energy[j], file_1[j], file_out, harm+2, PID, charge, EtaGap, prefix);
				//Analysis -> FlowVsPt_EtaSub(object.back(), centBinMin[c], centBinMax[c], rebin_vec);
	    		Analysis -> FlowVsPtDivideIntegralFLow_EtaSub_PID(object.back(), centBinMin[c], centBinMax[c], rebin_vec, rebin_vec_integral, 0);
	    		object.back() -> SetParametrsGraph(Form("data_v%i_%iGeV",harm+2,Arr_energy[j]),Form("#font[42]{#scale[1.3]{%s GeV}}", Energy_scan_text.at(Arr_energy[j])), color[j], style[j], size[j]);
				if(harm==1){
					for(Int_t n=0; n < object.back()->GetSizeVector(); n++){
						object.back()->СhangePointGraph(n,
														object.back()->GetPointXGraph(n), 
														object.back()->GetPointYGraph(n) * k[1],
														object.back()->GetPointXErrorGraph(n),
														object.back()->GetPointYErrorGraph(n)* k[1],
														object.back()->GetPointYSysErrorGraph(n) * k[1]);
					}
				}

			}

			flow.push_back(object);
			object.clear();

			if(c==4){
				if(harm==0) canvasFlow->SetTextInfo(Form("#font[42]{(v_{%i})}",harm+2), Form("V%i",harm+2), 0.06, 0.7*max_y_axis, 0.13, 0.0);
				if(harm==1) canvasFlow->SetTextInfo(Form("#font[42]{(v_{%i})}",harm+2), Form("V%i",harm+2), 0.06, 0.7*max_y_axis, 0.13, 0.0);
			}

			canvasFlow->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.85*max_y_axis, 0.13, 0.0);
			b++;
		}
	std::cout<<"good\n";

	}

	canvasFlow->SetTextInfo("#font[42]{#frac{v_{n}}{v_{n}^{int}}}","Y",0.25,0.6,0.5,0.0);
	canvasFlow->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %s}}", particle[GetNumberParticle(PID,charge)].Data()),"Info",1.3,0.85*max_y_axis,0.13,0.0);
	//canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{ STAR Preliminary}}"),"Prel",.4,0.85*max_y_axis,0.13,0.0);
	canvasFlow->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasFlow->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.13);
	canvasFlow->SetDrawObject(flow[0],"OnePad");
	canvasFlow->SetLegend(0.12,0.5,0.98,0.82,3);
	TCanvas *result = canvasFlow->CanvasNxM(1000,1280, 5, 2, 0.4, 0.575, flow,1,0,2);

	result->SaveAs(Form("%s/picture/3_PID/Sys/int/intSysBES_vn_%s%s_%s_pt_energy_0_60.pdf",outPath, PID.Data(), charge.Data(), EtaGap.Data()));
	result->SaveAs(Form("%s/picture/3_PID/Sys/int/intSysBES_vn_%s%s_%s_pt_energy_0_60.png",outPath, PID.Data(), charge.Data(), EtaGap.Data()));
	result->SaveAs(Form("%s/picture/3_PID/Sys/int/intSysBES_vn_%s%s_%s_pt_energy_0_60.C",outPath, PID.Data(), charge.Data(), EtaGap.Data()));
	delete result;

	for(Int_t i=0; i<6; i++){
		file_1[i]->Close();
		delete file_1[i];
	}

}




/*  С ратио  и без
void DrawSravOldNew27(Int_t Energy, Int_t harmonic, TString EtaGap, TString PID, Double_t min_x_axis, Double_t max_x_axis, Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_dif, Double_t max_y_axis_dif, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio){

	VEC rebin_vec;
	TString prefix;
	
	if( strncmp(PID, "Hadrons",7)==0){
		prefix="";
    }else{
    	prefix="TPCandTOF";
    }
	
	Int_t centBinMin[2][4]={{7,4,0},{7,2,2}};
	Int_t centBinMax[2][4]={{8,6,3},{8,6,8}};
	TString charge_arr[]={"Pos", "Neg"};
	TString data[]={"Old 27 GeV", "New 27 GeV"};


	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<std::vector<Draw_Object *>> flow;
	std::vector<std::vector<Draw_Object *>> Ratio;
	std::vector<Draw_Object*> object;
	std::vector<Draw_Object*> object_ratio;

	const Int_t color[]={2, 4, 1};
	const Int_t style[]={33, 24, 22};
	Float_t size[3]={1.8,1.8,1.8};

	Draw_Picture_new *canvasFlow = new Draw_Picture_new;
	Draw_Picture_new *canvasRatio = new Draw_Picture_new;

	TFile *file_read[2];
	file_read[0] = new TFile(Form("%s/file_pid/Flow_Old%iGeV_PID_10binYesTofdca11_oldPID_BadRunOld.root",outPath, Energy),"READ");;
	file_read[1] = new TFile(Form("%s/file_pid/Flow_New%iGeV_PID_10binYesTofdca11_oldPID_BadRunOld.root",outPath, Energy),"READ");;
	TFile *file_out = new TFile("OUT.root","RECREATE");	
	Int_t b=0;
	Int_t b2=0;

	Float_t 	max_y_axis_2=1.23;

	for(Int_t ch=0; ch<2;ch++){
		for(Int_t c=0; c<3;c++){
			for (Int_t f = 0; f < 2; f++){
			    //Data line
				rebin_vec = Rebin2( harmonic, Energy, PID, charge_arr[ch], centBinMin[harmonic-2][c], centBinMax[harmonic-2][c]);			    
			    object.push_back(new Draw_Object);
			    Analysis -> SetParametrs(Energy, file_read[f],file_out,harmonic,PID,charge_arr[ch],EtaGap,prefix);
			    object.back() -> SetParametrsGraph("data",Form("#font[42]{#scale[1.2]{%s}}",data[f].Data()), color[f], style[f], size[f]);
			    Analysis -> FlowVsPt_EtaSub(object.back(), centBinMin[harmonic-2][c], centBinMax[harmonic-2][c], rebin_vec);
			}

			object_ratio.push_back(new Draw_Object);
			Analysis -> RatioGraphPointToPoint(object_ratio.back(),object[0],object[1],"");
			object_ratio.back() -> SetParametrsGraph("data",Form("#font[42]{#scale[1.2]{#frac{Old 27 GeV}{New 27 GeV}}}"), color[2], style[2], size[2]);
			Ratio.push_back(object_ratio);
			object_ratio.clear();

			flow.push_back(object);
			object.clear();
			
			canvasRatio->SetNumberPad(Form("#font[42]{%s}", bukva[b2].Data()), "Pad", 0.04, 0.95*max_y_axis_2, 0.09, 0.0);
			b2++;

			if(ch==0){
				canvasFlow->SetTextInfo(Form("#font[42]{%i - %i %%}", cent[centBinMax[harmonic-2][c]+1], cent[centBinMin[harmonic-2][c]]), Form("C%i",c), 0.04, 0.73*max_y_axis, 0.09, 0.0);
			}

			canvasFlow->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.85*max_y_axis, 0.09, 0.0);
			b++;

		}

	}
	object.clear();
	for(Int_t c=0; c<3;c++){
		for (Int_t f = 0; f < 2; f++){
		    //Data line
			rebin_vec = Rebin2( harmonic, Energy, PID, "", centBinMin[harmonic-2][c], centBinMax[harmonic-2][c]);			    
		    object.push_back(new Draw_Object);
		    Analysis -> SetParametrs(Energy, file_read[f],file_out,harmonic,PID,"",EtaGap,prefix);
		    object.back() -> SetParametrsGraph("data",Form("#font[42]{%s}",data[f].Data()), color[f], style[f], size[f]);
		    Analysis -> DifferentParAnt(object.back(), centBinMin[harmonic-2][c], centBinMax[harmonic-2][c], rebin_vec);
		}


		object_ratio.push_back(new Draw_Object);
		Analysis -> RatioGraphPointToPoint(object_ratio.back(),object[0],object[1],"");
		object_ratio.back() -> SetParametrsGraph("data",Form("#font[42]{#scale[1.2]{#frac{Old 27 GeV}{New 27 GeV}}}"), color[2], style[2], size[2]);
		Ratio.push_back(object_ratio);
		object_ratio.clear();


		canvasRatio->SetNumberPad(Form("#font[42]{%s}", bukva[b2].Data()), "Pad", 0.04, 0.95*max_y_axis_2, 0.09, 0.0);
		b2++;

		flow.push_back(object);
		object.clear();

		canvasFlow->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.75*max_y_axis_dif, 0.09, 0.0);
		b++;

	}

	object.clear();
	for(Int_t c=0; c<3;c++){
		for (Int_t f = 0; f < 2; f++){
		    //Data line
			rebin_vec = Rebin2( harmonic, Energy, PID, "", centBinMin[harmonic-2][c], centBinMax[harmonic-2][c]);			    
		    object.push_back(new Draw_Object);
		    Analysis -> SetParametrs(Energy, file_read[f],file_out,harmonic,PID,"",EtaGap,prefix);
		    object.back() -> SetParametrsGraph("data",Form("#font[42]{%s}",data[f].Data()), color[f], style[f], size[f]);
		    Analysis -> RatioParAnt(object.back(), centBinMin[harmonic-2][c], centBinMax[harmonic-2][c], rebin_vec);
		}

		object_ratio.push_back(new Draw_Object);
		Analysis -> RatioGraphPointToPoint(object_ratio.back(),object[0],object[1],"");
		object_ratio.back() -> SetParametrsGraph("data",Form("#font[42]{#scale[1.2]{#frac{Old 27 GeV}{New 27 GeV}}}"), color[2], style[2], size[2]);
		Ratio.push_back(object_ratio);
		object_ratio.clear();

		canvasRatio->SetNumberPad(Form("#font[42]{%s}", bukva[b2].Data()), "Pad", 0.04, 0.95*max_y_axis_2, 0.09, 0.0);
		b2++;

		flow.push_back(object);
		object.clear();

		canvasFlow->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.95*max_y_axis_ratio, 0.09, 0.0);
		b++;

	}


    
    canvasFlow->SetTextInfo(Form("#font[42]{ v_{%i}(%s)}", harmonic, particle[GetNumberParticle(PID,"Pos")].Data() ),"Y",0.5, 0.8,0.3,90);
	canvasFlow->SetTextInfo(Form("#font[42]{ v_{%i}(%s)}", harmonic, particle[GetNumberParticle(PID,"Neg")].Data() ),"Y",0.5, 0.6,0.3,90);
	canvasFlow->SetTextInfo(Form("#font[42]{ v_{%i}(%s) - v_{%i}(%s) }", harmonic, particle[GetNumberParticle(PID,"Pos")].Data(), harmonic, particle[GetNumberParticle(PID,"Neg")].Data() ),"Y",0.5, 0.27,0.3,90);
	canvasFlow->SetTextInfo(Form("#font[42]{ #frac{v_{%i}(%s)}{v_{%i}(%s)} }", harmonic, particle[GetNumberParticle(PID,"Pos")].Data(), harmonic, particle[GetNumberParticle(PID,"Neg")].Data() ),"Y",0.5, 0.07,0.3,90);
	
	canvasFlow->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, #sqrt{s_{NN}}=27 GeV}}"),"Info",0.3,0.85*max_y_axis,0.09,0.0);
	canvasFlow->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasFlow->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasFlow->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.07);
	canvasFlow->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.07);
	canvasFlow->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_dif, max_y_axis_dif, 0.07);
	canvasFlow->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_ratio, max_y_axis_ratio, 0.07);
	canvasFlow->SetDrawObject(flow[0],"OnePad");
	canvasFlow->SetLegend(0.3,0.7,0.98,0.95,1);
	TCanvas *result = canvasFlow->CanvasNxM(1000,1280, 4, 3, 0.6, 0.46, flow, 4,0,3);

	result->SaveAs(Form("%s/picture/3_PID/New27Gev/v%i_%s_%s_cent.pdf",outPath, harmonic, PID.Data(), EtaGap.Data() ));
	delete result;


	 canvasRatio->SetTextInfo(Form("#font[42]{v_{%i}(%s)}", harmonic, particle[GetNumberParticle(PID,"Pos")].Data() ),"Y",0.5, 0.8,0.3,90);
	canvasRatio->SetTextInfo(Form("#font[42]{ v_{%i}(%s)}", harmonic, particle[GetNumberParticle(PID,"Neg")].Data() ),"Y",0.5, 0.6,0.3,90);
	canvasRatio->SetTextInfo(Form("#font[42]{ v_{%i}(%s) - v_{%i}(%s) }", harmonic, particle[GetNumberParticle(PID,"Pos")].Data(), harmonic, particle[GetNumberParticle(PID,"Neg")].Data() ),"Y",0.5, 0.27,0.3,90);
	canvasRatio->SetTextInfo(Form("#font[42]{ #frac{v_{%i}(%s)}{v_{%i}(%s)} }", harmonic, particle[GetNumberParticle(PID,"Pos")].Data(), harmonic, particle[GetNumberParticle(PID,"Neg")].Data() ),"Y",0.5, 0.07,0.3,90);
	
	canvasRatio->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasRatio->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, #sqrt{s_{NN}}=27 GeV}}"),"Info",0.3,0.85*max_y_axis,0.09,0.0);
	canvasRatio->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasRatio->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	if(harmonic==2){
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  0.77, 1.23, 0.07);
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  0.77, 1.23, 0.07);
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  0.13, 1.87, 0.07);
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  0.77, 1.23, 0.07);
	}else{
		canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  0.57, 1.43, 0.07);
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  0.57, 1.43, 0.07);
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  0.07, 2.97, 0.07);
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  0.57, 1.43, 0.07);
	}

	canvasRatio->SetDrawObject(Ratio[0],"OnePad");
	canvasRatio->SetLegend(0.2,0.7,0.8,0.98,1);
	TCanvas *result2 = canvasRatio->CanvasNxM(1000,1280, 4, 3, 0.6, 0.46, Ratio, 4,0,0);

	result2->SaveAs(Form("%s/picture/3_PID/New27Gev/Ratio_v%i_%s_%s_cent.pdf",outPath, harmonic, PID.Data(), EtaGap.Data() ));
	delete result2;


}
*/


void DrawSravOldNew27(Int_t Energy, Int_t harmonic, TString EtaGap, TString PID, Double_t min_x_axis, Double_t max_x_axis, Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_dif, Double_t max_y_axis_dif, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio){

	VEC rebin_vec;
	TString prefix;
	
	if( strncmp(PID, "Hadrons",7)==0){
		prefix="";
    }else{
    	prefix="TPCandTOF";
    }
	
	Int_t centBinMin[2][4]={{7,4,0},{7,2,2}};
	Int_t centBinMax[2][4]={{8,6,3},{8,6,8}};
	TString charge_arr[]={"Pos", "Neg"};
	TString data[]={"Old 27 GeV", "New 27 GeV"};


	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<std::vector<Draw_Object *>> flow;
	std::vector<std::vector<Draw_Object *>> Ratio;
	std::vector<Draw_Object*> object;
	std::vector<Draw_Object*> object_ratio;

	const Int_t color[]={2, 4, 1};
	const Int_t style[]={33, 24, 22};
	Float_t size[3]={1.8,1.8,1.8};

	Draw_Picture_new *canvasFlow = new Draw_Picture_new;
	Draw_Picture_new *canvasRatio = new Draw_Picture_new;

	TFile *file_read[2];
	file_read[0] = new TFile(Form("%s/file_pid/Flow_%iGeVRun10_PID_10binYesTofdca11_sys.root",outPath, Energy),"READ");;
	file_read[1] = new TFile(Form("%s/file_pid/Flow_%iGeV_PID_10binYesTofdca11_sys.root",outPath, Energy),"READ");;
	TFile *file_out = new TFile("OUT.root","RECREATE");	
	Int_t b=0;
	Int_t b2=0;

	Float_t 	max_y_axis_2=1.23;

	for(Int_t ch=0; ch<2;ch++){
		for(Int_t c=0; c<3;c++){
			for (Int_t f = 0; f < 2; f++){
			    //Data line
				rebin_vec = Rebin2( harmonic, Energy, PID, charge_arr[ch], centBinMin[harmonic-2][c], centBinMax[harmonic-2][c]);			    
			    object.push_back(new Draw_Object);
			    Analysis -> SetParametrs(Energy, file_read[f],file_out,harmonic,PID,charge_arr[ch],EtaGap,prefix);
			    object.back() -> SetParametrsGraph("data",Form("#font[42]{#scale[1.2]{%s}}",data[f].Data()), color[f], style[f], size[f]);
			    Analysis -> FlowVsPt_EtaSub(object.back(), centBinMin[harmonic-2][c], centBinMax[harmonic-2][c], rebin_vec);
			}

			object_ratio.push_back(new Draw_Object);
			Analysis -> RatioGraphPointToPoint(object_ratio.back(),object[0],object[1],"");
			object_ratio.back() -> SetParametrsGraph("data",Form("#font[42]{#scale[1.2]{#frac{Old 27 GeV}{New 27 GeV}}}"), color[2], style[2], size[2]);
			Ratio.push_back(object_ratio);
			object_ratio.clear();

			flow.push_back(object);
			object.clear();
			
			if(harmonic==2) canvasRatio->SetNumberPad(Form("#font[42]{%s}", bukva[b2].Data()), "Pad", 0.04, 0.95*max_y_axis_2, 0.09, 0.0);
			if(harmonic==3) canvasRatio->SetNumberPad(Form("#font[42]{%s}", bukva[b2].Data()), "Pad", 0.04, 0.95*1.43, 0.09, 0.0);
			b2++;

			if(ch==0){
				canvasFlow->SetTextInfo(Form("#font[42]{%i - %i %%}", cent[centBinMax[harmonic-2][c]+1], cent[centBinMin[harmonic-2][c]]), Form("C%i",c), 0.04, 0.73*max_y_axis, 0.09, 0.0);
			}

			canvasFlow->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.85*max_y_axis, 0.09, 0.0);
			b++;

		}

	}

    
    canvasFlow->SetTextInfo(Form("#font[42]{ v_{%i}(%s)}", harmonic, particle[GetNumberParticle(PID,"Pos")].Data() ),"Y",0.5, 0.7,0.3,90);
	canvasFlow->SetTextInfo(Form("#font[42]{ v_{%i}(%s)}", harmonic, particle[GetNumberParticle(PID,"Neg")].Data() ),"Y",0.5, 0.2,0.3,90);

	canvasFlow->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.5,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, #sqrt{s_{NN}}=27 GeV}}"),"Info",0.3,0.85*max_y_axis,0.09,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{STAR Preliminary}}"),"Prel",0.3,0.85*max_y_axis,0.09,0.0);
	canvasFlow->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasFlow->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasFlow->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.086);
	canvasFlow->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.086);
	canvasFlow->SetDrawObject(flow[0],"OnePad");
	canvasFlow->SetLegend(0.3,0.7,0.98,0.95,1);
	TCanvas *result = canvasFlow->CanvasNxM(1200,800, 2, 3, 0.6, 0.25, flow, 2,0,3);

	result->SaveAs(Form("%s/picture/3_PID/New27Gev/ver2_v%i_%s_%s_cent.pdf",outPath, harmonic, PID.Data(), EtaGap.Data() ));
	result->SaveAs(Form("%s/picture/3_PID/New27Gev/ver2_v%i_%s_%s_cent.png",outPath, harmonic, PID.Data(), EtaGap.Data() ));
	delete result;


	canvasRatio->SetTextInfo(Form("#font[42]{v_{%i}(%s)}", harmonic, particle[GetNumberParticle(PID,"Pos")].Data() ),"Y",0.5, 0.7,0.3,90);
	canvasRatio->SetTextInfo(Form("#font[42]{ v_{%i}(%s)}", harmonic, particle[GetNumberParticle(PID,"Neg")].Data() ),"Y",0.5, 0.2,0.3,90);
	canvasRatio->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.5,0.0);
	canvasRatio->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, #sqrt{s_{NN}}=27 GeV}}"),"Info",0.3,0.85*max_y_axis,0.09,0.0);
	canvasRatio->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasRatio->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	if(harmonic==2){
	canvasRatio->SetTextInfo(Form("#font[42]{ #scale[1.0]{ STAR Preliminary}}"),"Prel",0.3,0.96*1.23,0.09,0.0);
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  0.77, 1.23, 0.085);
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  0.77, 1.23, 0.085);
	}else{
	canvasRatio->SetTextInfo(Form("#font[42]{ #scale[1.0]{ STAR Preliminary}}"),"Prel",0.3,0.94*1.43,0.09,0.0);
		canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  0.57, 1.43, 0.086);
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  0.57, 1.43, 0.086);
	}

	canvasRatio->SetDrawObject(Ratio[0],"OnePad");
	canvasRatio->SetLegend(0.2,0.7,0.8,0.98,1);
	TCanvas *result2 = canvasRatio->CanvasNxM(1200,800, 2, 3, 0.6, 0.25, Ratio, 2,0,1);

	result2->SaveAs(Form("%s/picture/3_PID/New27Gev/ver2_Ratio_v%i_%s_%s_cent.pdf",outPath, harmonic, PID.Data(), EtaGap.Data() ));
	result2->SaveAs(Form("%s/picture/3_PID/New27Gev/ver2_Ratio_v%i_%s_%s_cent.png",outPath, harmonic, PID.Data(), EtaGap.Data() ));
	delete result2;


}





void FlowVsPtForBES(Int_t harmonic, Int_t CentBinMin, Int_t CentBinMax, Double_t ptMin, Double_t ptMax, TString EtaGap, TString PID, TString charge, Double_t max_y_axis,Double_t max_y_axis2, Double_t max_x_axis,Double_t MinYratio, Double_t MaxYratio, Double_t MinYratio2, Double_t MaxYratio2){

	VEC rebin_vec;
	VEC rebin_vec_integral;
	TString prefix;
	
	if( strncmp(PID, "Hadrons",7)==0){
		prefix="";
		rebin_vec_integral={0.2,ptMin,ptMax,5.2};
    }else{
    	prefix="TPCandTOF";
		rebin_vec_integral={0.15,ptMin,ptMax,5.2};
    }
	
	Int_t Arr_energy[7]={7,11,14,19,27,39,62};
	const Int_t color[]={2, 1, 4, 6, 4, 2, 46};
	const Int_t style[]={24, 33, 23, 29, 34, 47, 8, 29};
	const Double_t size[]={1.3,1.4,1.3,1.4,1.3,1.3,1.3};

	//Analysis_Function *gr = new Analysis_Function;
	TFile *file_1[7];
	for(Int_t i=0; i<7; i++){
		if( strncmp(PID, "Hadrons",7)==0){
			file_1[i] = new TFile(Form("%s/file_hadrons/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2.root",outPath,Arr_energy[i]),"READ");	
	    }else{
			file_1[i] = new TFile(Form("%s/file/Flow_%iGeV_PID_10binYesTofdca11_oldPID_BadRunOld.root",outPath, Arr_energy[i]),"READ");
	    }
	}
	TFile *file_out = new TFile("OUT.root","RECREATE");	

	Analysis_Function *Analysis = new Analysis_Function;
	
	std::vector<Draw_Object *> flow;
	std::vector<Draw_Object *> ratio;
	std::vector<Draw_Object *> flow_integral;
	std::vector<Draw_Object *> ratio_integral;

	for(Int_t j=6; j>=0; j--){

		rebin_vec = Rebin2( harmonic, Arr_energy[j], PID, charge, CentBinMin, CentBinMax);
	   	
	   	//Data line
	    flow.push_back(new Draw_Object);
	    Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harmonic,PID,charge,EtaGap,prefix);
	    flow.back() -> SetParametrsGraph("data",Form("#font[42]{%s GeV}",Energy_scan_text.at(Arr_energy[j])), color[j], style[j], size[j]);
	    Analysis -> FlowVsPt_EtaSub(flow.back(), CentBinMin, CentBinMax, rebin_vec);
		
	    flow_integral.push_back(new Draw_Object);
	    flow_integral.back() -> SetParametrsGraph("data",Form("#font[42]{%s GeV}",Energy_scan_text.at(Arr_energy[j])), color[j], style[j], size[j]);
	    Analysis -> FlowVsPtDivideIntegralFLow_EtaSub(flow_integral.back(), CentBinMin, CentBinMax, rebin_vec, rebin_vec_integral, 1);

		rebin_vec.clear();
	}
	//ratio line

	for(Int_t i=1; i<(int)flow.size(); i++){
	    ratio.push_back(new Draw_Object);
	    Analysis -> RatioGraphEVAL(ratio[i-1],flow[i],flow[0],"");
	    
	    ratio_integral.push_back(new Draw_Object);
	    Analysis -> RatioGraphEVAL(ratio_integral[i-1],flow_integral[i],flow_integral[0],"");
	
	}

	//
	Draw_Picture_new *canvas = new Draw_Picture_new;
   
    canvas->SetTLine(0., 0., max_x_axis, 0., 1, 2, 2, "Up");
    canvas->SetTLine(0., 0., max_x_axis, 0., 1, 2, 2, "Up");
    canvas->SetTLine(0., 0.9, max_x_axis, 0.9, 1, 2, 2, "Down");
    canvas->SetTLine(0., 1.1, max_x_axis, 1.1, 1, 2, 2, "Down");
    canvas->SetTLine(0., 1., max_x_axis, 1., flow[0] -> GetColor(), 2, 1, "Down");
	canvas->SetAxisToCanvsWithBottomPanel(0., max_x_axis, -0.01, max_y_axis, MinYratio, MaxYratio);
	canvas->SetDrawObject(flow,"Up");
	canvas->SetDrawObject(ratio,"Down");
	canvas->SetLegend(0.1,0.3,0.3,0.89,2);
	TCanvas *result = canvas->CanvasWithBottomPanel(Form("#font[42]{ #scale[1.2]{ #splitline{Au+Au #sqrt{s_{NN}} =7.7-62.4 GeV}{%s, centrality %i-%i %%}}}", particle[GetNumberParticle(PID,charge)].Data(),cent[CentBinMax+1], cent[CentBinMin]),0.3*max_x_axis,0.8*max_y_axis,
													"#font[42]{p_{T} [GeV/c]}",Form("v_{%i}",harmonic),Form("#font[42]{#frac{ v_{%i}(Energy)}{ v_{%i}(62.4 GeV)}}", harmonic, harmonic),"");
	result->SaveAs(Form("%s/picture/energy/v%i_%s%s_%i_%i_BES_%s.png",outPath, harmonic, PID.Data(), charge.Data(), cent[CentBinMax+1], cent[CentBinMin], EtaGap.Data()));
	delete result;

	// Integral
	Draw_Picture_new *canvas2 = new Draw_Picture_new;

    canvas2->SetTLine(0., 0., max_x_axis, 0., 1, 2, 2, "Up");
    canvas2->SetTLine(0., 0., max_x_axis, 0., 1, 2, 2, "Up");
    canvas2->SetTLine(0., 0.9, max_x_axis, 0.9, 1, 2, 2, "Down");
    canvas2->SetTLine(0., 1.1, max_x_axis, 1.1, 1, 2, 2, "Down");
    canvas2->SetTLine(0., 1., max_x_axis, 1., flow[0] -> GetColor(), 2, 1, "Down");
	canvas2->SetAxisToCanvsWithBottomPanel(0., max_x_axis, -0.01, max_y_axis2, MinYratio2, MaxYratio2);
	canvas2->SetDrawObject(flow_integral,"Up");
	canvas2->SetDrawObject(ratio_integral,"Down");
	canvas2->SetLegend(0.1,0.3,0.3,0.89,2);
	result = canvas2->CanvasWithBottomPanel(Form("#font[42]{ #scale[1.2]{ #splitline{Au+Au #sqrt{s_{NN}} =7.7-62.4 GeV}{%s, %.1f<p_{T}<%.1f,  %i-%i %%}}}", particle[GetNumberParticle(PID,charge)].Data(), ptMin,ptMax,cent[CentBinMax+1], cent[CentBinMin]),0.3*max_x_axis,0.8*max_y_axis2,
													"#font[42]{p_{T} [GeV/c]}",Form("#frac{v_{%i}}{v_{%i}^{int}}",harmonic, harmonic),Form("#font[42]{#frac{ Obs(Energy)}{ Obs(62.4 GeV)}}"),"");
	result->SaveAs(Form("%s/picture/energy/aIntegral_v%i_%s%s_%i_%i_BES_%s.png",outPath, harmonic, PID.Data(), charge.Data(), cent[CentBinMax+1], cent[CentBinMin], EtaGap.Data()));
	delete result;

	for(Int_t i=0; i<7; i++){
		file_1[i]->Close();
		delete file_1[i];
	}
}

void FlowVsRapidity(Int_t Energy, Int_t harmonic, Int_t CentBinMin, Int_t CentBinMax, Double_t ptMin, Double_t ptMax, TString EtaGap, TString PID, TString charge,  Double_t min_y_axis, Double_t max_y_axis){

	VEC rebin_vec;
	if(harmonic==2){
		rebin_vec = {-1.0,-0.75,-0.6,-0.45,-0.3,-0.15,0,0.15,0.3,0.45,0.6,0.75,1.0};
		//rebin_vec = {-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.15,0,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
	}else{
		rebin_vec = {-1.0,-0.8,-0.6,-0.4,-0.15,0,0.15,0.4,0.6,0.8,1.0};
	}
	
	TString prefix;
	TString nameAxisX;
	if( strncmp(PID, "Hadrons",7)==0){
		prefix="";
		nameAxisX="#font[42]{#eta}";
    }else{
    	prefix="TPCandTOF";
    	nameAxisX="#font[42]{y}";
    }
	
	//const Int_t color[]={2, 1, 2, 2, 1, 2, 46};
	//const Int_t style[]={26, 22, 25, 21, 24, 20, 8, 29};
	//const Double_t size[]={1.4,1.4,1.4,1.4,1.3,1.3,1.3};

	const Int_t color[]={1, 2, 4};
	const Int_t style_par[]={34, 21, 22};
	const Int_t style_ant[]={28, 25, 26};
	const Double_t size[]={2.0,2.0,2.0};
	Int_t cent_bin_Min[3]={7,2,2};
	Int_t cent_bin_Max[3]={8,6,8};

	TFile *file_read;
	if( strncmp(PID, "Hadrons",7)==0){
		file_read = new TFile(Form("%s/file_hadrons/with_rapidity/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2_new.root",outPath,Energy),"READ");	
    }else{
		file_read = new TFile(Form("%s/file_hadrons/with_rapidity/Flow_%iGeV_PID_10binYesTofdca11_with_rapidity.root",outPath, Energy),"READ");
    }
	
	TFile *file_out = new TFile("OUT.root","RECREATE");	

	std::vector<Draw_Object *> flow;

	Analysis_Function *Analysis = new Analysis_Function;
	
	Int_t p=0;
	//Data line
	
	for(Int_t c=0; c<3; c++){
		if( strncmp(PID, "Hadrons",7)!=0){
			Analysis -> SetParametrs(Energy, file_read,file_out,harmonic,PID,"Pos",EtaGap,prefix);
		    flow.push_back(new Draw_Object);
		    flow.back() -> SetParametrsGraph("data",Form("#font[42]{ #scale[2.0]{%s ,%i-%i %%}}",  particle[GetNumberParticle(PID,"Pos")].Data(),cent[cent_bin_Max[c]+1], cent[cent_bin_Min[c]]), color[p], style_par[p], size[p]);
		    Analysis -> FlowVsEta_EtaSub(flow.back(), cent_bin_Min[c], cent_bin_Max[c], ptMin, ptMax, rebin_vec);
		    
		    DrawEtaPID(file_read,cent_bin_Min[c], cent_bin_Max[c], ptMin, ptMax, Energy, PID, "Pos", EtaGap);

		    Analysis -> SetParametrs(Energy, file_read,file_out,harmonic,PID,"Neg",EtaGap,prefix);
		    flow.push_back(new Draw_Object);
		    flow.back() -> SetParametrsGraph("data",Form("#font[42]{ #scale[2.0]{%s ,%i-%i %%}}",  particle[GetNumberParticle(PID,"Neg")].Data(),cent[cent_bin_Max[c]+1], cent[cent_bin_Min[c]]), color[p], style_ant[p], size[p]);
		    Analysis -> FlowVsEta_EtaSub(flow.back(), cent_bin_Min[c], cent_bin_Max[c], ptMin, ptMax, rebin_vec);
		    
		    DrawEtaPID(file_read,cent_bin_Min[c], cent_bin_Max[c], ptMin, ptMax, Energy, PID, "Neg", EtaGap);

		    p++;
		
		}else{
			Analysis -> SetParametrs(Energy, file_read,file_out,harmonic,PID,"",EtaGap,prefix);
		    flow.push_back(new Draw_Object);
		    flow.back() -> SetParametrsGraph("data",Form("#font[42]{ #scale[2.0]{%s ,%i-%i %%}}",  particle[GetNumberParticle(PID,"")].Data(),cent[cent_bin_Max[c]+1], cent[cent_bin_Min[c]]), color[p], style_par[p], size[p]);
		    Analysis -> FlowVsEta_EtaSub(flow.back(), cent_bin_Min[c], cent_bin_Max[c], ptMin, ptMax, rebin_vec);
		    //DrawEta(file_read, cent_bin_Min[c], cent_bin_Max[c], Energy);
		}
	}

	Draw_Picture_new *canvas = new Draw_Picture_new;
	canvas->SetAxisToCanvsOne( -1.05, 1.05, min_y_axis, max_y_axis);
	canvas->SetDrawObject(flow,"Up");
	canvas->SetLegend(0.5,0.7,0.89,0.9,1);
	TCanvas *result = canvas->CanvasOne(Form("#font[42]{ #scale[0.7]{ #splitline{Au+Au #sqrt{s_{NN}} =%s GeV}{ %.1f<p_{T}<%.1f GeV/c }}}", Energy_scan_text.at(Energy),ptMin, ptMax),-0.98*1.05,0.86*max_y_axis,
											 nameAxisX,Form("v_{%i}",harmonic),"");
	if(strncmp(PID, "Hadrons",7)==0){
		result->SaveAs(Form("%s/picture/2_Charge_Hadrons/%i_GeV_v%i_eta_%s%s_%i_%i_%s.png",outPath, Energy, harmonic, PID.Data(), charge.Data(), cent[CentBinMax+1], cent[CentBinMin], EtaGap.Data()));
	}else{
		result->SaveAs(Form("%s/picture/2_Charge_Hadrons/pid/%i_GeV_v%i_eta_%s%s_%i_%i_%s.png",outPath, Energy, harmonic, PID.Data(), charge.Data(), cent[CentBinMax+1], cent[CentBinMin], EtaGap.Data()));
	}
	delete result;

	file_read->Close();

}

void FlowVsRapiditySymmetry(Int_t Energy, Int_t harmonic, Double_t ptMin, Double_t ptMax, TString EtaGap, TString PID, TString charge,  Double_t min_y_axis, Double_t max_y_axis){

	VEC rebin_vec;
	if(harmonic==2){
		//rebin_vec = {-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.15,0,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
		rebin_vec = {-1.0,-0.75,-0.6,-0.45,-0.3,-0.15,0,0.15,0.3,0.45,0.6,0.75,1.0};
	}else{
		rebin_vec = {-1.0,-0.7,-0.4,-0.15,0,0.15,0.4,0.7,1.0};
	}
	
	TString prefix;
	TString nameAxisX;
	if( strncmp(PID, "Hadrons",7)==0){
		prefix="";
		nameAxisX="#font[42]{#eta}";
    }else{
    	prefix="TPCandTOF";
    	nameAxisX="#font[42]{y}";
    }
	
	const Int_t color[]={1, 2, 4};
	const Int_t style[]={34, 21, 22};
	const Int_t style_sym[]={28, 25, 26};
	const Double_t size[]={1.4,1.4,1.4};
	Int_t cent_bin_Min[3]={7,2,2};
	Int_t cent_bin_Max[3]={8,6,8};

	TFile *file_read;
	if( strncmp(PID, "Hadrons",7)==0){
		file_read = new TFile(Form("%s/file_hadrons/with_rapidity/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2_new.root",outPath,Energy),"READ");	
    }else{
		file_read = new TFile(Form("%s/file/Flow_%iGeV_PID_10binYesTofdca11_oldPID_BadRunOld.root",outPath, Energy),"READ");
    }
	
	TFile *file_out = new TFile("OUT.root","RECREATE");	
	
	Analysis_Function *Analysis = new Analysis_Function;
	Analysis -> SetParametrs(Energy, file_read,file_out,harmonic,PID,charge,EtaGap,prefix);
	
	std::vector<Draw_Object *> flow;
	std::vector<Draw_Object *> ratio;

	for(Int_t c=0; c<3; c++){
	    flow.push_back(new Draw_Object);
   		flow.back() -> SetParametrsGraph("data",Form("#font[42]{%i-%i %%}",cent[cent_bin_Max[c]+1], cent[cent_bin_Min[c]]), color[c], style[c], size[c]);
    	Analysis -> FlowVsEta_EtaSub(flow.back(), cent_bin_Min[c], cent_bin_Max[c], ptMin, ptMax, rebin_vec);

    	flow.push_back(new Draw_Object);
   		flow.back() -> SetParametrsGraph("data",Form("#font[42]{%i-%i %%,symm}",cent[cent_bin_Max[c]+1], cent[cent_bin_Min[c]]), color[c], style_sym[c], size[c]);
    	Analysis -> FlowVsEta_EtaSub(flow.back(), cent_bin_Min[c], cent_bin_Max[c], ptMin, ptMax, rebin_vec);
    	for(Int_t p=0; p < flow.back()->GetSizeVector(); p++){
    		flow.back()->СhangePointGraph(p,(-1)*flow.back()->GetPointXGraph(p), flow.back()->GetPointYGraph(p), 0, flow.back()->GetPointYErrorGraph(p), 0.);
    	}
    	//DrawEta(file_read,cent_bin_Min[c], cent_bin_Max[c], Energy);
	}
	ratio.push_back(new Draw_Object);
   	Analysis -> RatioGraphEVAL(ratio.back(),flow[1],flow[0],"");
   	ratio.push_back(new Draw_Object);
   	Analysis -> RatioGraphEVAL(ratio.back(),flow[3],flow[2],"");
   	ratio.push_back(new Draw_Object);
   	Analysis -> RatioGraphEVAL(ratio.back(),flow[5],flow[4],"");
	/*
	canvas->SetDrawObject(flow,"Up");
	canvas->SetDrawObject(ratio,"Down");
	canvas->SetLegend(0.1,0.55,0.4,0.88);
	TCanvas *result = canvas->CanvasWithBottomPanel(Form("#font[42]{ #scale[1.2]{ #splitline{Au+Au #sqrt{s_{NN}} = %s GeV}{%s, centrality %i-%i %%}}}",Energy_scan_text.at(Energy), particle[GetNumberParticle(PID,charge)].Data(),cent[CentBinMax+1], cent[CentBinMin]),0.36*max_x,0.8*max_y,
													"#font[42]{p_{T} [GeV/c]}","v_{2}","#font[42]{#frac{Data}{Published}}",PhysRef(PID));
	result->SaveAs(Form("%s/picture/1_comparison_with_article/%iGeV/Publ_pt_%s%s_%i_%i_%iGeV_%s.png",outPath, Energy,PID.Data(), charge.Data(), cent[CentBinMax+1], cent[CentBinMin], Energy, EtaGap.Data()));
	delete result;
	*/
   	//Data line

	Draw_Picture_new *canvas = new Draw_Picture_new;
	//canvas->SetAxisToCanvsOne( -1.05, 1.05, min_y_axis, max_y_axis);
	canvas->SetTLine(-1.05, 0.9, 1.05, 0.9, 1, 2, 2, "Down");
    canvas->SetTLine(-1.05, 1.1, 1.05, 1.1, 1, 2, 2, "Down");
    canvas->SetTLine(-1.05, 1., 1.05, 1., 1, 2, 1, "Down");
	canvas->SetAxisToCanvsWithBottomPanel(-1.05, 1.05,  min_y_axis, max_y_axis, 0.8, 1.2);
	canvas->SetDrawObject(flow,"Up");
	canvas->SetDrawObject(ratio,"Down");
	canvas->SetLegend(0.6,0.65,0.89,0.89,1);
	TCanvas *result = canvas->CanvasWithBottomPanel(Form("#font[42]{ #scale[1.0]{ #splitline{Au+Au #sqrt{s_{NN}} =%s GeV, %s}{ %.1f<p_{T}<%.1f GeV/c }}}", Energy_scan_text.at(Energy), particle[GetNumberParticle(PID,charge)].Data(),ptMin, ptMax),-0.98*1.05,0.86*max_y_axis,
											 nameAxisX,Form("v_{%i}",harmonic),"#font[42]{#frac{Symmetrically}{Data}}","");
	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/Symm%i_GeV_v%i_eta_%s%s_%s_%.2f_%.2f.png",outPath, Energy, harmonic, PID.Data(), charge.Data(), EtaGap.Data(), ptMin, ptMax));
	delete result;

	file_read->Close();

}

void FlowVsRapidityForBES(Int_t harmonic, Double_t ptMin, Double_t ptMax, TString EtaGap, TString PID, TString charge,  Double_t min_y_axis, Double_t max_y_axis){

	VEC rebin_vec;
	if(harmonic==2){
		//rebin_vec = {-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.15,0,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
		rebin_vec = {-1.0,-0.75,-0.6,-0.45,-0.3,-0.15,0,0.15,0.3,0.45,0.6,0.75,1.0};
	}else{
		rebin_vec = {-1.0,-0.7,-0.4,-0.15,0,0.15,0.4,0.7,1.0};
	}
	
	TString prefix;
	TString nameAxisX;
	if( strncmp(PID, "Hadrons",7)==0){
		prefix="";
		nameAxisX="#font[42]{#eta}";
    }else{
    	prefix="TPCandTOF";
    	nameAxisX="#font[42]{y}";
    }
	
	Int_t Arr_energy[7]={7,11,14,19,27,39,62};
	const Int_t color[]={1, 4, 1, 1, 1, 4, 2};
	const Int_t style[]={26, 22, 27, 23, 25, 20, 21, 29};
	const Double_t size[]={1.4,1.4,1.6,1.4,1.4,1.4,1.4};
	Int_t cent_bin_Min[3]={7,2,2};
	Int_t cent_bin_Max[3]={8,6,8};

	TFile *file_read[7];
	for(Int_t i=0; i<7; i++){
		if( strncmp(PID, "Hadrons",7)==0){
			file_read[i] = new TFile(Form("%s/file_hadrons/with_rapidity/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2_new.root",outPath,Arr_energy[i]),"READ");	
	    }else{
			file_read[i] = new TFile(Form("%s/file/Flow_%iGeV_PID_10binYesTofdca11_oldPID_BadRunOld.root",outPath, Arr_energy[i]),"READ");
	    }
	}
	TFile *file_out = new TFile("OUT.root","RECREATE");	

	Analysis_Function *Analysis = new Analysis_Function;
	
	std::vector<std::vector<Draw_Object *>> flow_draw;
	std::vector<std::vector<Draw_Object *>> ratio_draw;
	std::vector<TString> TitlePad;

	std::vector<Draw_Object *> flow;
	std::vector<Draw_Object *> ratio;

	for(Int_t c=0; c<3; c++){
		for(Int_t g=0; g<7; g++){
			Analysis -> SetParametrs(Arr_energy[g], file_read[g],file_out,harmonic,PID,charge,EtaGap,prefix);
		    flow.push_back(new Draw_Object);
	   		flow.back() -> SetParametrsGraph("data",Form("#scale[1.2]{#font[42]{%s GeV}}",Energy_scan_text.at(Arr_energy[g])), color[g], style[g], size[g]);
	    	Analysis -> FlowVsEta_EtaSub(flow.back(), cent_bin_Min[c], cent_bin_Max[c], ptMin, ptMax, rebin_vec);
	    }

	    for(Int_t g=0; g<7; g++){
	    	    ratio.push_back(new Draw_Object);
	   			Analysis -> RatioGraphEVAL(ratio[g],flow[g],flow[5],"");
	    }

	   	flow_draw.push_back(flow);
	   	ratio_draw.push_back(ratio);
	   	flow.clear();
	   	ratio.clear();
	   	TitlePad.push_back(Form("#font[42]{%i-%i %%}",cent[cent_bin_Max[c]+1], cent[cent_bin_Min[c]]));

	}
	
	Draw_Picture_new *canvas = new Draw_Picture_new;
	
    canvas->SetTLine(-1.05, 1., 1.05, 1., flow_draw[0][0] -> GetColor(), 2, 2, "Down");
	canvas->SetAxisToCanvsNM(-1.05, 1.05,  min_y_axis, max_y_axis, 0.07);
	canvas->SetDrawObject(flow_draw[0],"OnePad");
	canvas->SetLegend(0.02,0.02,0.97,0.3,1);
	TCanvas *result = canvas->CanvasNxMRatio( 3, 0.35, 0.25, flow_draw, ratio_draw, Form("#font[42]{ #scale[1.0]{ #splitline{Au+Au, %s}{ %.1f<p_{T}<%.1f GeV/c }}}", particle[GetNumberParticle(PID,charge)].Data(),ptMin, ptMax),-0.98*1.05,0.86*max_y_axis,
											TitlePad, 0.9*1.05,0.86*max_y_axis,
											nameAxisX,Form("v_{%i}",harmonic),Form("#font[42]{#frac{v_{%i}(Energy)}{v_{%i}(39 GeV)}}",harmonic, harmonic),"");

	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/BES_v%i_eta_%s%s_%s_%.2f_%.2f.png",outPath, harmonic, PID.Data(), charge.Data(), EtaGap.Data(), ptMin, ptMax));
	/*delete result;

	for(Int_t i=0; i<7; i++){
		file_read[i]->Close();
		delete file_read[i];
	}
	*/
}


void FlowVsEnergyBESupdate(Double_t ptMin, Double_t ptMax, TString EtaGap, TString PID, TString charge, Double_t min_y_axis, Double_t max_y_axis, Double_t min_x_axis, Double_t max_x_axis,Double_t MinYratio, Double_t MaxYratio){

	VEC rebin_vec = {0.15,ptMin,ptMax,5.2};;
	Int_t point = 1; 
	TString prefix="";
	
	Int_t Arr_energy[6]={11,14,19,27,39,62};
	Double_t x_energy[6]={11.5,14.5,19.6,27,39,62.4};
	Int_t cent_bin_Min[9]={8,7,6,5,4,2,2,1,0};
	Int_t cent_bin_Max[9]={8,7,6,5,4,3,2,1,0};

	TFile *file_1[6];
	for(Int_t i=0; i<6; i++){
		if( strncmp(PID, "Hadrons",7)==0){
			file_1[i] = new TFile(Form("%s/file_hadrons/Sys/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2_sys.root",outPath,Arr_energy[i]),"READ");	
	    }else{
			file_1[i] = new TFile(Form("%s/file/Flow_%iGeV_PID_10binYesTofdca11_sys.root",outPath, Arr_energy[i]),"READ");
	    }
	}
	TFile *file_out = new TFile("OUT.root","RECREATE");	

	Analysis_Function *Analysis = new Analysis_Function;
	Draw_Object *object = new Draw_Object;
	std::vector<std::vector<Draw_Object *>> flow ;
	std::vector<std::vector<Draw_Object *>> ratio ;
	std::vector<Double_t> flow_sys_err;
	std::vector<std::vector<std::vector<Double_t>>> ratio_sys_err3;
	std::vector<std::vector<Double_t>> ratio_sys_err2;

	Double_t x,y,ex,ey;
	Double_t den,den_err;
	std::vector<TString> TitlePad;

	Draw_Picture_new *canvasFlow = new Draw_Picture_new;
	Draw_Picture_new *canvasRatio = new Draw_Picture_new;
	Int_t b=0;
	Int_t b2=0;

	std::vector<TString> sys_arr={"Eta15","Eta03","Eta05","Eta07"};

	for(Int_t c=0; c<6; c++){
		flow.push_back({new Draw_Object, new Draw_Object});
		ratio.push_back({new Draw_Object, new Draw_Object});
    	std::cout<<"good1\n";
		for(Int_t harm=0; harm<2; harm++){
			for(Int_t j=0; j<6; j++){ 	
			   	//Data line
			    Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,PID,charge,EtaGap,prefix);
			    Analysis -> FlowVsPt_EtaSub(object, cent_bin_Min[c], cent_bin_Max[c], rebin_vec);
			    object->GetPointDrawObject(point,x,y,ex,ey);
			    flow[c][harm]->SetPointGraph(x_energy[j],y,0.,ey);
			    ratio[c][harm]->SetPointGraph(x_energy[j],y,0.,ey);
			    object->ClearObject();
			}

			////////
		    for(Int_t s=0; s<sys_arr.size();s++){
				for(Int_t j=0; j<6; j++){ 	
				   	//Data line
				    Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,PID,charge,sys_arr[s],prefix);
				    Analysis -> FlowVsPt_EtaSub(object, cent_bin_Min[c], cent_bin_Max[c], rebin_vec);
				    object->GetPointDrawObject(point,x,y,ex,ey);
				    flow_sys_err.push_back(y);
				    object->ClearObject();
				}
				ratio_sys_err2.push_back(flow_sys_err);
				flow[c][harm]->SetSysArrey(flow_sys_err);
				flow_sys_err.clear();
		    }
		    ratio_sys_err3.push_back(ratio_sys_err2);
		    flow[c][harm]->SetSysErrors();
			///////

		}
		
		canvasFlow->SetNumberPad(Form("#font[42]{#scale[0.85]{%s %i-%i%%}}", bukva[b].Data(), cent[cent_bin_Max[c]+1], cent[cent_bin_Min[c]]), "Pad", 8, 0.85*max_y_axis, 0.09, 0.0);
		b++;

		canvasRatio->SetNumberPad(Form("#font[42]{#scale[0.85]{%s %i-%i%%}}", bukva[b].Data(), cent[cent_bin_Max[c]+1], cent[cent_bin_Min[c]]), "Pad", 8, 0.92*MaxYratio, 0.09, 0.0);
		b2++;

		for(Int_t harm=0; harm<2; harm++){
			flow[c][harm]->GetPointDrawObject(5,x,den,ex,den_err);
			Analysis->RatioGraphToOnePoint(ratio[c][harm],den, den_err);

			///////
		    for(Int_t s=0; s<sys_arr.size();s++){
		    	for(Int_t p=0; p<ratio_sys_err3[harm][s].size(); p++){
		    		Double_t den_sys = ratio_sys_err3[harm][s][5];
					ratio_sys_err3[harm][s][p] = ratio_sys_err3[harm][s][p] / den_sys;
				}
				ratio[c][harm]->SetSysArrey(ratio_sys_err3[harm][s]);
			}
			ratio[c][harm]->SetSysErrors();
			//////


		}
		
		flow[c][0] -> SetParametrsGraph(Form("data_v2_Hadrons_int_cent_%i_%i",cent[cent_bin_Max[c]+1], cent[cent_bin_Min[c]]),"#font[42]{v_{2}}",4,23,1.8);
		ratio[c][0] -> SetParametrsGraph(Form("ratio_data_v2_Hadrons_int_cent_%i_%i",cent[cent_bin_Max[c]+1], cent[cent_bin_Min[c]]),"#font[42]{v_{2}}",4,23,1.8);
		
		flow[c][1] -> SetParametrsGraph(Form("data_v3_Hadrons_int_cent_%i_%i",cent[cent_bin_Max[c]+1], cent[cent_bin_Min[c]]),"#font[42]{v_{3}}",2,22,1.8);
		ratio[c][1] -> SetParametrsGraph(Form("ratio_data_v3_Hadrons_int_cent_%i_%i",cent[cent_bin_Max[c]+1], cent[cent_bin_Min[c]]),"#font[42]{v_{3}}",2,22,1.8);
	}
	//

	canvasFlow->SetTextInfo(Form("#font[42]{#scale[1.7]{v_{n}}}"),"Y",0.5, 0.7,0.3,0);
	canvasFlow->SetTextInfo("#font[42]{#scale[0.85]{#sqrt{s_{NN}} [GeV]}}","X",0.5,0.4,0.5,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[0.85]{ #splitline{Au+Au, %s}{ %.1f<p_{T}<%.1f GeV/c }}}", particle[GetNumberParticle(PID,charge)].Data(),ptMin, ptMax),"Info",0.35*(max_x_axis - min_x_axis),0.85*max_y_axis,0.09,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[0.85]{ STAR Preliminary}}"),"Prel",0.43*(max_x_axis - min_x_axis),0.85*max_y_axis,0.09,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[0.85]{#splitline{Not corrected for p_{T}}{dependent efficiency}}}"),"2Prel",0.1*(max_x_axis - min_x_axis),0.6*max_y_axis,0.09,0.0);
	canvasFlow->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasFlow->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasFlow->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.086);
	canvasFlow->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.086);
	canvasFlow->SetDrawObject(flow[0],"OnePad");
	canvasFlow->SetLegend(0.5,0.7,0.98,0.95,1);
	TCanvas *result = canvasFlow->CanvasNxM(1000,1280, 3, 2, 0.44, 0.35, flow, 2,0,1);

	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/sys/vn_integral_%s_%s_%i_%i_cent.pdf",outPath, PID.Data(), EtaGap.Data(),cent[cent_bin_Max[5]+1], cent[cent_bin_Min[5]] ));
	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/sys/vn_integral_%s_%s_%i_%i_cent.png",outPath, PID.Data(), EtaGap.Data(),cent[cent_bin_Max[5]+1], cent[cent_bin_Min[5]] ));
	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/sys/vn_integral_%s_%s_%i_%i_cent.C",outPath, PID.Data(), EtaGap.Data(),cent[cent_bin_Max[5]+1], cent[cent_bin_Min[5]] ));
	delete result;

	canvasRatio->SetTextInfo(Form("#font[42]{#scale[1.2]{#frac{ v_{n}}{ v_{n}(62.4 GeV)}}}"),"Y",0.5, 0.4,0.3,90);
	canvasRatio->SetTextInfo("#font[42]{#scale[0.85]{#sqrt{s_{NN}} [GeV]}}","X",0.5,0.4,0.5,0.0);
	canvasRatio->SetTextInfo(Form("#font[42]{ #scale[0.85]{ #splitline{Au+Au, %s}{ %.1f<p_{T}<%.1f GeV/c }}}", particle[GetNumberParticle(PID,charge)].Data(),ptMin, ptMax),"Info",0.35*(max_x_axis - min_x_axis),MinYratio + 0.15*(MaxYratio - MinYratio),0.09,0.0);
	canvasRatio->SetTextInfo(Form("#font[42]{ #scale[0.85]{ STAR Preliminary}}"),"Prel",0.43*(max_x_axis - min_x_axis),MinYratio + 0.15*(MaxYratio - MinYratio),0.09,0.0);
	canvasRatio->SetTextInfo(Form("#font[42]{ #scale[0.85]{ Not corrected for pT dependent efficiency}}"),"2Prel",0.43*(max_x_axis - min_x_axis),MinYratio + 0.15*(MaxYratio - MinYratio),0.09,0.0);
	canvasRatio->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasRatio->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  MinYratio, MaxYratio, 0.086);
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  MinYratio, MaxYratio, 0.086);
	canvasRatio->SetDrawObject(ratio[0],"OnePad");
	canvasRatio->SetLegend(0.5,0.1,0.98,0.4,1);
	TCanvas *result2 = canvasRatio->CanvasNxM(1000,1280, 3, 2, 0.44, 0.35, ratio, 2,0,1);

	result2->SaveAs(Form("%s/picture/2_Charge_Hadrons/sys/vn_integral_ratio_%s_%s_%i_%i_cent.pdf",outPath, PID.Data(), EtaGap.Data(),cent[cent_bin_Max[5]+1], cent[cent_bin_Min[5]] ));
	result2->SaveAs(Form("%s/picture/2_Charge_Hadrons/sys/vn_integral_ratio_%s_%s_%i_%i_cent.png",outPath, PID.Data(), EtaGap.Data(),cent[cent_bin_Max[5]+1], cent[cent_bin_Min[5]] ));
	result2->SaveAs(Form("%s/picture/2_Charge_Hadrons/sys/vn_integral_ratio_%s_%s_%i_%i_cent.C",outPath, PID.Data(), EtaGap.Data(),cent[cent_bin_Max[5]+1], cent[cent_bin_Min[5]] ));
	delete result2;


  /* 	
   	canvas->SetTLine(0., 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvas->SetAxisToCanvsNM(0., max_x_axis,  min_y_axis, max_y_axis, 0.07);
	canvas->SetDrawObject(flow[0],"OnePad");
	canvas->SetLegend(0.1,0.5,0.93,0.7);
	TCanvas *result = canvas->CanvasNxM(3, 3, 0.35, 0.35, flow, Form("#font[42]{ #scale[1.5]{ #splitline{Au+Au, %s}{ %.1f<p_{T}<%.1f GeV/c }}}", particle[GetNumberParticle(PID,charge)].Data(),ptMin, ptMax), 5,0.8*max_y_axis,
											TitlePad, 60, 0.86*max_y_axis,
											"#font[42]{#scale[0.65]{#sqrt{s_{NN}[GeV]}}}",Form("#scale[1.7]{v_{n}}"),"",0);

	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/vn_pt_onePad_%.1f_%.1f_%s%s_BES_%s.png",outPath, ptMin, ptMax,PID.Data(), charge.Data(), EtaGap.Data()));
	delete result;

*/
	/*
	Draw_Picture_new *canvas2 = new Draw_Picture_new;
   	Float_t max_ratio_y =1.1;
   	canvas2->SetTLine(0., 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvas2->SetAxisToCanvsNM(0., max_x_axis,  0.5 , max_ratio_y, 0.07);
	canvas2->SetDrawObject(flow[0],"OnePad");
	canvas2->SetLegend(0.1,0.5,0.93,0.7);
	TCanvas *result2 = canvas2->CanvasNxM(3, 3, 0.35, 0.35, ratio, Form("#font[42]{ #scale[1.5]{ #splitline{Au+Au, %s}{ %.1f<p_{T}<%.1f GeV/c }}}", particle[GetNumberParticle(PID,charge)].Data(),ptMin, ptMax), 5,0.92*max_ratio_y,
											TitlePad, 60, 0.86*max_ratio_y,
											"#font[42]{#scale[0.65]{#sqrt{s_{NN}[GeV]}}}",Form("#font[42]{#scale[1.0]{#frac{ v_{n}}{ v_{n}(62.4 GeV)}}}"),"",1);

	result2->SaveAs(Form("%s/picture/2_Charge_Hadrons/vn_pt_onePad_ratio_%.1f_%.1f_%s%s_BES_%s.png",outPath, ptMin, ptMax,PID.Data(), charge.Data(), EtaGap.Data()));
	

	delete result;
	*/

	for(Int_t i=0; i<6; i++){
		file_1[i]->Close();
		delete file_1[i];
	}
}

void FlowDifferentEtaGap(Int_t Energy, Int_t harmonic, Int_t CentBinMin, Int_t CentBinMax, TString PID, TString charge, Double_t max_y_axis, Double_t max_x_axis, Double_t MinYratio, Double_t MaxYratio){

	VEC	rebin_vec = Rebin2( harmonic, Energy, PID, charge, CentBinMin, CentBinMax);
	TString prefix;
	
	if( strncmp(PID, "Hadrons",7)==0){
		prefix="";
    }else{
    	prefix="TPCandTOF";
    }
	
	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<Draw_Object *> flow;
	std::vector<Draw_Object *> ratio;
	
	std::vector<TString> NameEtaHadrons = {"Eta15","Eta04","Eta07","Eta10"};
	std::vector<Double_t> EtaVecHadrons = {0.075, 0.2, 0.35, 0.5};
	std::vector<TString> NameEtaPID = {"Eta01","Eta03","Eta05","Eta07"};
	std::vector<Double_t> EtaVecPID = {0.05, 0.15, 0.25, 0.35};
	TString suf = "def";

	const Int_t color[]={2, 1, 4, 6, 4, 2, 46};
	const Int_t style[]={23, 22, 25, 23, 34, 47, 8, 29};
	const Double_t size[]={1.5,1.5,1.5,1.5,1.5,1.5,1.5};

	TFile *file_read;

	if( strncmp(PID, "Hadrons",7)==0){
		file_read = new TFile(Form("%s/file_hadrons/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2.root",outPath,Energy),"READ");	
    }else{
    	if(Energy==27){
			file_read = new TFile(Form("%s/file_pid/Flow_%iGeV_PID_10binYesTofdca11_sys.root",outPath, Energy),"READ");
    	}else{
			file_read = new TFile(Form("%s/file_pid/Flow_%iGeV_PID_10binYesTofdca11_sys.root",outPath, Energy),"READ");
    	}
    }

	TFile *file_out = new TFile("OUT.root","RECREATE");

	for(Int_t eta=0; eta<4; eta++){
	    flow.push_back(new Draw_Object);
    	if(strncmp(PID, "Hadrons",7)==0){
    		Analysis -> SetParametrs(Energy, file_read,file_out,harmonic,PID,charge,NameEtaHadrons[eta],prefix);
    		flow.back() -> SetParametrsGraph(Form("data%.3f",EtaVecHadrons[eta]),Form("#font[42]{#eta-gap=%.1f %s}",2.0*EtaVecHadrons[eta], suf.Data()), color[eta], style[eta], size[eta]);
    	}else{
    		Analysis -> SetParametrs(Energy, file_read,file_out,harmonic,PID,charge,NameEtaPID[eta],prefix);
    		flow.back() -> SetParametrsGraph(Form("data%.3f",EtaVecPID[eta]),Form("#font[42]{#eta-gap=%.1f %s}",2.0*EtaVecPID[eta], suf.Data()), color[eta], style[eta], size[eta]);
    	}
    	Analysis -> FlowVsPt_EtaSub(flow.back(), CentBinMin, CentBinMax, rebin_vec);
    	suf="";

	}
    //Data line

    //ratio line
    for(Int_t i=1; i<flow.size();i++){
    	ratio.push_back(new Draw_Object);
    	Analysis -> RatioGraphEVAL(ratio[i-1],flow[i],flow[0],"");
	}
    //canvas setening
    Draw_Picture_new *canvas = new Draw_Picture_new;
    
    canvas->SetTLine(0., 0., max_x_axis, 0., 1, 2, 2, "Up");
    canvas->SetTLine(0., 0., max_x_axis, 0., 1, 2, 2, "Up");
    canvas->SetTLine(0., 0.9, max_x_axis, 0.9, 1, 2, 2, "Down");
    canvas->SetTLine(0., 1.1, max_x_axis, 1.1, 1, 2, 2, "Down");
    canvas->SetTLine(0., 1., max_x_axis, 1., flow[0] -> GetColor(), 2, 1, "Down");
	canvas->SetAxisToCanvsWithBottomPanel(0., max_x_axis, -0.01, max_y_axis, MinYratio, MaxYratio);
	canvas->SetDrawObject(flow,"Up");
	canvas->SetDrawObject(ratio,"Down");
	canvas->SetLegend(0.1,0.55,0.4,0.88,1);
	TCanvas *result = canvas->CanvasWithBottomPanel(Form("#font[42]{ #scale[1.2]{ #splitline{Au+Au #sqrt{s_{NN}} = %s GeV}{%s, centrality %i-%i %%}}}",Energy_scan_text.at(Energy), particle[GetNumberParticle(PID,charge)].Data(),cent[CentBinMax+1], cent[CentBinMin]),0.36*max_x_axis,0.8*max_y_axis,
													"#font[42]{p_{T} [GeV/c]}",Form("v_{%i}",harmonic),"#font[42]{#frac{Data}{Def}}","");
	
	if( strncmp(PID, "Hadrons",7)==0){
		result->SaveAs(Form("%s/picture/2_Charge_Hadrons/different_eta_gap/DifferentEtaGap_%iGeV_v%i_pt_%s%s_%i_%i.png",outPath, Energy,harmonic,PID.Data(), charge.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	}else{
		result->SaveAs(Form("%s/picture/3_PID/different_eta_gap/DifferentEtaGap_%iGeV_v%i_pt_%s%s_%i_%i.png",outPath, Energy,harmonic,PID.Data(), charge.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	}

	delete result;
	file_read->Close();
	file_out->Close();

}

void FlowDifferentMethodPID(Int_t Energy, Int_t harmonic, Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, TString PID, TString charge, Double_t max_y_axis, Double_t max_x_axis, Double_t MinYratio, Double_t MaxYratio){

	VEC	rebin_vec = Rebin2( harmonic, Energy, PID, charge, CentBinMin, CentBinMax);

	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<Draw_Object *> flow;
	std::vector<Draw_Object *> ratio;

	std::vector<Draw_Object *> flow_tpc;
	std::vector<Draw_Object *> ratio_tpc;
	
	TString suf = "(Def)";

	Double_t grXmin=0.012;
	//rebin_vec={0.2,0.4,0.6,0.8,1.0,1.2,2.0,2.2,5.0};
	//rebin_vec = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,5.0};
	rebin_vec = {0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.4,5.2};
	//rebin_vec = {0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.4,2.6,2.8,3.0,3.2};
	if(strncmp(PID, "Proton",6)==0){
		grXmin=0.352;
		//rebin_vec={0.4,0.6,0.8,1.0,1.2,1.4,1.6,5.0};
		rebin_vec = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,5.0};
	}


	const Int_t color[]={2, 1, 4, 6, 4, 2, 46};
	const Int_t style[]={23, 22, 25, 23, 34, 47, 8, 29};
	const Double_t size[]={2.,2.,2.,2.0,2.,2.,2.};

	TString method[]={"TPCandTOF","CombPID1","TPC","TPCandTOF","TPC"};
	//TString method_text[]={"TPC and TOF","Comb PID","TPC n#sigma=1.5","TPC n#sigma=2.0","TPC n#sigma=2.5"};
	TString method_text[]={"TPC and TOF","Comb PID","Comb PID (p=90)","TPC and TOF (cut 2)","TPC n#sigma=1.5","TPC n#sigma=2.0","TPC n#sigma=2.5"};
	TString nSigma_TPC[]={"15","20","25"};

	TFile *file_read;
	TFile *file_read2 = new TFile(Form("%s/file_pid/Flow_27GeV_comb_newCut.root",outPath, Energy),"READ");
	TFile *file_read_nSigma[3];


	if(Energy==39){
		file_read = new TFile(Form("%s/file_pid/nSigma/Flow_%iGeVper1_PID_10binYesTofdca11_sys_nSigma20.root",outPath, Energy),"READ");
		for(Int_t i=0; i<3; i++){
			file_read_nSigma[i] = new TFile(Form("%s/file_pid/nSigma/Flow_%iGeVper1_PID_10binYesTofdca11_sys_nSigma%s.root",outPath, Energy, nSigma_TPC[i].Data()),"READ");
		}
	}else{
		//file_read = new TFile(Form("%s/file_pid/Flow_%iGeV_CombPIDall_new.root",outPath, Energy),"READ");
		//file_read = new TFile(Form("%s/file_pid/Flow_27GeV_comb_good.root",outPath, Energy),"READ");
		//file_read = new TFile(Form("%s/file_pid/test27GeV10.root",outPath, Energy),"READ");
		//file_read = new TFile(Form("%s/file_pid/Flow_27GeV_comb_m55.root",outPath, Energy),"READ");
		//file_read = new TFile(Form("%s/macroPID/fitGausRun18/Flow_27GeV_Run18_comb.root",outPath, Energy),"READ");
		file_read = new TFile(Form("%s/macroPID/fitGausRun18/Flow_27GeV18_11_all_new.root",outPath, Energy),"READ");
		//file_read = new TFile(Form("%s/file_pid/test27GeV10.root",outPath, Energy),"READ");
		//file_read = new TFile(Form("%s/file_pid/Flow_27GeV_comb_newCut.root",outPath, Energy),"READ");
		//file_read = new TFile(Form("%s/file_pid/Flow_%iGeV_PID_10binYesTofdca11_sys.root",outPath, Energy),"READ");
	}

	TFile *file_out = new TFile("OUT.root","RECREATE");

	for(Int_t eta=0; eta<2; eta++){
	    flow.push_back(new Draw_Object);

	    if(eta==4){
			if(strncmp(PID, "Pion",4)==0){
				rebin_vec={0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,5.0};
			}
			if(strncmp(PID, "Kaon",4)==0){
				rebin_vec={0.2,0.25,0.3,0.35,0.4,0.45,0.5,5.0};
			}
			if(strncmp(PID, "Proton",6)==0){
				rebin_vec={0.5,0.55,0.6,0.7,0.8,0.9,1.0,1.1,5.0};
			}
	    }

	    if(eta<4){
    		Analysis -> SetParametrs(Energy, file_read,file_out,harmonic,PID,charge,EtaGap,method[eta]);
    		if(eta==2)Analysis -> SetParametrs(Energy, file_read,file_out,harmonic,PID,charge,EtaGap,"TPC");
    		if(eta==3)Analysis -> SetParametrs(Energy, file_read2,file_out,harmonic,PID,charge,EtaGap,"TPCandTOF");
    	}else{
    		Analysis -> SetParametrs(Energy, file_read_nSigma[eta-2],file_out,harmonic,PID,charge,EtaGap,method[eta]);
    	}
    	flow.back() -> SetParametrsGraph(Form("data%s",method[eta].Data()),Form("#font[42]{%s %s}", method_text[eta].Data(), suf.Data()), color[eta], style[eta], size[eta]);
    	Analysis -> FlowVsPt_EtaSub(flow.back(), CentBinMin, CentBinMax, rebin_vec);
    	suf="";

	}
    //Data line

    //ratio line
    for(Int_t i=1; i<flow.size();i++){
    	ratio.push_back(new Draw_Object);
    	Analysis -> RatioGraphEVAL(ratio[i-1],flow[i],flow[0],"");
    	std::cout<<"good\n";
	}

		MinYratio=0.945;
		MaxYratio=1.055;

    //canvas setening
    Draw_Picture_new *canvas = new Draw_Picture_new;
    max_x_axis = 3.13;
   	max_y_axis = max_y_axis*1.5;
    canvas->SetTLine(grXmin, 0., max_x_axis, 0., 1, 2, 2, "Up");
    canvas->SetTLine(grXmin, 0., max_x_axis, 0., 1, 2, 2, "Up");
    canvas->SetTLine(grXmin, 0.9, max_x_axis, 0.9, 1, 2, 2, "Down");
    canvas->SetTLine(grXmin, 1.1, max_x_axis, 1.1, 1, 2, 2, "Down");
    canvas->SetTLine(grXmin, 1., max_x_axis, 1., flow[0] -> GetColor(), 2, 1, "Down");
	canvas->SetAxisToCanvsWithBottomPanel(grXmin, max_x_axis, -0.01, max_y_axis, MinYratio, MaxYratio);
	canvas->SetDrawObject(flow,"Up");
	canvas->SetDrawObject(ratio,"Down");
	canvas->SetLegend(0.1,0.5,0.4,0.88,1);
	TCanvas *result = canvas->CanvasWithBottomPanel(Form("#font[42]{ #scale[1.2]{ #splitline{Au+Au #sqrt{s_{NN}} = %s GeV}{%s, centrality %i-%i %%}}}",Energy_scan_text.at(Energy), particle[GetNumberParticle(PID,charge)].Data(),cent[CentBinMax+1], cent[CentBinMin]),(0.8*grXmin) +0.36*max_x_axis,0.8*max_y_axis,
													"#font[42]{p_{T} [GeV/c]}",Form("v_{%i}",harmonic),"#font[42]{#frac{Data}{Def}}","");
	
	result->SaveAs(Form("/home/demanov/New_work/MyAnalysis/macroPID/fitGausRun18/pict2/1/testDifferenPIDcomb_%iGeVRun18_v%i_pt_%s%s_%i_%i_new.png", Energy,harmonic,PID.Data(), charge.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	//result->SaveAs(Form("%s/picture/3_PID/different_PID/testDifferenPIDcomb_%iGeVRun18_v%i_pt_%s%s_%i_%i_new.png",outPath, Energy,harmonic,PID.Data(), charge.Data(), cent[CentBinMax+1], cent[CentBinMin]));

	delete result;
	/*

	if(strncmp(PID, "Pion",4)==0){
		rebin_vec={0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,5.0};
	}
	if(strncmp(PID, "Kaon",4)==0){
		rebin_vec={0.2,0.25,0.3,0.35,0.4,0.45,0.5,5.0};
	}
	if(strncmp(PID, "Proton",6)==0){
		rebin_vec={0.5,0.55,0.6,0.7,0.8,0.9,1.0,1.1,5.0};
	}

	for(Int_t eta=0; eta<3; eta++){
		if(eta==1){
    		suf="Def";
    	}
	    flow_tpc.push_back(new Draw_Object);
    	Analysis -> SetParametrs(Energy, file_read_nSigma[eta],file_out,harmonic,PID,charge,EtaGap,"TPC");
    	flow_tpc.back() -> SetParametrsGraph(Form("dataTPC%i",eta),Form("#font[42]{%s %s}", method_text[eta+2].Data(), suf.Data()), color[eta], style[eta], size[eta]);
    	Analysis -> FlowVsPt_EtaSub(flow_tpc.back(), CentBinMin, CentBinMax, rebin_vec);
    	suf="";

	}
    //Data line

    //ratio line
    ratio_tpc.push_back(new Draw_Object);
    Analysis -> RatioGraphPointToPoint(ratio_tpc.back(),flow_tpc[0],flow_tpc[1],"");
    ratio_tpc.push_back(new Draw_Object);
    Analysis -> RatioGraphPointToPoint(ratio_tpc.back(),flow_tpc[2],flow_tpc[1],"");

    //canvas setening
    Draw_Picture_new *canvas2 = new Draw_Picture_new;
    
    canvas2->SetTLine(grXmin, 0., max_x_axis, 0., 1, 2, 2, "Up");
    canvas2->SetTLine(grXmin, 0., max_x_axis, 0., 1, 2, 2, "Up");
    canvas2->SetTLine(grXmin, 0.9, max_x_axis, 0.9, 1, 2, 2, "Down");
    canvas2->SetTLine(grXmin, 1.1, max_x_axis, 1.1, 1, 2, 2, "Down");
    canvas2->SetTLine(grXmin, 1., max_x_axis, 1., flow[1] -> GetColor(), 2, 1, "Down");
	canvas2->SetAxisToCanvsWithBottomPanel(grXmin, max_x_axis, -0.01, max_y_axis, 0.85, 1.15);
	canvas2->SetDrawObject(flow_tpc,"Up");
	canvas2->SetDrawObject(ratio_tpc,"Down");
	canvas2->SetLegend(0.1,0.55,0.4,0.88,1);
	TCanvas *result2 = canvas2->CanvasWithBottomPanel(Form("#font[42]{ #scale[1.2]{ #splitline{Au+Au #sqrt{s_{NN}} = %s GeV}{%s, centrality %i-%i %%}}}",Energy_scan_text.at(Energy), particle[GetNumberParticle(PID,charge)].Data(),cent[CentBinMax+1], cent[CentBinMin]),grXmin + 0.36*max_x_axis,0.8*max_y_axis,
													"#font[42]{p_{T} [GeV/c]}",Form("v_{%i}",harmonic),"#font[42]{#frac{Data}{Def}}","");
	
	result2->SaveAs(Form("%s/picture/3_PID/different_PID/DifferenPIDonlyTPC_%iGeV_v%i_pt_%s%s_%i_%i.png",outPath, Energy,harmonic,PID.Data(), charge.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result2;
	*/
	file_read->Close();
	file_out->Close();

}

void ParAntFlowVsPtForBES( Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, TString PID, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio, Double_t min_y_axis_dif, Double_t max_y_axis_dif){

	VEC rebin_vec;
	Int_t point = 1; 
	TString prefix = "TPCandTOF";
	
	Int_t Arr_energy[6]={11,14,19,27,39,62};
	const Int_t color[]={2, 4, 1};
	const Int_t style_par[]={33, 21, 22};
	const Int_t style_ant[]={27, 25, 26};
	Float_t size[3]={2.2,2.2,2.2};

	TFile *file_1[6];
	for(Int_t i=0; i<6; i++){
		file_1[i] = new TFile(Form("%s/file_pid/Flow_%iGeV_PID_10binYesTofdca11_sys.root",outPath, Arr_energy[i]),"READ");
	}
	TFile *file_out = new TFile("OUT.root","RECREATE");

	Draw_Picture_new *canvasFlow = new Draw_Picture_new;
	Draw_Picture_new *canvasDifference = new Draw_Picture_new;
	Draw_Picture_new *canvasRatio = new Draw_Picture_new;
	
	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<std::vector<Draw_Object *>> flow;
	std::vector<std::vector<Draw_Object *>> ratio;
	std::vector<std::vector<Draw_Object *>> difference;
	std::vector<Draw_Object *> object;
	
	//Double_t coord_y_draw_ratio = min_y_axis_ratio + (max_y_axis_ratio - min_y_axis_ratio)*0.8;
	Double_t coord_y_draw_ratio = min_y_axis_ratio + (max_y_axis_ratio - min_y_axis_ratio)*0.87;
	Double_t coord_y_draw_dif = min_y_axis_dif + (max_y_axis_dif - min_y_axis_dif)*0.8;
	
	Int_t b=0;

	for(Int_t j=4; j<6; j++){
		
		canvasFlow->SetTextInfo(Form("#font[42]{%s GeV}", Energy_scan_text.at(Arr_energy[j])), "Pad", 0.25, 0.75*max_y_axis, 0.13, 0.0);
		canvasDifference->SetTextInfo(Form("#font[42]{%s GeV}", Energy_scan_text.at(Arr_energy[j])), "Pad", 0.25, coord_y_draw_dif, 0.07, 0.0);
		canvasRatio->SetTextInfo(Form("#font[42]{%s GeV}", Energy_scan_text.at(Arr_energy[j])), "Pad", 0.25, coord_y_draw_ratio, 0.07, 0.0);
		
		for(Int_t harm=0; harm<1; harm++){

			rebin_vec = Rebin2K( harm+2, Arr_energy[j], PID, "", CentBinMin, CentBinMax);
			
			object.push_back(new Draw_Object);
			Analysis -> SetParametrs(Arr_energy[j], file_1[j], file_out, harm+2, PID, "Pos", EtaGap, prefix);
			Analysis -> FlowVsPt_EtaSub(object.back(), CentBinMin, CentBinMax, rebin_vec);
    		object.back() -> SetParametrsGraph(Form("data_v%i",harm+2),Form("#font[42]{#scale[1.2]{%s}}", particle[GetNumberParticle(PID,"Pos")].Data()), color[0], style_par[0], size[0]);
			
			object.push_back(new Draw_Object);
			Analysis -> SetParametrs(Arr_energy[j], file_1[j], file_out, harm+2, PID, "Neg", EtaGap, prefix);
			Analysis -> FlowVsPt_EtaSub(object.back(), CentBinMin, CentBinMax, rebin_vec);
    		object.back() -> SetParametrsGraph(Form("data_v%i",harm+2),Form("#font[42]{#scale[1.2]{%s}}", particle[GetNumberParticle(PID,"Neg")].Data()), color[1], style_ant[0], size[0]);
			
			flow.push_back(object);
			object.clear();

			object.push_back(new Draw_Object);
			Analysis -> DifferentParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
    		object.back() -> SetParametrsGraph(Form("data_v%i",harm+2),Form("#font[42]{#scale[1.2]{%s}}", particle[GetNumberParticle(PID,"Pos")].Data()), color[0], style_par[0], size[0]);
    		difference.push_back(object);
    		object.clear();

			object.push_back(new Draw_Object);
			Analysis -> RatioParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
    		object.back() -> SetParametrsGraph(Form("data_v%i",harm+2),Form("#font[42]{#scale[1.2]{%s}}", particle[GetNumberParticle(PID,"Pos")].Data()), color[0], style_par[0], size[0]);
    		ratio.push_back(object);
    		object.clear();

			if(j==5){
				canvasFlow->SetTextInfo(Form("#font[42]{n=%i}", harm+2), Form("V%i",harm+2), 0.85*max_x_axis, 0.75*max_y_axis, 0.13, 0.0);
				canvasDifference->SetTextInfo(Form("#font[42]{n=%i}", harm+2), Form("V%i",harm+2), 0.83*max_x_axis,coord_y_draw_dif, 0.07, 0.0);
				canvasRatio->SetTextInfo(Form("#font[42]{n=%i}", harm+2), Form("V%i",harm+2), 0.85*max_x_axis, coord_y_draw_ratio, 0.07, 0.0);
			}

			canvasFlow->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.75*max_y_axis, 0.13, 0.0);
			canvasDifference->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, coord_y_draw_dif, 0.07, 0.0);
			canvasRatio->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, coord_y_draw_ratio, 0.07, 0.0);
			b++;

		}
	}
	canvasFlow->SetTextInfo("#font[42]{v_{n}}","Y",0.25,0.6,0.5,0.0);
	canvasFlow->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %i-%i %% }}", cent[CentBinMax+1], cent[CentBinMin]),"Info",1.1,0.75*max_y_axis,0.13,0.0);
	//canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{STAR Preliminary }}"),"Prel",1.1,0.75*max_y_axis,0.13,0.0);
	canvasFlow->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasFlow->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.13);
	canvasFlow->SetDrawObject(flow[0],"OnePad");
	canvasFlow->SetLegend(0.15,0.5,0.35,0.92,1);
	TCanvas *result = canvasFlow->CanvasNxM(1000,1280, 2, 1, 0.4, 0.2, flow,1,0,1);

	result->SaveAs(Form("%s/picture/3_PID/BES_vn_%s_%s_cent_%i_%i.pdf",outPath, PID.Data(), EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result->SaveAs(Form("%s/picture/3_PID/BES_vn_%s_%s_cent_%i_%i.png",outPath, PID.Data(), EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result;

	canvasRatio->SetTextInfo(Form("#font[42]{#frac{ v_{n}(%s) }{ v_{n}(%s) } }", particle[GetNumberParticle(PID,"Pos")].Data(), particle[GetNumberParticle(PID,"Neg")].Data()),"Y",0.07,0.7,0.37,0);
	canvasRatio->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasRatio->SetTextInfo(Form("#font[42]{ #scale[1.0]{#splitline{ Au+Au, %i-%i %% }{This analysis}}}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.7,coord_y_draw_ratio,0.07,0.0);
	//canvasRatio->SetTextInfo(Form("#font[42]{ #scale[1.0]{ STAR Preliminary}}"),"Prel",1.1,coord_y_draw_ratio,0.13,0.0);
	canvasRatio->SetTLine(min_x_axis, 1.05, max_x_axis, 1.05, 1, 2, 2, "Up");
	canvasRatio->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasRatio->SetTLine(min_x_axis, 0.95, max_x_axis, 0.95, 1, 2, 2, "Up");
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_ratio, max_y_axis_ratio, 0.07);
	canvasRatio->SetDrawObject(ratio[0],"OnePad");
	//canvasRatio->SetLegend(0.1,0.4,0.3,0.92,1,0,1);
	TCanvas *result2 = canvasRatio->CanvasNxM(1000,1280, 2, 1, 0.3, 0.25, ratio,1,0,1);

	result2->SaveAs(Form("%s/picture/3_PID/BES_vn_ratio_%s_%s_cent_%i_%i.pdf",outPath, PID.Data(), EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result2->SaveAs(Form("%s/picture/3_PID/BES_vn_ratio_%s_%s_cent_%i_%i.png",outPath, PID.Data(), EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result2;

	canvasDifference->SetTextInfo(Form("#font[42]{ v_{n}(%s) - v_{n}(%s) }",particle[GetNumberParticle(PID,"Pos")].Data(),particle[GetNumberParticle(PID,"Neg")].Data() ),"Y",0.4, 0.5,0.4,90);
	canvasDifference->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.37,0.0);
	canvasDifference->SetTextInfo(Form("#font[42]{ #scale[1.0]{#splitline{ Au+Au, %i-%i %% }{This analysis}}}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.7,coord_y_draw_dif,0.07,0.0);
	//canvasDifference->SetTextInfo(Form("#font[42]{ #scale[1.0]{ STAR Preliminary }}"),"Prel",1.1,coord_y_draw_dif,0.13,0.0);
	canvasDifference->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasDifference->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_dif, max_y_axis_dif, 0.07);
	canvasDifference->SetDrawObject(difference[0],"OnePad");
	//canvasDifference->SetLegend(0.1,0.7,0.3,0.92,1,0,1);
	TCanvas *result3 = canvasDifference->CanvasNxM(1150,1280, 2, 1, 0.3, 0.25, difference,1,0,1);

	result3->SaveAs(Form("%s/picture/3_PID/BES_vn_diff_%s_%s_cent_%i_%i.pdf",outPath, PID.Data(), EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result3->SaveAs(Form("%s/picture/3_PID/BES_vn_diff_%s_%s_cent_%i_%i.png",outPath, PID.Data(), EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result3;

	for(Int_t i=0; i<6; i++){
		file_1[i]->Close();
		delete file_1[i];
	}

}

void FlowVsPtForBESmultyPad( Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, TString charge, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio, Double_t min_y_axis_dif, Double_t max_y_axis_dif){

	VEC rebin_vec;
	Int_t point = 1; 
	TString prefix = "TPCandTOF";
	
	Int_t Arr_energy[6]={11,14,19,27,39,62};
	const Int_t color[]={2, 1, 4, 6, 4, 2, 46};
	const Int_t style[]={24, 22, 23, 8, 25, 34, 8, 29};
	const Double_t size[]={1.5,1.5,1.5,1.5,1.5,1.5,1.5};
	TString arr_par[]={"Pion","Kaon","Proton"};
	Float_t k[2]={1.0,2.0};

	TFile *file_1[6];
	for(Int_t i=0; i<6; i++){
		file_1[i] = new TFile(Form("%s/file_pid/Flow_%iGeV_PID_10binYesTofdca11_sys.root",outPath, Arr_energy[i]),"READ");
	}
	TFile *file_out = new TFile("OUT.root","RECREATE");

	Draw_Picture_new *canvasFlow = new Draw_Picture_new;
	
	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<std::vector<Draw_Object *>> flow;
	std::vector<Draw_Object *> object;
	
	Int_t b=0;

	for(Int_t par=0; par<3; par++){
		for(Int_t harm=0; harm<2; harm++){
			for(Int_t j=0; j<6; j++){

				rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], charge, CentBinMin, CentBinMax);
				
				object.push_back(new Draw_Object);
				Analysis -> SetParametrs(Arr_energy[j], file_1[j], file_out, harm+2, arr_par[par], charge, EtaGap, prefix);
				Analysis -> FlowVsPt_EtaSub(object.back(), CentBinMin, CentBinMax, rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i",harm+2),Form("#font[42]{#scale[1.3]{%s GeV}}", Energy_scan_text.at(Arr_energy[j])), color[j], style[j], size[j]);

				if(harm==1){
					for(Int_t n=0; n < object.back()->GetSizeVector(); n++){
						object.back()->СhangePointGraph(n,
														object.back()->GetPointXGraph(n), 
														object.back()->GetPointYGraph(n) * k[1],
														object.back()->GetPointXErrorGraph(n),
														object.back()->GetPointYErrorGraph(n) * k[1],
														object.back()->GetPointYSysErrorGraph(n) * k[1]);
					}
				}

				if(j==5 && par==2){
					if(harm==0) canvasFlow->SetTextInfo(Form("#font[42]{(v_{%i})}",harm+2), Form("V%i",harm+2), 0.6*max_x_axis, 0.85*max_y_axis, 0.07, 0.0);
					if(harm==1) canvasFlow->SetTextInfo(Form("#font[42]{(%.1f #times v_{%i})}", k[harm],harm+2), Form("V%i",harm+2), 0.6*max_x_axis, 0.85*max_y_axis, 0.07, 0.0);
				}
			}
			canvasFlow->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.85*max_y_axis, 0.07, 0.0);
			b++;

			flow.push_back(object);
			object.clear();
		}
	}

	canvasFlow->SetTextInfo(Form("#font[42]{v_{n}(%s)}",particle[GetNumberParticle("Pion",charge)].Data()),"Y",0.45,0.75,0.43,90);
	canvasFlow->SetTextInfo(Form("#font[42]{v_{n}(%s)}",particle[GetNumberParticle("Kaon",charge)].Data()),"Y",0.45,0.4,0.43,90);
	canvasFlow->SetTextInfo(Form("#font[42]{v_{n}(%s)}",particle[GetNumberParticle("Proton",charge)].Data()),"Y",0.45,0.1,0.43,90);
	canvasFlow->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.45,0.6,0.43,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %i-%i %% }}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.2,0.85*max_y_axis,0.07,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{ STAR Preliminary }}"),"Prel",0.2,0.85*max_y_axis,0.07,0.0);
	canvasFlow->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasFlow->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.07);
	canvasFlow->SetDrawObject(flow[0],"OnePad");
	canvasFlow->SetLegend(0.1,0.75,0.98,0.98,3);
	TCanvas *result = canvasFlow->CanvasNxM(1200,1280, 3, 2, 0.4, 0.35, flow,1,0,2);

	result->SaveAs(Form("%s/picture/3_PID/multPad_BES_vn_%s_%s_cent_%i_%i.pdf",outPath, charge.Data(),EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result->SaveAs(Form("%s/picture/3_PID/multPad_BES_vn_%s_%s_cent_%i_%i.png",outPath, charge.Data(),EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result->SaveAs(Form("%s/picture/3_PID/multPad_BES_vn_%s_%s_cent_%i_%i.C",outPath, charge.Data(),EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result;

	for(Int_t i=0; i<6; i++){
		file_1[i]->Close();
		delete file_1[i];
	}

}

void FlowVsPtForBESmultyPad_Int( Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, TString charge, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio, Double_t min_y_axis_dif, Double_t max_y_axis_dif){

	VEC rebin_vec; 
	TString prefix = "TPCandTOF";

	VEC rebin_vec_integral;

	Int_t Arr_energy[6]={11,14,19,27,39,62};
	const Int_t color[]={2, 1, 4, 6, 4, 2, 46};
	const Int_t style[]={24, 22, 23, 8, 25, 34, 8, 29};
	const Double_t size[]={1.5,1.5,1.5,1.5,1.5,1.5,1.5};
	TString arr_par[]={"Pion","Kaon","Proton"};
	Float_t k[2]={1.0,1.0};

	TFile *file_1[6];
	for(Int_t i=0; i<6; i++){
		file_1[i] = new TFile(Form("%s/file_pid/Flow_%iGeV_PID_10binYesTofdca11_sys.root",outPath, Arr_energy[i]),"READ");
	}
	TFile *file_out = new TFile("OUT.root","RECREATE");

	Draw_Picture_new *canvasFlow = new Draw_Picture_new;
	
	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<std::vector<Draw_Object *>> flow;
	std::vector<Draw_Object *> object;
	
	Int_t b=0;

	for(Int_t par=0; par<3; par++){

		rebin_vec_integral={0.2,2.0,5.2};
		if(strncmp(arr_par[par],"Proton",6)==0){
			rebin_vec_integral={0.5,2.0,5.2};	
		}
		
		for(Int_t harm=0; harm<2; harm++){
			for(Int_t j=0; j<6; j++){

				rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], charge, CentBinMin, CentBinMax);
				
				object.push_back(new Draw_Object);
				Analysis -> SetParametrs(Arr_energy[j], file_1[j], file_out, harm+2, arr_par[par], charge, EtaGap, prefix);
	    		Analysis -> FlowVsPtDivideIntegralFLow_EtaSub_PID(object.back(), CentBinMin, CentBinMax, rebin_vec, rebin_vec_integral, 0);
				//Analysis -> FlowVsPt_EtaSub(object.back(), CentBinMin, CentBinMax, rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i_%iGeV_%s",harm+2, Arr_energy[j],arr_par[par].Data()),Form("#font[42]{#scale[1.3]{%s GeV}}", Energy_scan_text.at(Arr_energy[j])), color[j], style[j], size[j]);

				if(harm==1){
					for(Int_t n=0; n < object.back()->GetSizeVector(); n++){
						object.back()->СhangePointGraph(n,
														object.back()->GetPointXGraph(n), 
														object.back()->GetPointYGraph(n) * k[1],
														object.back()->GetPointXErrorGraph(n),
														object.back()->GetPointYErrorGraph(n) * k[1],
														object.back()->GetPointYSysErrorGraph(n) * k[1]);
					}
				}

				if(j==5 && par==2){
					if(harm==0) canvasFlow->SetTextInfo(Form("#font[42]{(v_{%i})}",harm+2), Form("V%i",harm+2), 0.85*max_x_axis, 0.85*max_y_axis, 0.07, 0.0);
					if(harm==1) canvasFlow->SetTextInfo(Form("#font[42]{(v_{%i})}",harm+2), Form("V%i",harm+2), 0.85*max_x_axis, 0.85*max_y_axis, 0.07, 0.0);
				}
			}
			canvasFlow->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.85*max_y_axis, 0.07, 0.0);
			b++;

			flow.push_back(object);
			object.clear();
		}
	}

	canvasFlow->SetTextInfo(Form("#font[42]{#frac{v_{n}(%s)}{v_{n}^{int}}}",particle[GetNumberParticle("Pion",charge)].Data()),"Y",0.25,0.75,0.33,0);
	canvasFlow->SetTextInfo(Form("#font[42]{#frac{v_{n}(%s)}{v_{n}^{int}}}",particle[GetNumberParticle("Kaon",charge)].Data()),"Y",0.25,0.4,0.33,0);
	canvasFlow->SetTextInfo(Form("#font[42]{#frac{v_{n}(%s)}{v_{n}^{int}}}",particle[GetNumberParticle("Proton",charge)].Data()),"Y",0.25,0.1,0.33,0);
	canvasFlow->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.45,0.6,0.43,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %i-%i %% }}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.2,0.85*max_y_axis,0.07,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{ STAR Preliminary }}"),"Prel",0.2,0.85*max_y_axis,0.07,0.0);
	canvasFlow->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasFlow->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.07);
	canvasFlow->SetDrawObject(flow[0],"OnePad");
	canvasFlow->SetLegend(0.1,0.75,0.98,0.98,3);
	TCanvas *result = canvasFlow->CanvasNxM(1200,1280, 3, 2, 0.4, 0.35, flow,1,0,2);

	result->SaveAs(Form("%s/picture/3_PID/Sys/int/Int_multPad_BES_vn_%s_%s_cent_%i_%i.pdf",outPath, charge.Data(),EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result->SaveAs(Form("%s/picture/3_PID/Sys/int/Int_multPad_BES_vn_%s_%s_cent_%i_%i.png",outPath, charge.Data(),EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result->SaveAs(Form("%s/picture/3_PID/Sys/int/Int_multPad_BES_vn_%s_%s_cent_%i_%i.C",outPath, charge.Data(),EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result;

	for(Int_t i=0; i<6; i++){
		file_1[i]->Close();
		delete file_1[i];
	}

}
/*
void ParAntFlowVsPtForBESupdate( Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, Double_t min_x_axis, Double_t max_x_axis){
	/// это общая версия, нужная
	VEC rebin_vec;
	Int_t point = 1; 
	TString prefix = "TPCandTOF";

	Int_t Arr_energy[7]={27,11,14,19,27,39,62};
	const Int_t color[]={2, 1, 4, 6, 4, 2, 46};
	const Int_t style[]={22, 24, 23, 8, 25, 34, 8, 29};
	const Double_t size[]={1.3,1.4,1.3,1.4,1.3,1.3,1.3};
	TString arr_par[]={"Pion","Kaon","Proton"};

	TFile *file_1[7];
	for(Int_t i=0; i<6; i++){
		file_1[i] = new TFile(Form("%s/file_pid/Flow_%iGeV_PID_10binYesTofdca11_sys.root",outPath, Arr_energy[i]),"READ");
	}
	TFile *file_out = new TFile("OUT.root","RECREATE");

	Draw_Picture_new *canvasDif = new Draw_Picture_new;
	Draw_Picture_new *canvasRatio = new Draw_Picture_new;
	
	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<std::vector<Draw_Object *>> flow_ratio;
	std::vector<std::vector<Draw_Object *>> flow_difference;
	std::vector<Draw_Object *> object;
	
	//Double_t min_y_axis_ratio[3]={0.78,0.68,0.8};
	Double_t min_y_axis_ratio[3]={0.93,0.86,0.8};
	//Double_t max_y_axis_ratio[3]={1.06,1.18,2.23};
	Double_t max_y_axis_ratio[3]={1.03,1.18,2.23};
	
	Double_t min_y_axis_dif[3]={-0.0042, -0.0032, -0.007};
	//Double_t min_y_axis_dif[3]={-0.006, -0.011, -0.007};
	//Double_t min_y_axis_dif[3]={-0.022, -0.022, -0.012};
	//Double_t max_y_axis_dif[3]={0.006, 0.006, 0.023};
	Double_t max_y_axis_dif[3]={0.0036, 0.0042, 0.023};
	//Double_t max_y_axis_dif[3]={0.017, 0.022, 0.044};

	Double_t coord_y_draw_dif[3];
	Double_t coord_y_draw_ratio[3];

	Int_t b=0;
	Int_t b2=0;

	for(Int_t par=0; par<3; par++){

		coord_y_draw_dif[par]=min_y_axis_dif[par] + (max_y_axis_dif[par] - min_y_axis_dif[par])*0.9;
		coord_y_draw_ratio[par]=min_y_axis_ratio[par] + (max_y_axis_ratio[par] - min_y_axis_ratio[par])*0.9;
		if(par==0){
			//coord_y_draw_ratio[par]=min_y_axis_ratio[par] + (max_y_axis_ratio[par] - min_y_axis_ratio[par])*0.99;
		}

		for(Int_t harm=0; harm<2; harm++){
		
			for(Int_t j=0; j<1; j++){

				rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], "", CentBinMin, CentBinMax);
	    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"",EtaGap,prefix);
				object.push_back(new Draw_Object);
				Analysis -> DifferentParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i_%iGeV_%s",harm+2,Arr_energy[j],arr_par[par].Data()),Form("#font[42]{#scale[1.0]{%s GeV Run18}}", Energy_scan_text.at(Arr_energy[j])), color[par], style[j], size[j]);
			}
			flow_difference.push_back(object);
			object.clear();
			
			if(par==2){
				canvasDif->SetTextInfo(Form("#font[42]{v_{%i}}", harm+2), Form("V%i",harm+2), 0.85*max_x_axis,0.95*coord_y_draw_dif[0], 0.08, 0.0);
			}
			canvasDif->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.95*coord_y_draw_dif[par], 0.08, 0.0);
			b++;			

		}

		for(Int_t harm=0; harm<2; harm++){
		
			for(Int_t j=0; j<1; j++){

				rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], "", CentBinMin, CentBinMax);
	    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"",EtaGap,prefix);
				object.push_back(new Draw_Object);
				Analysis -> RatioParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i_%iGeV_%s",harm+2,Arr_energy[j],arr_par[par].Data()),Form("#font[42]{#scale[1.0]{%s GeV Run18}}", Energy_scan_text.at(Arr_energy[j])), color[par], style[j], size[j]);
			}
			flow_ratio.push_back(object);
			object.clear();
			if(par==2){
				canvasRatio->SetTextInfo(Form("#font[42]{v_{%i}}", harm+2), Form("V%i",harm+2), 0.85*max_x_axis,0.99*coord_y_draw_ratio[0], 0.08, 0.0);
			}
			canvasRatio->SetNumberPad(Form("#font[42]{%s}", bukva[b2].Data()), "Pad", 0.04, 0.99*coord_y_draw_ratio[par], 0.08, 0.0);
			b2++;
		}
	}

	canvasDif->SetTextInfo(Form("#font[42]{ v_{n}(%s) - v_{n}(%s) }",particle[GetNumberParticle("Pion","Pos")].Data(),particle[GetNumberParticle("Pion","Neg")].Data() ),"Y",0.4, 0.7,0.4,90);
	canvasDif->SetTextInfo(Form("#font[42]{ v_{n}(%s) - v_{n}(%s) }",particle[GetNumberParticle("Kaon","Pos")].Data(),particle[GetNumberParticle("Kaon","Neg")].Data() ),"Y",0.4, 0.4,0.4,90);
	canvasDif->SetTextInfo(Form("#font[42]{ v_{n}(%s) - v_{n}(%s) }",particle[GetNumberParticle("Proton","Pos")].Data(),particle[GetNumberParticle("Proton","Neg")].Data() ),"Y",0.4, 0.07,0.4,90);
	canvasDif->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasDif->SetTextInfo(Form("#font[42]{ #scale[1.0]{#splitline{Au+Au #sqrt{s_{NN}}=27 GeV}{Centrality %i-%i%%}}}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.3,0.8*coord_y_draw_dif[0],0.08,0.0);//0.8
	//canvasDif->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %i-%i %% }}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.6,0.95*coord_y_draw_dif[0],0.08,0.0);
	canvasDif->SetTextInfo(Form("#font[42]{ #scale[1.0]{ STAR Preliminary }}"),"Prel",0.6,0.95*coord_y_draw_dif[2],0.08,0.0);
	canvasDif->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasDif->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasDif->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_dif[0], max_y_axis_dif[0], 0.06);
	canvasDif->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_dif[1], max_y_axis_dif[1], 0.06);
	canvasDif->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_dif[2], max_y_axis_dif[2], 0.06);
	canvasDif->SetDrawObject(flow_difference[0],"OnePad");
	canvasDif->SetLegend(0.07,0.08,0.98,0.2,3);
	TCanvas *result = canvasDif->CanvasNxM(1000,1280, 3, 2, 0.4, 0.35, flow_difference,3,0,1);

	result->SaveAs(Form("%s/picture/3_PID/Sys/test/only27BESupdateDif_vn_%s_cent_%i_%i_new.pdf",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result->SaveAs(Form("%s/picture/3_PID/Sys/test/only27BESupdateDif_vn_%s_cent_%i_%i_new.png",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result->SaveAs(Form("%s/picture/3_PID/Sys/test/only27BESupdateDif_vn_%s_cent_%i_%i_fin.C",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result;

	canvasRatio->SetTextInfo(Form("#font[42]{#frac{ v_{n}(%s) }{ v_{n}(%s) } }", particle[GetNumberParticle("Pion","Pos")].Data(), particle[GetNumberParticle("Pion","Neg")].Data()),"Y",0.07,0.8,0.3,0);
	canvasRatio->SetTextInfo(Form("#font[42]{#frac{ v_{n}(%s) }{ v_{n}(%s) } }", particle[GetNumberParticle("Kaon","Pos")].Data(), particle[GetNumberParticle("Kaon","Neg")].Data()),"Y",0.07,0.5,0.3,0);
	canvasRatio->SetTextInfo(Form("#font[42]{#frac{ v_{n}(%s) }{ v_{n}(%s) } }", particle[GetNumberParticle("Proton","Pos")].Data(), particle[GetNumberParticle("Proton","Neg")].Data()),"Y",0.07,0.15,0.3,0);
	canvasRatio->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasRatio->SetTextInfo(Form("#font[42]{ #scale[1.0]{#splitline{Au+Au #sqrt{s_{NN}}=27 GeV}{Centrality %i-%i%%}}}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.3,0.99*coord_y_draw_ratio[0],0.08,0.0);//0.98//0.8
	//canvasRatio->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %i-%i %% }}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.6,0.98*coord_y_draw_ratio[0],0.08,0.0);
	canvasRatio->SetTextInfo(Form("#font[42]{ #scale[1.0]{ STAR Preliminary }}"),"Prel",0.6,0.96*coord_y_draw_ratio[2],0.08,0.0);
	canvasRatio->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasRatio->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_ratio[0], max_y_axis_ratio[0], 0.06);
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_ratio[1], max_y_axis_ratio[1], 0.06);
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_ratio[2], max_y_axis_ratio[2], 0.06);
	canvasRatio->SetDrawObject(flow_ratio[0],"OnePad");
	canvasRatio->SetLegend(0.04,0.7,0.98,0.85,3);
	TCanvas *result2 = canvasRatio->CanvasNxM(1000,1280, 3, 2, 0.4, 0.35, flow_ratio,3,0,0);

	result2->SaveAs(Form("%s/picture/3_PID/Sys/test/only27BESupdateRatio_vn_%s_cent_%i_%i_new.pdf",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result2->SaveAs(Form("%s/picture/3_PID/Sys/test/only27BESupdateRatio_vn_%s_cent_%i_%i_new.png",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result2->SaveAs(Form("%s/picture/3_PID/Sys/test/only27BESupdateRatio_vn_%s_cent_%i_%i_fin.C",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	//result2->SaveAs(Form("%s/picture/3_PID/Sys/test/only27BESupdateRatio_vn_%s_cent_%i_%i_new.cpp",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result2;


	for(Int_t i=0; i<6; i++){
		file_1[i]->Close();
		delete file_1[i];
	}

}
*/

void ParAntFlowVsPtForBESupdate( Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, Double_t min_x_axis, Double_t max_x_axis){

	VEC rebin_vec;

	Int_t point = 1; 
	TString prefix = "TPCandTOF";
	//TString prefix2 = "TPCandTOF";
	TString prefix2 = "CombPID1";
	//TString prefix3 = "CombPID1";
	//TString prefix4 = "CombPID2";
	//TString prefix5 = "CombPID6";
	//TString prefix6 = "CombPID5";
	//TString prefix7 = "CombPID6";

	Int_t Arr_energy[7]={27,11,14,19,27,39,62};
	const Int_t color[]={2, 1, 4, 6, 4, 2, 46};
	const Int_t style[]={22, 28, 23, 8, 25, 34, 8, 29};

	const Int_t color2[]={4, 6, 4, 2, 46};
	const Int_t style2[]={32, 23, 25, 34, 8, 29};
	
	const Int_t colorall[]={2, 1, 4, 6, 2, 1, 4, 7};
	const Int_t styleall[]={34, 22, 23, 20, 26, 32, 28, 21};


	//const Double_t size[]={1.3,1.4,1.3,1.4,1.3,1.3,1.3};
	const Double_t size[]={2.0,2.1,2.0,2.1,2.0,2.0,2.0};
	TString arr_par[]={"Pion","Kaon","Proton"};

	TFile *file_1[7];
	for(Int_t i=0; i<6; i++){
		file_1[i] = new TFile(Form("%s/file_pid/Flow_%iGeV_PID_10binYesTofdca11_sys.root",outPath, Arr_energy[i]),"READ");
	}
	//file_1[0] = new TFile(Form("%s/file_pid/Flow_27GeV_comb_m55.root",outPath),"READ");
	//TFile *file_111 = new TFile(Form("%s/file_pid/Flow_27GeV_comb_newCut.root",outPath),"READ");
	//file_1[0] = new TFile(Form("%s/macroPID/fitGausRun18/Flow_27GeV_Run18_comb.root",outPath),"READ");
	file_1[1] = new TFile(Form("%s/macroPID/fitGausRun18/Flow_27GeV18_11_all.root",outPath),"READ");
	file_1[0] = new TFile(Form("%s/macroPID/fitGausRun18/Flow_27GeV18_11_all_new.root",outPath),"READ");
	//file_1[0] = new TFile(Form("%s/macroPID/fitGausRun18/StRuns15.root",outPath),"READ");
	//file_1[0] = new TFile(Form("%s/file_pid/test27GeV10.root",outPath),"READ");

	TFile *file_out = new TFile("OUT.root","RECREATE");

	Draw_Picture_new *canvasDif = new Draw_Picture_new;
	Draw_Picture_new *canvasRatio = new Draw_Picture_new;
	
	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<std::vector<Draw_Object *>> flow_ratio;
	std::vector<std::vector<Draw_Object *>> flow_difference;
	std::vector<Draw_Object *> object;
	
	//Double_t min_y_axis_ratio[3]={0.78,0.68,0.8};
	Double_t min_y_axis_ratio[3]={0.947,0.91,0.8};
	//Double_t max_y_axis_ratio[3]={1.06,1.18,2.23};
	Double_t max_y_axis_ratio[3]={1.03,1.09,2.23};
	
	Double_t min_y_axis_dif[3]={-0.0042, -0.0092, -0.004};
	//Double_t min_y_axis_dif[3]={-0.006, -0.011, -0.007};
	//Double_t min_y_axis_dif[3]={-0.022, -0.022, -0.012};
	//Double_t max_y_axis_dif[3]={0.006, 0.006, 0.023};
	Double_t max_y_axis_dif[3]={0.0026, 0.0062, 0.023};
	//Double_t max_y_axis_dif[3]={0.017, 0.022, 0.044};

	Double_t coord_y_draw_dif[3];
	Double_t coord_y_draw_ratio[3];

	Int_t b=0;
	Int_t b2=0;

	for(Int_t par=0; par<1; par++){

		coord_y_draw_dif[par]=min_y_axis_dif[par] + (max_y_axis_dif[par] - min_y_axis_dif[par])*0.9;
		coord_y_draw_ratio[par]=min_y_axis_ratio[par] + (max_y_axis_ratio[par] - min_y_axis_ratio[par])*0.9;
		if(par==0){
			//coord_y_draw_ratio[par]=min_y_axis_ratio[par] + (max_y_axis_ratio[par] - min_y_axis_ratio[par])*0.99;
		}

		for(Int_t harm=0; harm<2; harm++){
		
			for(Int_t j=0; j<1; j++){

				//rebin_vec={0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.4};
				
				rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], "", CentBinMin, CentBinMax);
	    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"",EtaGap,prefix);
				object.push_back(new Draw_Object);
				Analysis -> DifferentParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i_%iGeV_%s",harm+2,Arr_energy[j],arr_par[par].Data()),Form("#font[42]{#scale[1.0]{TPC and TOF}}", Energy_scan_text.at(Arr_energy[j])), color[par+1], style[j], size[j]);
			}
			for(Int_t j=0; j<1; j++){

				//rebin_vec={0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.4};
				
				rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], "", CentBinMin, CentBinMax);
	    		Analysis -> SetParametrs(Arr_energy[j], file_1[1],file_out,harm+2,arr_par[par],"",EtaGap,prefix2);
				object.push_back(new Draw_Object);
				Analysis -> DifferentParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i_%iGeV_%s",harm+2,Arr_energy[j],arr_par[par].Data()),Form("#font[42]{#scale[1.0]{Comb PID}}", Energy_scan_text.at(Arr_energy[j])), color[par], style[j+1], size[j]);
			}

			/////////////////////////////////////////////////////////////
			/*
			for(Int_t j=0; j<1; j++){

				//rebin_vec={0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.4};
				
				rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], "", CentBinMin, CentBinMax);
	    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"",EtaGap,prefix3);
				object.push_back(new Draw_Object);
				Analysis -> DifferentParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i_%iGeV_%s",harm+2,Arr_energy[j],arr_par[par].Data()),Form("#font[42]{#scale[1.0]{Comb PID 2}}", Energy_scan_text.at(Arr_energy[j])), colorall[2], styleall[2], size[j]);
			}
			
			for(Int_t j=0; j<1; j++){

				//rebin_vec={0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.4};
				
				rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], "", CentBinMin, CentBinMax);
	    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"",EtaGap,prefix4);
				object.push_back(new Draw_Object);
				Analysis -> DifferentParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i_%iGeV_%s",harm+2,Arr_energy[j],arr_par[par].Data()),Form("#font[42]{#scale[1.0]{Comb PID 3}}", Energy_scan_text.at(Arr_energy[j])), colorall[3], styleall[3], size[j]);
			}
			
			for(Int_t j=0; j<1; j++){

				//rebin_vec={0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.4};
				
				rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], "", CentBinMin, CentBinMax);
	    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"",EtaGap,prefix5);
				object.push_back(new Draw_Object);
				Analysis -> DifferentParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i_%iGeV_%s",harm+2,Arr_energy[j],arr_par[par].Data()),Form("#font[42]{#scale[1.0]{Comb PID 4}}", Energy_scan_text.at(Arr_energy[j])), colorall[4], styleall[4], size[j]);
			}
			
			for(Int_t j=0; j<1; j++){

				//rebin_vec={0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.4};
				
				rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], "", CentBinMin, CentBinMax);
	    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"",EtaGap,prefix6);
				object.push_back(new Draw_Object);
				Analysis -> DifferentParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i_%iGeV_%s",harm+2,Arr_energy[j],arr_par[par].Data()),Form("#font[42]{#scale[1.0]{Comb PID 5}}", Energy_scan_text.at(Arr_energy[j])), colorall[5], styleall[5], size[j]);
			}

			for(Int_t j=0; j<1; j++){

				//rebin_vec={0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.4};
				
				rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], "", CentBinMin, CentBinMax);
	    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"",EtaGap,prefix7);
				object.push_back(new Draw_Object);
				Analysis -> DifferentParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i_%iGeV_%s",harm+2,Arr_energy[j],arr_par[par].Data()),Form("#font[42]{#scale[1.0]{Comb PID 6}}", Energy_scan_text.at(Arr_energy[j])), colorall[6], styleall[6], size[j]);
			}
			*/
			////////////////////////////////////////////////////////////////////////////////////////
			
			flow_difference.push_back(object);
			object.clear();


			if(par==2){
				canvasDif->SetTextInfo(Form("#font[42]{v_{%i}}", harm+2), Form("V%i",harm+2), 0.85*max_x_axis,0.95*coord_y_draw_dif[2], 0.08, 0.0);
			}
			canvasDif->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.95*coord_y_draw_dif[par], 0.08, 0.0);
			b++;			

		}

		for(Int_t harm=0; harm<2; harm++){
		
			for(Int_t j=0; j<1; j++){

				//rebin_vec={0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.4};
				
				rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], "", CentBinMin, CentBinMax);
	    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"",EtaGap,prefix);
				object.push_back(new Draw_Object);
				Analysis -> RatioParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i_%iGeV_%s",harm+2,Arr_energy[j],arr_par[par].Data()),Form("#font[42]{#scale[1.0]{TPC and TOF}}", Energy_scan_text.at(Arr_energy[j])), color[par+1], style[j], size[j]);
			}
			for(Int_t j=0; j<1; j++){

				//rebin_vec={0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.4};
				
				rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], "", CentBinMin, CentBinMax);
	    		Analysis -> SetParametrs(Arr_energy[j], file_1[1],file_out,harm+2,arr_par[par],"",EtaGap,prefix2);
				object.push_back(new Draw_Object);
				Analysis -> RatioParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i_%iGeV_%s",harm+2,Arr_energy[j],arr_par[par].Data()),Form("#font[42]{#scale[1.0]{Comb PID}}", Energy_scan_text.at(Arr_energy[j])), color[par], style[j+1], size[j]);
			}
			/*
			////////////////////////////////////
			
			for(Int_t j=0; j<1; j++){

				//rebin_vec={0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.4};
				
				rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], "", CentBinMin, CentBinMax);
	    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"",EtaGap,prefix3);
				object.push_back(new Draw_Object);
				Analysis -> RatioParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i_%iGeV_%s",harm+2,Arr_energy[j],arr_par[par].Data()),Form("#font[42]{#scale[1.0]{Comb PID 2}}", Energy_scan_text.at(Arr_energy[j])), colorall[2], styleall[2], size[j]);
			}
			
			for(Int_t j=0; j<1; j++){

				//rebin_vec={0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.4};
				
				rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], "", CentBinMin, CentBinMax);
	    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"",EtaGap,prefix4);
				object.push_back(new Draw_Object);
				Analysis -> RatioParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i_%iGeV_%s",harm+2,Arr_energy[j],arr_par[par].Data()),Form("#font[42]{#scale[1.0]{Comb PID 3}}", Energy_scan_text.at(Arr_energy[j])), colorall[3], styleall[3], size[j]);
			}
			
			for(Int_t j=0; j<1; j++){

				//rebin_vec={0.2,0.4,0.6,0.8,1.0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.4};
				
				rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], "", CentBinMin, CentBinMax);
	    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"",EtaGap,prefix5);
				object.push_back(new Draw_Object);
				Analysis -> RatioParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i_%iGeV_%s",harm+2,Arr_energy[j],arr_par[par].Data()),Form("#font[42]{#scale[1.0]{Comb PID 4}}", Energy_scan_text.at(Arr_energy[j])), colorall[4], styleall[4], size[j]);
			}
			*/
			////////////////////////////////////

			flow_ratio.push_back(object);
			object.clear();
			if(par==2){
				canvasRatio->SetTextInfo(Form("#font[42]{v_{%i}}", harm+2), Form("V%i",harm+2), 0.85*max_x_axis,0.99*coord_y_draw_ratio[2], 0.08, 0.0);
			}
			canvasRatio->SetNumberPad(Form("#font[42]{%s}", bukva[b2].Data()), "Pad", 0.04, 0.99*coord_y_draw_ratio[par], 0.08, 0.0);
			b2++;
		}
	}





	std::vector<Draw_Object *> object_sys;
	//ratio line
    for(Int_t i=0; i<flow_difference[0].size();i++){
    	object_sys.push_back(new Draw_Object);
    	Analysis -> RatioGraphEVAL(object_sys[i],flow_difference[0][i],flow_difference[0][1],"");
    	std::cout<<"good\n";
	}
	    //canvas setening
    Draw_Picture_new *canvas = new Draw_Picture_new;
    canvas->SetTLine(0.01, 0., max_x_axis, 0., 1, 2, 2, "Up");
    canvas->SetTLine(0.01, 0., max_x_axis, 0., 1, 2, 2, "Up");
    canvas->SetTLine(0.01, 0.9, max_x_axis, 0.9, 1, 2, 2, "Down");
    canvas->SetTLine(0.01, 1.1, max_x_axis, 1.1, 1, 2, 2, "Down");
    canvas->SetTLine(0.01, 1., max_x_axis, 1., flow_difference[0][1]->GetColor(), 2, 1, "Down");
	canvas->SetAxisToCanvsWithBottomPanel(0.01, max_x_axis, min_y_axis_dif[0]+0.002,  max_y_axis_dif[0], 0.76, 1.24);
	canvas->SetDrawObject(flow_difference[0],"Up");
	canvas->SetDrawObject(object_sys,"Down");
	//canvas->SetLegend(0.1,0.1,0.4,0.45,1);
	canvas->SetLegend(0.1,0.74,0.4,0.89,1);
	TCanvas *result33 = canvas->CanvasWithBottomPanel(Form("#font[42]{ #scale[1.2]{ #splitline{Au+Au #sqrt{s_{NN}} = 27 GeV}{Pion, centrality %i-%i %%}}}",cent[CentBinMax+1], cent[CentBinMin]), 0.36*max_x_axis, 0.88*coord_y_draw_dif[0],
													"#font[42]{p_{T} [GeV/c]}",
													Form("#font[42]{ v_{n}(%s) - v_{n}(%s) }",particle[GetNumberParticle("Pion","Pos")].Data(),particle[GetNumberParticle("Pion","Neg")].Data() ),
													"#font[42]{#frac{Data}{Def}}","");
	
	result33->SaveAs(Form("%s/picture/3_PID/Sys/test/run18/Aprez_%iGeVRun18_v%i_pt_Pion_%i_%i_new.png",outPath, 27,2, cent[CentBinMax+1], cent[CentBinMin]));


	std::vector<Draw_Object *> object_sys2;
	//ratio line
    for(Int_t i=0; i<flow_ratio[0].size();i++){
    	object_sys2.push_back(new Draw_Object);
    	Analysis -> RatioGraphEVAL(object_sys2[i],flow_ratio[0][i],flow_ratio[0][1],"");
    	std::cout<<"good\n";
	}
	    //canvas setening
    Draw_Picture_new *canvas32 = new Draw_Picture_new;
    canvas32->SetTLine(0.01, 1., max_x_axis, 1., 1, 2, 2, "Up");
    canvas32->SetTLine(0.01, 1., max_x_axis, 1., 1, 2, 2, "Up");
    canvas32->SetTLine(0.01, 0.9, max_x_axis, 0.9, 1, 2, 2, "Down");
    canvas32->SetTLine(0.01, 1.1, max_x_axis, 1.1, 1, 2, 2, "Down");
    canvas32->SetTLine(0.01, 1., max_x_axis, 1., flow_ratio[0][1]->GetColor(), 2, 1, "Down");
	canvas32->SetAxisToCanvsWithBottomPanel(0.01, max_x_axis, min_y_axis_ratio[0],  max_y_axis_ratio[0], 0.96, 1.04);
	canvas32->SetDrawObject(flow_ratio[0],"Up");
	canvas32->SetDrawObject(object_sys2,"Down");
	canvas32->SetLegend(0.1,0.78,0.4,0.89,1);
	//canvas32->SetLegend(0.1,0.1,0.4,0.45,1);
	TCanvas *result32 = canvas32->CanvasWithBottomPanel(Form("#font[42]{ #scale[1.2]{ #splitline{Au+Au #sqrt{s_{NN}} = 27 GeV}{Pion, centrality %i-%i %%}}}",cent[CentBinMax+1], cent[CentBinMin]), 0.36*max_x_axis, 0.99*coord_y_draw_ratio[0],
													"#font[42]{p_{T} [GeV/c]}",
													Form("#font[42]{#frac{ v_{n}(%s) }{ v_{n}(%s) } }", particle[GetNumberParticle("Pion","Pos")].Data(), particle[GetNumberParticle("Pion","Neg")].Data()),
													"#font[42]{#frac{Data}{Def}}","");
	
	result32->SaveAs(Form("%s/picture/3_PID/Sys/test/run18/Aprez_Rat_%iGeVRun18_v%i_pt_Pion_%i_%i_new.png",outPath, 27,2, cent[CentBinMax+1], cent[CentBinMin]));
	






	canvasDif->SetTextInfo(Form("#font[42]{ v_{n}(%s) - v_{n}(%s) }",particle[GetNumberParticle("Pion","Pos")].Data(),particle[GetNumberParticle("Pion","Neg")].Data() ),"Y",0.4, 0.4,0.4,90);
	//canvasDif->SetTextInfo(Form("#font[42]{ v_{n}(%s) - v_{n}(%s) }",particle[GetNumberParticle("Kaon","Pos")].Data(),particle[GetNumberParticle("Kaon","Neg")].Data() ),"Y",0.4, 0.4,0.4,90);
	//canvasDif->SetTextInfo(Form("#font[42]{ v_{n}(%s) - v_{n}(%s) }",particle[GetNumberParticle("Proton","Pos")].Data(),particle[GetNumberParticle("Proton","Neg")].Data() ),"Y",0.4, 0.3,0.35,90);
	canvasDif->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.3,0.35,0.8,0.0);
	canvasDif->SetTextInfo(Form("#font[42]{ #scale[1.0]{#splitline{Au+Au #sqrt{s_{NN}}=27 GeV}{Centrality %i-%i%%}}}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.3,0.9*coord_y_draw_dif[0],0.08,0.0);//0.8
	//canvasDif->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %i-%i %% }}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.6,0.95*coord_y_draw_dif[0],0.08,0.0);
	//canvasDif->SetTextInfo(Form("#font[42]{ #scale[1.0]{ STAR Preliminary }}"),"Prel",0.6,0.95*coord_y_draw_dif[2],0.08,0.0);
	canvasDif->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasDif->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasDif->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_dif[0], max_y_axis_dif[0], 0.06);
	//canvasDif->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_dif[1], max_y_axis_dif[1], 0.06);
	//canvasDif->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_dif[2], max_y_axis_dif[2], 0.06);
	canvasDif->SetDrawObject(flow_difference[0],"OnePad");
	canvasDif->SetLegend(0.07,0.08,0.68,0.25,1);
	TCanvas *result = canvasDif->CanvasNxM(1000,450, 1, 2, 0.4, 0.18, flow_difference,3,0,0);

	//result->SaveAs(Form("%s/picture/3_PID/Sys/test/NewComb_onlyProton27Run18_BESupdateDif_vn_%s_cent_%i_%i.png",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result->SaveAs(Form("%s/picture/3_PID/Sys/test/run18/Aprez_onlyPion27Run18_BESupdateDif_vn_%s_cent_%i_%i.png",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	//result->SaveAs(Form("%s/picture/3_PID/Sys/test/NewComb2_onlyPionRun18_27BESupdateDif_vn_%s_cent_%i_%i.png",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	//result->SaveAs(Form("%s/picture/3_PID/Sys/test/Comb_onlyKaon27BESupdateDif_vn_%s_cent_%i_%i.C",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result;
	
	canvasRatio->SetTextInfo(Form("#font[42]{#frac{ v_{n}(%s) }{ v_{n}(%s) } }", particle[GetNumberParticle("Pion","Pos")].Data(), particle[GetNumberParticle("Pion","Neg")].Data()),"Y",0.07,0.5,0.3,0);
	//canvasRatio->SetTextInfo(Form("#font[42]{#frac{ v_{n}(%s) }{ v_{n}(%s) } }", particle[GetNumberParticle("Kaon","Pos")].Data(), particle[GetNumberParticle("Kaon","Neg")].Data()),"Y",0.07,0.5,0.3,0);
	//canvasRatio->SetTextInfo(Form("#font[42]{#frac{ v_{n}(%s) }{ v_{n}(%s) } }", particle[GetNumberParticle("Proton","Pos")].Data(), particle[GetNumberParticle("Proton","Neg")].Data()),"Y",0.07,0.5,0.3,0);
	canvasRatio->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.3,0.35,0.8,0.0);
	canvasRatio->SetTextInfo(Form("#font[42]{ #scale[1.0]{#splitline{Au+Au #sqrt{s_{NN}}=27 GeV}{Centrality %i-%i%%}}}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.3,0.99*coord_y_draw_ratio[0],0.08,0.0);//0.98//0.8
	//canvasRatio->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %i-%i %% }}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.6,0.98*coord_y_draw_ratio[0],0.08,0.0);
	//canvasRatio->SetTextInfo(Form("#font[42]{ #scale[1.0]{ STAR Preliminary }}"),"Prel",0.6,0.99*coord_y_draw_ratio[2],0.08,0.0);
	canvasRatio->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasRatio->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_ratio[0], max_y_axis_ratio[0], 0.06);
	//canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_ratio[1], max_y_axis_ratio[1], 0.06);
	//canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_ratio[2], max_y_axis_ratio[2], 0.06);
	canvasRatio->SetDrawObject(flow_ratio[0],"OnePad");
	canvasRatio->SetLegend(0.04,0.05,0.68,0.25,1);
	TCanvas *result2 = canvasRatio->CanvasNxM(1000,450, 1, 2, 0.4, 0.18, flow_ratio,3,0,0);

	//result2->SaveAs(Form("%s/picture/3_PID/Sys/test/NewComb_onlyProton27Run18_BESupdateRatio_vn_%s_cent_%i_%i.png",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result2->SaveAs(Form("%s/picture/3_PID/Sys/test/run18/Aprez_onlyPion27Run18_BESupdateRatio_vn_%s_cent_%i_%i.png",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	//result2->SaveAs(Form("%s/picture/3_PID/Sys/test/NewComb2_onlyPion27Run18_BESupdateRatio_vn_%s_cent_%i_%i.png",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	//result2->SaveAs(Form("%s/picture/3_PID/Sys/test/Comb_onlyKaon27BESupdateRatio_vn_%s_cent_%i_%i.C",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	//result2->SaveAs(Form("%s/picture/3_PID/Sys/test/only27BESupdateRatio_vn_%s_cent_%i_%i_new.cpp",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result2;
	

	for(Int_t i=0; i<6; i++){
		file_1[i]->Close();
		delete file_1[i];
	}

}

/*
void ParAntVsEnergyMultyPad( Int_t CentBinMin, Int_t CentBinMax, Double_t ptMin, Double_t ptMax, TString EtaGap,  Double_t min_x_axis, Double_t max_x_axis, Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio, Double_t min_y_axis_dif, Double_t max_y_axis_dif){
	
	std::vector<VEC> rebin_vec = {{0.15,0.2,1.5,5.2},
								  {0.15,0.2,1.5,5.2},
								  {0.15,0.5,2.0,5.2}};
	Int_t point = 1; 
	TString prefix="TPCandTOF";

	const Int_t color[3]={2, 4, 1};
	const Int_t style[3]={23, 21, 22};
	Float_t size[3]={1.5,1.5,1.5};
	//const Int_t style_ant[]={27, 25, 26};
	std::vector<TString> sys_arr={"Eta01","Eta03","Eta05","Eta07"};
	std::vector<Double_t> flow_sys;

	Int_t Arr_energy[6]={11,14,19,27,39,62};
	Double_t x_energy[6]={11.5,14.5,19.6,27,39,62.4};
	TString arr_par[]={"Pion","Kaon","Proton"};

	TFile *file_1[6];
	for(Int_t i=0; i<6; i++){
		file_1[i] = new TFile(Form("%s/file_pid/Flow_%iGeV_PID_10binYesTofdca11_sys.root",outPath, Arr_energy[i]),"READ");
	}
	TFile *file_out = new TFile("OUT.root","RECREATE");	

	Analysis_Function *Analysis = new Analysis_Function;
	Draw_Object * object = new Draw_Object;
	std::vector<Draw_Object *> res;
	std::vector<std::vector<Draw_Object *>> flow ;
	Draw_Picture_new *canvasDif = new Draw_Picture_new;
	Int_t b=0;

	for(Int_t harm=0; harm<2; harm++){
		for(Int_t par=0; par<3; par++){

			res.push_back(new Draw_Object);
			
			for(Int_t j=0; j<6; j++){
				if(par==2){
					point=0;
				}else{
					point=1;
				}
	    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"",EtaGap,prefix);
				Analysis -> DifferentParAnt(object, CentBinMin, CentBinMax, rebin_vec[par]);
				res.back()->SetPointGraph(x_energy[j],object->GetPointYGraph(point),0.,object->GetPointYErrorGraph(point));
				res.back()->SetPointSysErrorGraph(object->GetPointYSysErrorGraph(point));
				object->ClearObject();
			}
	    	res.back() -> SetParametrsGraph(Form("data_v%i_%s",harm+2,arr_par[par].Data()),Form("#font[42]{#scale[1.1]{%s - %s}}", particle[GetNumberParticle(arr_par[par],"Pos")].Data(), particle[GetNumberParticle(arr_par[par],"Neg")].Data()), color[par], style[par], size[par]);
	    	
	    	/////////////////////// sys закоментить это надо
	    	for(Int_t s=0; s<sys_arr.size(); s++){
				for(Int_t j=0; j<6; j++){
					if(par==2){
						point=0;
					}else{
						point=1;
					}
		    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"",sys_arr[s],prefix);
					Analysis -> DifferentParAnt(object, CentBinMin, CentBinMax, rebin_vec[par]);
					flow_sys.push_back(object->GetPointYGraph(point));
					object->ClearObject();
				}
				res.back()->SetSysArrey(flow_sys);
				flow_sys.clear();
			}
			res.back()->SetSysErrors();
	    	///////////////////////
	    	

		}

		flow.push_back(res);
		res.clear();
		canvasDif->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.85*max_y_axis, 0.06, 0.0);
		canvasDif->SetTextInfo(Form("#font[42]{n = %i}", harm+2), Form("V%i",harm+2), 0.85*max_x_axis,0.85*max_y_axis, 0.06, 0.0);
			
		b++;
	}


	canvasDif->SetTextInfo(Form("#font[42]{ v_{n}(X) - v_{n}(#bar{X}) }"),"Y",0.4, 0.4,0.3,90);
	canvasDif->SetTextInfo("#font[42]{#sqrt{s_{NN}} [GeV]}","X",0.5,0.4,0.6,0.0);
	canvasDif->SetTextInfo(Form("#font[42]{ #scale[1.2]{ #splitline{Au+Au, %i-%i %%}{STAR Preliminary} }}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.2*max_x_axis,0.85*max_y_axis,0.06,0.0);
	canvasDif->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasDif->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasDif->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.06);
	canvasDif->SetDrawObject(flow[0],"OnePad");
	canvasDif->SetLegend(0.1,0.7,0.4,0.98,1);
	TCanvas *result = canvasDif->CanvasNxM(1000,650, 1, 2, 0.4, 0.15, flow,1,0,1);

	result->SaveAs(Form("%s/picture/3_PID/Sys/test/ParAntVsEnergy_vn_%s_cent_%i_%i_pt_%1.f_%1.f_k.pdf",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin], ptMin, ptMax));
	result->SaveAs(Form("%s/picture/3_PID/Sys/test/ParAntVsEnergy_vn_%s_cent_%i_%i_pt_%1.f_%1.f_k.png",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin], ptMin, ptMax));
	result->SaveAs(Form("%s/picture/3_PID/Sys/test/ParAntVsEnergy_vn_%s_cent_%i_%i_pt_%1.f_%1.f_k.C",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin], ptMin, ptMax));
	delete result;

	std::cout<<"good\n";

	for(Int_t i=0; i<6; i++){
		file_1[i]->Close();
		delete file_1[i];
	}
}
*/
void ParAntVsEnergyMultyPad( Int_t CentBinMin, Int_t CentBinMax, Double_t ptMin, Double_t ptMax, TString EtaGap,  Double_t min_x_axis, Double_t max_x_axis, Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio, Double_t min_y_axis_dif, Double_t max_y_axis_dif){
	
	std::vector<VEC> rebin_vec = {{0.15,0.2,1.5,5.2},
								  {0.15,0.2,1.5,5.2},
								  {0.15,0.5,2.0,5.2}};
	Int_t point = 1; 
	TString prefix="TPCandTOF";

	const Int_t color[3]={2, 4, 1};
	const Int_t style[3]={23, 21, 22};
	Float_t size[3]={1.5,1.5,1.5};
	//const Int_t style_ant[]={27, 25, 26};
	std::vector<TString> sys_arr={"Eta01","Eta03","Eta05","Eta07"};
	std::vector<Double_t> flow_sys;

	Int_t Arr_energy[6]={11,14,19,27,39,62};
	Double_t x_energy[6]={11.5,14.5,19.6,27,39,62.4};
	TString arr_par[]={"Pion","Kaon","Proton"};

	TFile *file_1[6];
	for(Int_t i=0; i<6; i++){
		file_1[i] = new TFile(Form("%s/file_pid/Flow_%iGeV_PID_10binYesTofdca11_sys.root",outPath, Arr_energy[i]),"READ");
	}
	TFile *file_out = new TFile("OUT.root","RECREATE");	

	Analysis_Function *Analysis = new Analysis_Function;
	Draw_Object * object = new Draw_Object;
	std::vector<Draw_Object *> res;
	std::vector<std::vector<Draw_Object *>> flow ;
	Draw_Picture_new *canvasDif = new Draw_Picture_new;
	Int_t b=0;

	TString pttext[]={"0.2 < p_{T} < 1.5 GeV/c","0.2 < p_{T} < 1.5 GeV/c","0.5 < p_{T} < 2.0 GeV/c"};

	for(Int_t harm=1; harm<2; harm++){
		for(Int_t par=0; par<3; par++){

			res.push_back(new Draw_Object);
			
			for(Int_t j=0; j<6; j++){
				if(par==2){
					point=0;
				}else{
					point=1;
				}
	    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"",EtaGap,prefix);
				Analysis -> DifferentParAnt(object, CentBinMin, CentBinMax, rebin_vec[par]);
				res.back()->SetPointGraph(x_energy[j],object->GetPointYGraph(point),0.,object->GetPointYErrorGraph(point));
				res.back()->SetPointSysErrorGraph(object->GetPointYSysErrorGraph(point));
				object->ClearObject();
			}
	    	res.back() -> SetParametrsGraph(Form("data_v%i_%s",harm+2,arr_par[par].Data()),Form("#font[42]{#scale[1.1]{%s - %s,  %s}}", particle[GetNumberParticle(arr_par[par],"Pos")].Data(), particle[GetNumberParticle(arr_par[par],"Neg")].Data(), pttext[par].Data()), color[par], style[par], size[par]);
	    	/*
	    	/////////////////////// sys
	    	for(Int_t s=0; s<sys_arr.size(); s++){
				for(Int_t j=0; j<6; j++){
					if(par==2){
						point=0;
					}else{
						point=1;
					}
		    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"",sys_arr[s],prefix);
					Analysis -> DifferentParAnt(object, CentBinMin, CentBinMax, rebin_vec[par]);
					flow_sys.push_back(object->GetPointYGraph(point));
					object->ClearObject();
				}
				res.back()->SetSysArrey(flow_sys);
				flow_sys.clear();
			}
			res.back()->SetSysErrors();
	    	///////////////////////
	    	*/

		}

		flow.push_back(res);
		res.clear();
		canvasDif->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.85*max_y_axis, 0.05, 0.0);
		canvasDif->SetTextInfo(Form("#font[42]{n = %i}", harm+2), Form("V%i",harm+2), 0.85*max_x_axis,0.85*max_y_axis, 0.05, 0.0);
			
		b++;
	}


	canvasDif->SetTextInfo(Form("#font[42]{ v_{3}(X) - v_{3}(#bar{X}) }"),"Y",0.6, 0.4,0.5,90);
	canvasDif->SetTextInfo("#font[42]{#sqrt{s_{NN}} [GeV]}","X",0.5,0.4,0.6,0.0);
	canvasDif->SetTextInfo(Form("#font[42]{ #scale[1.2]{ #splitline{Au+Au, %i-%i %%}{STAR Preliminary} }}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.2*max_x_axis,0.85*max_y_axis,0.05,0.0);
	canvasDif->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasDif->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasDif->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.05);
	canvasDif->SetDrawObject(flow[0],"OnePad");
	canvasDif->SetLegend(0.35,0.65,0.94,0.85,1);
	TCanvas *result = canvasDif->CanvasNxM(1000,900, 1, 1, 0.3, 0.14, flow,1,0,0);

	result->SaveAs(Form("%s/picture/3_PID/Sys/test/ParAntVsEnergy_v3/ver3/ParAntVsEnergy_v3_%s_cent_%i_%i.pdf",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	//result->SaveAs(Form("%s/picture/3_PID/Sys/test/ParAntVsEnergy_v3_%s_cent_%i_%i_pt_%1.f_%1.f.pdf",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin], ptMin, ptMax));
	result->SaveAs(Form("%s/picture/3_PID/Sys/test/ParAntVsEnergy_v3/ver3/ParAntVsEnergy_v3_%s_cent_%i_%i.png",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	//result->SaveAs(Form("%s/picture/3_PID/Sys/test/ParAntVsEnergy_v3_%s_cent_%i_%i_pt_%1.f_%1.f.png",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin], ptMin, ptMax));
	result->SaveAs(Form("%s/picture/3_PID/Sys/test/ParAntVsEnergy_v3/ver3/ParAntVsEnergy_v3_%s_cent_%i_%i.C",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	//result->SaveAs(Form("%s/picture/3_PID/Sys/test/ParAntVsEnergy_v3_%s_cent_%i_%i_pt_%1.f_%1.f.C",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin], ptMin, ptMax));
	delete result;

	std::cout<<"good\n";

	for(Int_t i=0; i<6; i++){
		file_1[i]->Close();
		delete file_1[i];
	}
}


int GetNumberParticle(TString PID, TString charge){
	
	Int_t nPID;
    
    if( strncmp(PID, "_proton",7)==0 || strncmp(PID, "Proton",6)==0){
        if(strncmp(charge, "Pos",3)==0){
          nPID=4;
        }
        else{
          nPID=5;
        }
        if(strncmp(charge, "",1)==0){
          nPID=8;
        }
    }
    if( strncmp(PID, "_kaon",5)==0 || strncmp(PID, "Kaon",4)==0){
        if(strncmp(charge, "Pos",3)==0){
          nPID=2;
        }
        else{
          nPID=3;
        }
        if(strncmp(charge, "",1)==0){
          nPID=7;
        }
    }
    if( strncmp(PID, "_pion",5)==0 || strncmp(PID, "Pion",4)==0){
        if(strncmp(charge, "Pos",3)==0){
          nPID=0;
        }
        else{
          nPID=1;
        }
        if(strncmp(charge, "",1)==0){
          nPID=6;
        }
    }
    if( strncmp(PID, "Hadrons",7)==0 ){
        nPID=9;
    }
    if( strncmp(PID, "",1)==0){
        nPID=10;
    }
    if( strncmp(PID, "PID",3)==0){
        nPID=11;
    }
    return nPID;
}

void DrawFlowVsEnergyBES(){
	
	//FlowVsEnergyBESupdate( 0.2, 3.2, "Eta15","Hadrons", "", 0.0, 0.115, 75, 0.5, 1.15);
	//FlowVsEnergyBESupdate( 0.2, 2.0, "Eta15","Hadrons", "", 0.0, 0.115, 4, 75, 0.5, 1.15);

}

void DrawFlowDifferentEtaGap(TString mod, Int_t Energy){
	if(strncmp(mod, "Hadrons",7)==0){
		FlowDifferentEtaGap(Energy, 2, 2, 6,"Hadrons","",0.28, 3.5,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 7, 8,"Hadrons","",0.28, 3.5,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 2, 8,"Hadrons","",0.28, 3.5,0.89,1.11);
		FlowDifferentEtaGap(Energy, 3, 2, 6,"Hadrons","",0.15, 3.5,0.8,1.2);
		FlowDifferentEtaGap(Energy, 3, 7, 8,"Hadrons","",0.15, 3.5,0.8,1.2);
		FlowDifferentEtaGap(Energy, 3, 2, 8,"Hadrons","",0.15, 3.5,0.8,1.2);
		//FlowDifferentEtaGap(Energy, 3, 2, 8,"Hadrons","",0.15, 3.5,0.8,1.2);
	}else{
		FlowDifferentEtaGap(Energy, 2, 2, 6,"Pion","Pos",0.18, 2.87,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 7, 8,"Pion","Pos",0.13, 2.87,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 2, 8,"Pion","Pos",0.19, 2.87,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 2, 6,"Pion","Neg",0.18, 2.87,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 7, 8,"Pion","Neg",0.13, 2.87,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 2, 8,"Pion","Neg",0.19, 2.87,0.89,1.11);

		FlowDifferentEtaGap(Energy, 2, 2, 6,"Kaon","Pos",0.18, 2.87,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 7, 8,"Kaon","Pos",0.13, 2.87,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 2, 8,"Kaon","Pos",0.19, 2.87,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 2, 6,"Kaon","Neg",0.18, 2.87,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 7, 8,"Kaon","Neg",0.13, 2.87,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 2, 8,"Kaon","Neg",0.19, 2.87,0.89,1.11);

		FlowDifferentEtaGap(Energy, 2, 2, 6,"Proton","Pos",0.18, 2.87,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 7, 8,"Proton","Pos",0.18, 2.87,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 2, 8,"Proton","Pos",0.21, 2.87,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 2, 6,"Proton","Neg",0.18, 2.87,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 7, 8,"Proton","Neg",0.15, 2.87,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 2, 8,"Proton","Neg",0.21, 2.87,0.89,1.11);

		FlowDifferentEtaGap(Energy, 3, 2, 6,"Pion","Pos",0.15, 2.87,0.84,1.16);
		FlowDifferentEtaGap(Energy, 3, 2, 8,"Pion","Pos",0.18, 2.87,0.84,1.16);
		FlowDifferentEtaGap(Energy, 3, 2, 6,"Pion","Neg",0.15, 2.87,0.84,1.16);
		FlowDifferentEtaGap(Energy, 3, 2, 8,"Pion","Neg",0.18, 2.87,0.84,1.16);

		FlowDifferentEtaGap(Energy, 3, 2, 6,"Kaon","Pos",0.15, 2.87,0.84,1.16);
		FlowDifferentEtaGap(Energy, 3, 2, 8,"Kaon","Pos",0.18, 2.87,0.84,1.16);
		FlowDifferentEtaGap(Energy, 3, 2, 6,"Kaon","Neg",0.15, 2.87,0.84,1.16);
		FlowDifferentEtaGap(Energy, 3, 2, 8,"Kaon","Neg",0.18, 2.87,0.84,1.16);

		FlowDifferentEtaGap(Energy, 3, 2, 6,"Proton","Pos",0.18, 2.87,0.84,1.16);
		FlowDifferentEtaGap(Energy, 3, 2, 8,"Proton","Pos",0.18, 2.87,0.84,1.16);
		FlowDifferentEtaGap(Energy, 3, 2, 6,"Proton","Neg",0.18, 2.87,0.84,1.16);
		FlowDifferentEtaGap(Energy, 3, 2, 8,"Proton","Neg",0.18, 2.87,0.84,1.16);
	
	}
}

void DrawFlowDifferentMethodPID(TString mod, Int_t Energy){

	FlowDifferentMethodPID(Energy, 2, 2, 6,"Eta01","Pion","Pos",0.18, 1.37,0.79,1.21);
	FlowDifferentMethodPID(Energy, 2, 7, 8,"Eta01","Pion","Pos",0.157, 1.37,0.79,1.21);
	FlowDifferentMethodPID(Energy, 2, 2, 8,"Eta01","Pion","Pos",0.157, 1.37,0.79,1.21);
	FlowDifferentMethodPID(Energy, 2, 2, 6,"Eta01","Pion","Neg",0.18, 1.37,0.79,1.21);
	FlowDifferentMethodPID(Energy, 2, 7, 8,"Eta01","Pion","Neg",0.157, 1.37,0.79,1.21);
	FlowDifferentMethodPID(Energy, 2, 2, 8,"Eta01","Pion","Neg",0.157, 1.37,0.79,1.21);

	FlowDifferentMethodPID(Energy, 2, 2, 6,"Eta01","Kaon","Pos",0.18, 1.37,0.79,1.21);
	FlowDifferentMethodPID(Energy, 2, 7, 8,"Eta01","Kaon","Pos",0.157, 1.37,0.79,1.21);
	FlowDifferentMethodPID(Energy, 2, 2, 8,"Eta01","Kaon","Pos",0.157, 1.37,0.72,1.28);
	FlowDifferentMethodPID(Energy, 2, 2, 6,"Eta01","Kaon","Neg",0.18, 1.37,0.79,1.21);
	FlowDifferentMethodPID(Energy, 2, 7, 8,"Eta01","Kaon","Neg",0.157, 1.37,0.79,1.21);
	FlowDifferentMethodPID(Energy, 2, 2, 8,"Eta01","Kaon","Neg",0.157, 1.37,0.72,1.28);

/*

	FlowDifferentMethodPID(Energy, 2, 2, 6,"Eta01","Proton","Pos",0.18, 1.37,0.79,1.21);
	FlowDifferentMethodPID(Energy, 2, 7, 8,"Eta01","Proton","Pos",0.157, 1.37,0.79,1.21);
	FlowDifferentMethodPID(Energy, 2, 2, 8,"Eta01","Proton","Pos",0.157, 1.37,0.72,1.28);
	FlowDifferentMethodPID(Energy, 2, 2, 6,"Eta01","Proton","Neg",0.18, 1.37,0.79,1.21);
	FlowDifferentMethodPID(Energy, 2, 7, 8,"Eta01","Proton","Neg",0.157, 1.37,0.79,1.21);
	FlowDifferentMethodPID(Energy, 2, 2, 8,"Eta01","Proton","Neg",0.157, 1.37,0.72,1.28);
*/	
	FlowDifferentMethodPID(Energy, 3, 2, 6,"Eta01","Pion","Pos",0.093, 1.37,0.79,1.21);
	FlowDifferentMethodPID(Energy, 3, 2, 8,"Eta01","Pion","Pos",0.093, 1.37,0.72,1.28);
	FlowDifferentMethodPID(Energy, 3, 7, 8,"Eta01","Pion","Pos",0.093, 1.37,0.72,1.28);
	FlowDifferentMethodPID(Energy, 3, 2, 6,"Eta01","Pion","Neg",0.093, 1.37,0.79,1.21);
	FlowDifferentMethodPID(Energy, 3, 2, 8,"Eta01","Pion","Neg",0.093, 1.37,0.72,1.28);
	FlowDifferentMethodPID(Energy, 3, 7, 8,"Eta01","Pion","Neg",0.093, 1.37,0.72,1.28);

	FlowDifferentMethodPID(Energy, 3, 2, 6,"Eta01","Kaon","Pos",0.093, 1.37,0.74,1.26);
	FlowDifferentMethodPID(Energy, 3, 2, 8,"Eta01","Kaon","Pos",0.093, 1.37,0.72,1.28);
	FlowDifferentMethodPID(Energy, 3, 7, 8,"Eta01","Kaon","Pos",0.093, 1.37,0.72,1.28);
	FlowDifferentMethodPID(Energy, 3, 2, 6,"Eta01","Kaon","Neg",0.093, 1.37,0.74,1.26);
	FlowDifferentMethodPID(Energy, 3, 2, 8,"Eta01","Kaon","Neg",0.093, 1.37,0.72,1.28);
	FlowDifferentMethodPID(Energy, 3, 7, 8,"Eta01","Kaon","Neg",0.093, 1.37,0.72,1.28);
/*
	FlowDifferentMethodPID(Energy, 3, 2, 6,"Eta01","Proton","Pos",0.1, 1.37,0.74,1.26);
	FlowDifferentMethodPID(Energy, 3, 2, 8,"Eta01","Proton","Pos",0.1, 1.37,0.72,1.28);
	FlowDifferentMethodPID(Energy, 3, 2, 6,"Eta01","Proton","Neg",0.1, 1.37,0.74,1.26);
	FlowDifferentMethodPID(Energy, 3, 2, 8,"Eta01","Proton","Neg",0.1, 1.37,0.72,1.28);
	*/

}

void DrawFlowVsRapidity(Int_t Energy){
 	
	//FlowVsRapidity(Energy,2,3,8,0.2,1.5,"Eta01","Pion", "",0.01,0.1);
	//FlowVsRapidity(Energy,2,3,8,0.2,1.5,"Eta01","Kaon", "",0.01,0.12);
	//FlowVsRapidity(Energy,2,3,8,0.5,1.5,"Eta01","Proton", "",0.01,0.15);

	//FlowVsRapidity(Energy,3,3,8,0.2,1.5,"Eta01","Pion", "",0.007,0.02);
	//FlowVsRapidity(Energy,3,3,8,0.2,1.5,"Eta01","Kaon", "",0.009,0.027);
	//FlowVsRapidity(Energy,3,3,8,0.5,1.5,"Eta01","Proton", "",0.007,0.03);
	
	//FlowVsRapidity(Energy,2,3,8,0.3,1.5,"Eta01","Pion", "",0.01,0.1);
	//FlowVsRapidity(Energy,2,3,8,0.3,1.5,"Eta01","Pion", "",0.01,0.1);
	 
	FlowVsRapidity(Energy,2,3,8,0.2,3.2,"Eta15","Hadrons", "",0.01,0.1);
 	//FlowVsRapiditySymmetry(Energy,2,0.2,2.,"Eta15","Hadrons", "",0.01,0.09);
 	//FlowVsRapiditySymmetry(Energy,2,0.2,3.2,"Eta15","Hadrons", "",0.01,0.09);
 	//FlowVsRapiditySymmetry(Energy,3,0.2,2.,"Eta15","Hadrons", "",0.007,0.023);
 	//FlowVsRapiditySymmetry(Energy,3,0.2,3.2,"Eta15","Hadrons", "",0.007,0.023);
 	//FlowVsRapidity(Energy,3,3,8,0.2,3.2,"Eta15","Hadrons", "",0.008,0.027);
}

void DrawFlowVsPtForBES(TString mod){
	if(strncmp(mod, "Hadrons",7)==0){
		FlowVsPtForBES( 3, 4, 6, 0.2, 3.2,"Eta15","Hadrons", "", 0.15, 7. ,3.5, 0.6, 1.1, 0.7, 1.3);
		FlowVsPtForBES( 2, 4, 6, 0.2, 3.2,"Eta15","Hadrons", "", 0.28, 5. ,3.5, 0.6, 1.1, 0.7, 1.3);
		FlowVsPtForBES( 3, 7, 8, 0.2, 3.2,"Eta15","Hadrons", "", 0.15, 7. ,3.5, 0.5, 1.1, 0.7, 1.3);
		FlowVsPtForBES( 2, 7, 8, 0.2, 3.2,"Eta15","Hadrons", "", 0.28, 5. ,3.5, 0.5, 1.1, 0.7, 1.3);
		FlowVsPtForBES( 3, 2, 6, 0.2, 3.2,"Eta15","Hadrons", "", 0.15, 7. ,3.5, 0.6, 1.1, 0.7, 1.3);
		FlowVsPtForBES( 2, 2, 6, 0.2, 3.2,"Eta15","Hadrons", "", 0.28, 5. ,3.5, 0.6, 1.1, 0.7, 1.3);
		FlowVsPtForBES( 3, 3, 6, 0.2, 3.2,"Eta15","Hadrons", "", 0.15, 7. ,3.5, 0.6, 1.1, 0.7, 1.3);
		FlowVsPtForBES( 2, 3, 6, 0.2, 3.2,"Eta15","Hadrons", "", 0.28, 5. ,3.5, 0.6, 1.1, 0.7, 1.3);
		FlowVsPtForBES( 3, 2, 8, 0.2, 3.2,"Eta15","Hadrons", "", 0.15, 7. ,3.5, 0.6, 1.1, 0.7, 1.3);
		FlowVsPtForBES( 2, 2, 8, 0.2, 3.2,"Eta15","Hadrons", "", 0.28, 5. ,3.5, 0.6, 1.1, 0.7, 1.3);

	}else{/*
		FlowVsPtForBES( 2, 2, 8, "Eta01","Kaon", "", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 2, 8, "Eta01","Pion", "", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 2, 8, "Eta01","Proton", "", 0.28,4.5,0.6, 1.1);

		FlowVsPtForBES( 2, 7, 8, "Eta01","Kaon", "", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 7, 8, "Eta01","Pion", "", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 7, 8, "Eta01","Proton", "", 0.28,4.5,0.6, 1.1);

		FlowVsPtForBES( 2, 2, 6, "Eta01","Kaon", "", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 2, 6, "Eta01","Pion", "", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 2, 6, "Eta01","Proton", "", 0.28,4.5,0.6, 1.1);

		FlowVsPtForBES( 2, 4, 6, "Eta01","Kaon", "", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 4, 6, "Eta01","Pion", "", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 4, 6, "Eta01","Proton", "", 0.28,4.5,0.6, 1.1);

		FlowVsPtForBES( 3, 2, 8, "Eta01","Kaon", "", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 2, 8, "Eta01","Pion", "", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 2, 8, "Eta01","Proton", "", 0.15,3.5,0.6, 1.1);

		FlowVsPtForBES( 3, 7, 8, "Eta01","Kaon", "", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 7, 8, "Eta01","Pion", "", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 7, 8, "Eta01","Proton", "", 0.15,3.5,0.6, 1.1);

		FlowVsPtForBES( 3, 2, 6, "Eta01","Kaon", "", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 2, 6, "Eta01","Pion", "", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 2, 6, "Eta01","Proton", "", 0.15,3.5,0.6, 1.1);

		FlowVsPtForBES( 3, 4, 6, "Eta01","Kaon", "", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 4, 6, "Eta01","Pion", "", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 4, 6, "Eta01","Proton", "", 0.15,3.5,0.6, 1.1);


		FlowVsPtForBES( 2, 2, 8, "Eta01","Kaon", "Pos", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 2, 8, "Eta01","Pion", "Pos", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 2, 8, "Eta01","Proton", "Pos", 0.28,4.5,0.6, 1.1);

		FlowVsPtForBES( 2, 7, 8, "Eta01","Kaon", "Pos", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 7, 8, "Eta01","Pion", "Pos", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 7, 8, "Eta01","Proton", "Pos", 0.28,4.5,0.6, 1.1);

		FlowVsPtForBES( 2, 2, 6, "Eta01","Kaon", "Pos", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 2, 6, "Eta01","Pion", "Pos", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 2, 6, "Eta01","Proton", "Pos", 0.28,4.5,0.6, 1.1);

		FlowVsPtForBES( 2, 4, 6, "Eta01","Kaon", "Pos", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 4, 6, "Eta01","Pion", "Pos", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 4, 6, "Eta01","Proton", "Pos", 0.28,4.5,0.6, 1.1);

		FlowVsPtForBES( 3, 2, 8, "Eta01","Kaon", "Pos", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 2, 8, "Eta01","Pion", "Pos", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 2, 8, "Eta01","Proton", "Pos", 0.15,3.5,0.6, 1.1);

		FlowVsPtForBES( 3, 7, 8, "Eta01","Kaon", "Pos", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 7, 8, "Eta01","Pion", "Pos", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 7, 8, "Eta01","Proton", "Pos", 0.15,3.5,0.6, 1.1);

		FlowVsPtForBES( 3, 2, 6, "Eta01","Kaon", "Pos", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 2, 6, "Eta01","Pion", "Pos", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 2, 6, "Eta01","Proton", "Pos", 0.15,3.5,0.6, 1.1);

		FlowVsPtForBES( 3, 4, 6, "Eta01","Kaon", "Pos", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 4, 6, "Eta01","Pion", "Pos", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 4, 6, "Eta01","Proton", "Pos", 0.15,3.5,0.6, 1.1);


		FlowVsPtForBES( 2, 2, 8, "Eta01","Kaon", "Neg", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 2, 8, "Eta01","Pion", "Neg", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 2, 8, "Eta01","Proton", "Neg", 0.28,4.5,0.6, 1.1);

		FlowVsPtForBES( 2, 7, 8, "Eta01","Kaon", "Neg", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 7, 8, "Eta01","Pion", "Neg", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 7, 8, "Eta01","Proton", "Neg", 0.28,4.5,0.6, 1.1);

		FlowVsPtForBES( 2, 2, 6, "Eta01","Kaon", "Neg", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 2, 6, "Eta01","Pion", "Neg", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 2, 6, "Eta01","Proton", "Neg", 0.28,4.5,0.6, 1.1);

		FlowVsPtForBES( 2, 4, 6, "Eta01","Kaon", "Neg", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 4, 6, "Eta01","Pion", "Neg", 0.28,4.5,0.6, 1.1);
		FlowVsPtForBES( 2, 4, 6, "Eta01","Proton", "Neg", 0.28,4.5,0.6, 1.1);

		FlowVsPtForBES( 3, 2, 8, "Eta01","Kaon", "Neg", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 2, 8, "Eta01","Pion", "Neg", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 2, 8, "Eta01","Proton", "Neg", 0.15,3.5,0.6, 1.1);

		FlowVsPtForBES( 3, 7, 8, "Eta01","Kaon", "Neg", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 7, 8, "Eta01","Pion", "Neg", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 7, 8, "Eta01","Proton", "Neg", 0.15,3.5,0.6, 1.1);

		FlowVsPtForBES( 3, 2, 6, "Eta01","Kaon", "Neg", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 2, 6, "Eta01","Pion", "Neg", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 2, 6, "Eta01","Proton", "Neg", 0.15,3.5,0.6, 1.1);

		FlowVsPtForBES( 3, 4, 6, "Eta01","Kaon", "Neg", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 4, 6, "Eta01","Pion", "Neg", 0.15,3.5,0.6, 1.1);
		FlowVsPtForBES( 3, 4, 6, "Eta01","Proton", "Neg", 0.15,3.5,0.6, 1.1);
		*/
	}

}

void makeplotstyle(){
  
    TStyle *mystyle = new TStyle("PlottingInStyle", "Style for Summary Plots");
    mystyle->SetPalette(1);
    mystyle->SetCanvasColor(10);
    mystyle->SetHistFillColor(10);
    mystyle->SetHistFillStyle(0);
    mystyle->SetOptTitle(0);
    mystyle->SetOptStat(0);
    mystyle->SetCanvasBorderMode(0);//removes the yellow frame around the canvas
    mystyle->SetPadLeftMargin(0.16);
    mystyle->SetPadBottomMargin(0.15);
    mystyle->SetPadTickX(1);
    mystyle->SetPadTickY(1);
    mystyle->SetAxisColor(1, "X");
    mystyle->SetAxisColor(1, "Y");
    mystyle->SetLabelColor(1, "X");
    mystyle->SetLabelColor(1, "Y");
    mystyle->SetTickLength(0.03, "X");
    mystyle->SetTickLength(0.03, "Y");
    mystyle->SetTitleXSize(0.05);
    mystyle->SetTitleYSize(0.05);
    mystyle->SetNdivisions(508, "X");
    mystyle->SetNdivisions(505, "Y");
    mystyle->SetTitleXOffset(1.2);
    mystyle->SetTitleYOffset(1.4);
    mystyle->SetLabelOffset(0.02, "X");
    mystyle->SetLabelOffset(0.02, "Y");
    mystyle->SetLabelSize(0.05, "X");
    mystyle->SetLabelSize(0.05, "Y");
    //mystyle->SetGridx();

  	TFile f("style.root", "RECREATE");
  	f.cd();
  	mystyle->Write();
  	f.Close();
}

VEC Rebin2K(Int_t harmonic, Int_t Energy, TString PID, TString charge, Int_t CentBinMin, Int_t CentBinMax){
	VEC rebin_vector;
	
	rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.2, 10.0};
	return rebin_vector;


}


VEC Rebin2(Int_t harmonic, Int_t Energy, TString PID, TString charge, Int_t CentBinMin, Int_t CentBinMax){

	VEC rebin_vector;

	/*
	if(strncmp(PID, "Hadrons",7)==0 ){
		rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 10.0};
		//rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 10.0};
	}

	*/
	///////////////////////////////////////////////////////
	/*
	if(harmonic==2){

		if(Energy==19 || Energy==27 || Energy==39 || Energy == 62){
			if(CentBinMin==2 && CentBinMax==8){ //0-60
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.2, 4.0, 10.0};
			}
			if(CentBinMin==7 && CentBinMax==8){ //0-10
				if( strncmp(PID, "Kaon",4)==0 ){
					rebin_vector = {0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.8,5.0};
				}else{
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.8, 2.2, 2.6, 3.2 ,5.0};
				}
			}
			if(CentBinMin==2 && CentBinMax==6){ //10-60
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 2.6, 3.4, 5.0};
			}
					if(CentBinMin==4 && CentBinMax==6){ //10-60
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 2.6, 3.4, 5.0};
			}
		}

		if(Energy==14){
			if(CentBinMin==2 && CentBinMax==8){ //0-60
				if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4, 5.0};
				}
			}
			if(CentBinMin==7 && CentBinMax==8){ //0-10
				if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
					rebin_vector = {0.2, 0.6, 1.0, 1.4, 2.0, 2.6, 5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.8, 1.6, 2.4, 10.0};
				}
			}
			if(CentBinMin==2 && CentBinMax==6){ //10-60
				if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 2.6, 5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.8, 2.2, 5.0};
				}
			}
					if(CentBinMin==4 && CentBinMax==6){ //10-60
				if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 2.6, 5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.8, 2.2, 5.0};
				}
			}
		}

		if(Energy==11){
			if(CentBinMin==2 && CentBinMax==8){ //0-60
				if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4, 5.0};
				}
			}
			if(CentBinMin==7 && CentBinMax==8){ //0-10
				if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
					rebin_vector = {0.2, 0.6, 1.0, 1.4, 2.0, 2.6, 5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.8, 1.6, 2.4, 10.0};
				}
			}
			if(CentBinMin==2 && CentBinMax==6){ //10-60
				if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 2.6, 5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.8, 2.2, 5.0};
				}
			}
					if(CentBinMin==4 && CentBinMax==6){ //10-60
				if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 2.6, 5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.8, 2.2, 5.0};
				}
			}
		}

		if(Energy==7){
			if(CentBinMin==2 && CentBinMax==8){ //0-60
				if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
					rebin_vector = {0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.8,5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 10.0};
				}
			}
			if(CentBinMin==7 && CentBinMax==8){ //0-10
				if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
					//rebin_vector = {0.2, 0.4, 0.6, 1.0, 1.4, 1.8, 5.0};
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.4, 1.8, 2.2,5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.6, 1.0, 2.0, 10.0};
				}
			}
			if(CentBinMin==2 && CentBinMax==6){ //10-60
				if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
					rebin_vector = {0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.8,5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.6, 1.0, 1.4, 2.0, 10.0};
				}
			}
			if(CentBinMin==4 && CentBinMax==6){ //10-40
				if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
					rebin_vector = {0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.8,5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.6, 1.0, 1.4, 2.0, 10.0};
				}
			}
		}
		
	}
	///// потом закоментить
		*/

	if(harmonic==2){

		if(Energy==19 || Energy==27 || Energy==39 || Energy == 62){
			if( (CentBinMin==2 && CentBinMax==6) || (CentBinMin==2 && CentBinMax==8)){ //10-60
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 2.6, 3.2, 5.0};
			}
			if((CentBinMin==4 && CentBinMax==6) || (CentBinMin==2 && CentBinMax==6)){ //10-60
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 2.6, 3.2, 5.0};
			}
		}

		if((CentBinMin==3 && CentBinMax==8) || (CentBinMin==2 && CentBinMax==8) || (CentBinMin==4 && CentBinMax==8) || (CentBinMin==2 && CentBinMax==6) || (CentBinMin==7 && CentBinMax==8)){ //0-60
			
			if(Energy>25){
				if(strncmp(PID, "Pion",4)==0 ){
					rebin_vector = {0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4, 3.2, 5.0};
					//rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4, 3.2, 5.0};
				}
				if(strncmp(PID, "Kaon",4)==0 ){
					///нужно///rebin_vector = {0.4, 0.6, 0.8, 1.0, 1.4}; // до 1.8
					//rebin_vector = {0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 5.0}; старое
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 2.8, 5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					//rebin_vector = {0.2, 0.5, 0.6, 0.8, 1.2, 1.8, 2.2, 2.8, 5.0};
					rebin_vector = {0.2, 0.5, 0.6, 0.8, 1.2, 1.8, 2.2, 2.8, 5.0};
				}
			}

			if(Energy==19){
				if(strncmp(PID, "Pion",4)==0 ){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.6, 2.0, 2.6, 3.2, 5.0};;
				}
				if(strncmp(PID, "Kaon",4)==0 ){
					//rebin_vector = {0.2, 0.4, 0.6, 1.0, 1.4, 1.8, 2.4, 3.0, 5.0};;
					rebin_vector = {0.2, 0.4, 0.6, 1.0, 1.4, 1.6, 5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.5, 0.6, 1.0, 1.4, 1.8, 2.2, 2.8, 10.0};
				}
			}

			if(Energy==14){
				if(strncmp(PID, "Pion",4)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.6, 5.0};
				}
				if(strncmp(PID, "Kaon",4)==0){
					//rebin_vector = {0.2, 0.6, 1.0, 1.4, 1.8, 2.4, 5.0};
					rebin_vector = {0.2, 0.6, 1.0, 1.4, 1.6, 5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.5, 0.8, 1.2, 1.8, 2.4, 10.0};
				}
			}

			if(Energy==11 || Energy==7){
				if(strncmp(PID, "Pion",4)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.6, 2.2, 2.8, 5.0};
				}
				if(strncmp(PID, "Kaon",4)==0){
					rebin_vector = {0.2, 0.6, 1.0, 1.4, 1.6, 5.0};
					//rebin_vector = {0.2, 0.6, 1.0, 1.4, 2.0, 2.8, 5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.5, 0.8, 1.4, 2.6, 10.0};
				}
			}

		}
	}

	if(harmonic==3){

		if( (CentBinMin==3 && CentBinMax==8) || (CentBinMin==2 && CentBinMax==8) || (CentBinMin==2 && CentBinMax==6) || (CentBinMin==4 && CentBinMax==8)){ //0-60
			
			if(Energy>25){
				if(strncmp(PID, "Pion",4)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.6, 2.2, 2.6, 5.0};
					//rebin_vector = {0.4, 0.6, 0.8, 1.2, 1.6, 2.2, 2.6, 5.0};
				}
				if(strncmp(PID, "Kaon",4)==0){
					//rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.4, 1.8, 2.2, 2.8, 5.0};
					
					//rebin_vector = {0.2, 0.4, 0.6, 1.0, 1.4, 1.8, 2.2, 2.8, 5.0};
					//rebin_vector = {0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.8, 5.0};
					rebin_vector = {0.2, 0.6, 1.0, 1.4, 1.6, 5.0};
					//rebin_vector = { 0.6, 1.0, 1.4, 1.8, 2.2, 2.8, 5.0};
				
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.5, 0.8, 1.2, 1.8, 2.6, 10.0};
				}
			}


			if(Energy==19){
				if(strncmp(PID, "Pion",4)==0 ){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.6, 2.4, 5.0};
				}
				if(strncmp(PID, "Kaon",4)==0 ){
					//rebin_vector = {0.2, 0.4, 0.8, 1.2, 1.6, 2.4, 5.0};
					rebin_vector = {0.2, 0.4, 0.8, 1.2, 1.6, 5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.5, 0.8, 1.2, 1.8, 2.6, 10.0};
				}
			}

			if(Energy==14){
				if(strncmp(PID, "Pion",4)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.6, 2.0, 5.0};
				}
				if(strncmp(PID, "Kaon",4)==0){
					//rebin_vector = {0.2, 0.6, 1.0, 1.4, 2.0, 5.0};
					rebin_vector = {0.2, 0.6, 1.0, 1.4, 1.6, 5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.5, 0.8, 1.6, 2.6, 10.0};
				}
			}

			if(Energy==11 || Energy==7){
				if(strncmp(PID, "Pion",4)==0){
					rebin_vector = {0.2, 0.6, 1.0, 1.6, 2.2, 5.0};
				}
				if(strncmp(PID, "Kaon",4)==0){
					//rebin_vector = {0.2, 0.6, 1.0, 1.6, 2.2, 5.0};
					rebin_vector = {0.2, 0.6, 1.0, 1.6, 5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.5, 0.8, 1.4, 2.6, 10.0};
				}
			}




			//////////////////////

			if(Energy==27 || Energy == 39 /*|| Energy==62*/){
				if(strncmp(PID, "Pion",4)==0){
					//rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 5.0};
					//rebin_vector = { 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2,2.6, 5.0};
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 5.0};
				}
				if(strncmp(PID, "Kaon",4)==0){
					/////rebin_vector = { 0.4, 0.6, 0.8, 1.2, 1.4, 1.8, 5.0};
					//rebin_vector = { 0.4, 0.6, 0.8, 1.2, 1.4, 1.8, 2.4, 5.0};
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.4, 1.8, 2.2, 2.8, 5.0};
					//////////rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.4, 1.6, 5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = { 0.5, 0.8, 1.0, 1.2, 1.4, 1.8, 2.4, 10.0};
				}
			}

			//////////////

		}


		if(CentBinMin==7 && CentBinMax==8){ //0-10
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.6, 1.0, 1.4, 1.8, 2.4, 5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.5, 1.0, 1.4, 2.2, 2.8, 10.0};
			}
		}
		if(CentBinMin==2 && CentBinMax==6){ //10-60
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.8,5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.5, 1.0, 1.4, 2.0, 10.0};
			}
			if(Energy==27){
				if(strncmp(PID, "Pion",4)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 5.0};
				}
				if(strncmp(PID, "Kaon",4)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.4, 1.8, 2.2, 2.6, 5.0};				
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = { 0.5, 0.8, 1.2, 1.4, 1.8, 2.4, 10.0};
				}
			}
		}
		
		if(CentBinMin==4 && CentBinMax==6){ //10-60
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 1.0, 1.4, 1.8, 2.2, 2.8,5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.5, 1.0, 1.4, 2.0, 10.0};
				if(Energy>15){
					rebin_vector = {0.2, 0.5, 0.8, 1.2, 1.6, 2.2, 2.8, 10.0};
				}
				if(Energy==27){
					if(strncmp(PID, "Pion",4)==0){
						rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 5.0};
					}
					if(strncmp(PID, "Kaon",4)==0){
						rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.4, 1.8, 2.2, 2.6, 5.0};
					}
					if(strncmp(PID, "Proton",6)==0 ){
						rebin_vector = {0.2, 0.5, 0.8, 1.2, 1.8, 2.6, 10.0};
					}
				}
			}
		}
		if(strncmp(PID, "Hadrons",7)==0 ){
			if(Energy==7){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.4, 3.0, 10.0};
			}
			if(Energy==11){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.6, 2.2, 2.8, 10.0};
				if(CentBinMin==2 && CentBinMax==3){
					rebin_vector = {0.2, 0.4, 0.6, 1.0, 1.4, 2.0, 2.6, 10.0};
				}
			}
			if(Energy==14){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 10.0};
			}
			if(Energy==19){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4, 3.0, 10.0};
			}
			if(Energy==27 || Energy==39 || Energy==62){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8 , 10.0};
			}
		}
	}
	
	if(Energy==27 && harmonic==2){
		if(CentBinMin==0 && CentBinMax==8){ //0-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5,  10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.15, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 10.0};
			}
		}
		if(CentBinMin==7 && CentBinMax==8){ //0-10
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.2, 2.5, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 10.0};
			}
		}
		if( (CentBinMin==4 && CentBinMax==6) || (CentBinMin==4 && CentBinMax==8)){ //10-40
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 10.0};
			}
		}
		if(CentBinMin==0 && CentBinMax==3){ //40-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 2.5, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.5, 10.0};
	
			}	
		}
	}
	
	if(strncmp(PID, "Hadrons",7)==0 ){
		if(harmonic==2){
			rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 10.0};
			if( Energy < 20){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.2, 2.8, 10.0};
			}
		}
	}
	/*
		}else{
			//rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 10.0};
			rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 5.0};
			if( Energy < 20){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.6, 2.2, 2.8, 5.0};
			}
		}
	}
	*/
	return rebin_vector;
}
/*
void DrawEta(TFile *file, Int_t centMin, Int_t centMax,Int_t Energy){
	std::cout<<centMin<<"\t"<<centMax<<"\n";
	TCanvas *can = new TCanvas("canvas","plot",640,480);
	std::vector<TH1D *> histo;
	for(Int_t c=centMin; c<=centMax; c++){
		histo.push_back((TH1D*)file->Get(Form("hEtaCentBin%i",c)));
	}
	
	TH1D result;
	result=histo[0];

	for(Int_t c=1; c<histo.size(); c++){
		result->Add(histo[c],1.0);
	}

	Double_t max = result->GetBinContent(result->GetMaximumBin());
	TLegend *legend = new TLegend(0.2,0.82,0.7,0.89);
	legend->AddEntry(result,Form("Centrality %i-%i %%",cent[centMax+1], cent[centMin]),"p");
	result->SetAxisRange(0 ,1.2*max,"Y");
	result->SetLineColor(4);
	result->SetLineWidth(2);
	result->Draw();
	legend->Draw("same");

	histo.clear();


	can->SaveAs(Form("/home/demanov/New_work/MyAnalysis/picture/2_Charge_Hadrons/eta/Eta_%iGeV_cent_%i_%i.png",Energy, cent[centMax+1], cent[centMin]));
}
*/
void DrawEtaPID(TFile *file, Int_t centMin, Int_t centMax, Double_t ptMin, Double_t ptMax, Int_t Energy, TString PID, TString charge, TString EtaGap){
	TCanvas *can = new TCanvas("canvas","plot",640,480);
	std::vector<TH2D *> histo;
	for(Int_t c=centMin; c<=centMax; c++){
		histo.push_back(  (TH2D*)file->Get(Form("h_rapidityEastTPC%s%s%sEtaTPCandTOFCentBin%i",PID.Data(),charge.Data(),EtaGap.Data(),c)) );
		histo.push_back(  (TH2D*)file->Get(Form("h_rapidityWestTPC%s%s%sEtaTPCandTOFCentBin%i",PID.Data(),charge.Data(),EtaGap.Data(),c)) );
	}
	for(Int_t c=1; c<histo.size(); c++){
		histo[0]->Add(histo[c],1.0);
	}
	TH1D *result;

	Int_t nBinsPt=100;
	Double_t pt_min = 0.15;
	Double_t pt_max = 5.15;
	Double_t CenaDeleniy = (pt_max - pt_min)/100.;
	Int_t ptBinMin = (int)((ptMin-0.15) / CenaDeleniy) + 1;
	Int_t ptBinMax = (int)((ptMax-0.15) / CenaDeleniy);

	result=(TH1D*)histo[0] -> ProjectionX("",ptBinMin,ptBinMax);

	Double_t max = result->GetBinContent(result->GetMaximumBin());
	TLegend *legend = new TLegend(0.2,0.82,0.7,0.89);
	legend->AddEntry(result,Form("Centrality %i-%i %%",cent[centMax+1], cent[centMin]),"p");
	result->SetAxisRange(0 ,1.2*max,"Y");
	result->SetAxisRange(-1.1 ,1.1,"X");
	result->SetLineColor(4);
	result->SetLineWidth(2);
	result->Draw();
	legend->Draw("same");

	can->SaveAs(Form("/home/demanov/New_work/MyAnalysis/picture/2_Charge_Hadrons/eta/pid/Eta_%iGeV_%s%s_%s_cent_%i_%i.png",Energy, PID.Data(), charge.Data(), EtaGap.Data(), cent[centMax+1], cent[centMin]));
	delete result;
	histo.clear();
}


void PrintPoint(Int_t Energy, Int_t harmonic ,Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, TString PID, TString charge){

	//VEC rebin_vec = Rebin2(harmonic, Energy, PID, charge, CentBinMin, CentBinMax);
	//VEC rebin_vec = Rebin2(harmonic, Energy, PID, charge, CentBinMin, CentBinMax);
	VEC rebin_vec = Rebin2K(harmonic, Energy, PID, charge, CentBinMin, CentBinMax);
	TString prefix;
	
	if( strncmp(PID, "Hadrons",7)==0){
		prefix="";
    }else{
    	prefix="TPCandTOF";
    }
	
	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<Draw_Object *> flow;
		
	TFile *file_read;
	TFile *file_out = new TFile("OUT.root","RECREATE");	
	if( strncmp(PID, "Hadrons",7)==0){
		file_read = new TFile(Form("%s/file_hadrons/Sys/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2_sys.root",outPath,Energy),"READ");	
    }else{
		//file_read = new TFile(Form("%s/file/Flow_New%iGeV_PID_10binYesTofdca11_oldPID_BadRunOld.root",outPath, Energy),"READ");
		file_read = new TFile(Form("%s/file_pid/Flow_%iGeV_PID_10binYesTofdca11_sys.root",outPath, Energy),"READ");
    }

    //Data line
    flow.push_back(new Draw_Object);
    Analysis -> SetParametrs(Energy, file_read,file_out,harmonic,PID,charge,EtaGap,prefix);
    flow.back() -> SetParametrsGraph("data","#font[42]{Data}", 2, 33, 1.4);
    Analysis -> FlowVsPt_EtaSub(flow.back(), CentBinMin, CentBinMax, rebin_vec);

    TString text_vn;
    text_vn="std::vector<Double_t> ";
    text_vn+=Energy;
    text_vn+="GeV_";
    text_vn+=PID;
    text_vn+=charge;
    text_vn+="_cent_";
    text_vn+=cent[CentBinMax+1];
    text_vn+="_";
    text_vn+=cent[CentBinMin];
    text_vn+="_v";
    text_vn+=harmonic;
    text_vn+=" = { ";

    TString text_vn_err;
    text_vn_err="std::vector<Double_t> ";
    text_vn_err+=Energy;
    text_vn_err+="GeV_";
    text_vn_err+=PID;
    text_vn_err+=charge;
    text_vn_err+="_cent_";
    text_vn_err+=cent[CentBinMax+1];
    text_vn_err+="_";
    text_vn_err+=cent[CentBinMin];
    text_vn_err+="_v";
    text_vn_err+=harmonic;
    text_vn_err+="_err";
    text_vn_err+=" = { ";

    TString text_pt;
    text_pt="std::vector<Double_t> ";
    text_pt+=Energy;
    text_pt+="GeV_";
    text_pt+=PID;
    text_pt+=charge;
    text_pt+="_cent_";
    text_pt+=cent[CentBinMax+1];
    text_pt+="_";
    text_pt+=cent[CentBinMin];
    text_pt+="_v";
    text_pt+=harmonic;
    text_pt+="_pt";
    text_pt+=" = { ";

    std::cout<<text_pt;
    for(Int_t i=0; i < flow[0]->GetSizeVector(); i++){
    	std::cout<<flow[0]->GetPointXGraph(i);
    	if(i<flow[0]->GetSizeVector()-1){
    		std::cout<<", ";
    	}
    }
    std::cout<<"};\n";

        std::cout<<text_vn;
    for(Int_t i=0; i < flow[0]->GetSizeVector(); i++){
    	std::cout<<flow[0]->GetPointYGraph(i);
    	if(i<flow[0]->GetSizeVector()-1){
    		std::cout<<", ";
    	}
    }
    std::cout<<"};\n";

        std::cout<<text_vn_err;
    for(Int_t i=0; i < flow[0]->GetSizeVector(); i++){
    	std::cout<<flow[0]->GetPointYErrorGraph(i);
    	if(i<flow[0]->GetSizeVector()-1){
    		std::cout<<", ";
    	}
    }
    std::cout<<"};\n\n";

       //canvas setening
    /*
    Draw_Picture_new *canvas = new Draw_Picture_new;
    
    canvas->SetTLine(0., 0., 3.0, 0., 1, 2, 2, "Up");
    canvas->SetTLine(0., 0., 3.0, 0., 1, 2, 2, "Up");
    canvas->SetTLine(0., 0.9, 3.0, 0.9, 1, 2, 2, "Down");
    canvas->SetTLine(0., 1.1, 3.0, 1.1, 1, 2, 2, "Down");
    canvas->SetTLine(0., 1., 3.0, 1., flow[0] -> GetColor(), 2, 1, "Down");
	canvas->SetAxisToCanvsWithBottomPanel(0., 3.0, -0.01, 0.25, 0.8, 1.2);
	canvas->SetDrawObject(flow,"Up");
	canvas->SetDrawObject(flow,"Down");
	canvas->SetLegend(0.1,0.55,0.4,0.88,1);
	TCanvas *result = canvas->CanvasWithBottomPanel(Form("#font[42]{ #scale[1.2]{ #splitline{Au+Au #sqrt{s_{NN}} =  GeV}{%s, centrality %i-%i %%}}}", particle[GetNumberParticle(PID,charge)].Data(),cent[CentBinMax+1], cent[CentBinMin]),0.36*3.0,0.8*0.25,
													"#font[42]{p_{T} [GeV/c]}","v_{2}","#font[42]{#frac{Data}{Published}}","");
	result->SaveAs(Form("%s/picture/print/v%i_%s%s_%i_%i_%iGeV_%s.png",outPath, harmonic,PID.Data(), charge.Data(), cent[CentBinMax+1], cent[CentBinMin], Energy, EtaGap.Data()));
	delete result;
	*/
    //canvas setening
	file_read->Close();
	file_out->Close();
}

void funPrintPoint(Int_t Energy){

	TString pid[] = {"Pion", "Kaon","Proton"};
	TString pid_ch[] = {"Pos","Neg"};
	
/*
	for(Int_t h=0; h<2; h++){
		PrintPoint(Energy, h+2, 7, 8, "Eta15", "Hadrons", "");
		for(Int_t c=8; c>0; c--){
			PrintPoint(Energy, h+2, c, c, "Eta15", "Hadrons", "");
		}
	}
*/
	for(Int_t h=0; h<2; h++){
		for(Int_t p=1; p<2; p++){
			for(Int_t ch=0; ch<2;ch++){

				PrintPoint(Energy, h+2, 7, 8, "Eta01", pid[p], pid_ch[ch]);
				PrintPoint(Energy, h+2, 4, 6, "Eta01", pid[p], pid_ch[ch]);
				PrintPoint(Energy, h+2, 4, 8, "Eta01", pid[p], pid_ch[ch]);
				PrintPoint(Energy, h+2, 2, 8, "Eta01", pid[p], pid_ch[ch]);

			}
		}
	}

}

VEC Rebin2H(Int_t harmonic, Int_t Energy, TString PID, TString charge, Int_t CentBinMin, Int_t CentBinMax){

	VEC rebin_vector;


	/*
	if(Energy==27 && harmonic==2){
		if(CentBinMin==0 && CentBinMax==8){ //0-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.2, 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.15, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.2, 4.0, 10.0};
			}
		}
		if(CentBinMin==7 && CentBinMax==8){ //0-10
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.2, 2.6, 3.2, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.6, 10.0};
			}
		}
		if(CentBinMin==4 && CentBinMax==6){ //10-40
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.6, 3.2 , 4.0, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.2 , 3.6, 10.0};
			}
		}
		if(CentBinMin==0 && CentBinMax==3){ //40-80
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.2, 2.6, 3.2, 10.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.6, 10.0};
	
			}	
		}
	}
	*/
	if(harmonic==2){

			
		if(Energy>25){
			if(strncmp(PID, "Pion",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4, 3.2, 5.0};
			}
			if(strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.8, 2.2, 2.8, 5.0};
			}
		}

		if(Energy==19){
			if(strncmp(PID, "Pion",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.6, 2.0, 2.6, 3.2, 5.0};;
			}
			if(strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 1.0, 1.4, 1.8, 2.4, 3.0, 5.0};;
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 1.0, 1.4, 1.8, 2.2, 2.8, 10.0};
			}
		}

		if(Energy==14){
			if(strncmp(PID, "Pion",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.6, 5.0};
			}
			if(strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.6, 1.0, 1.4, 1.8, 2.4, 5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.4, 0.8, 1.2, 1.8, 2.4, 10.0};
			}
		}

		if(Energy==11 || Energy==7){
			if(strncmp(PID, "Pion",4)==0){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.6, 2.2, 2.8, 5.0};
			}
			if(strncmp(PID, "Kaon",4)==0){
				rebin_vector = {0.2, 0.6, 1.0, 1.4, 2.0, 2.8, 5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.4, 0.8, 1.4, 2.6, 10.0};
			}
		}

	}

	if(harmonic==3){
 	
		if(strncmp(PID, "Pion",4)==0  ){
			if(Energy==7){
				rebin_vector = {0.2, 0.4, 0.8, 1.2, 1.6, 2.2, 2.8, 10.0};
			}
			if(Energy==11){
					rebin_vector = {0.2, 0.4, 0.6, 1.0, 1.4, 2.0, 2.6, 10.0};

				//rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.6, 2.2, 2.8, 10.0};
				if(CentBinMin==2 && CentBinMax==3){
					rebin_vector = {0.2, 0.6, 1.0, 1.4, 2.0, 2.6, 10.0};
				}
			}
			if(Energy==14){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.6, 2.0, 2.6, 10.0};
			}
			if(Energy==19){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 10.0};
			}
			if(Energy==27 || Energy==39 || Energy==62){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.2, 2.8, 10.0};
			}
		}

		if(strncmp(PID, "Kaon",4)==0  ){
			if(Energy==7){
				rebin_vector = {0.2, 0.4, 0.8, 1.2, 1.6, 2.2, 2.8, 10.0};
			}
			if(Energy==11){
					rebin_vector = {0.2, 0.6, 1.0, 1.4, 2.0, 2.6, 10.0};

				//rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.6, 2.2, 2.8, 10.0};
				if(CentBinMin==2 && CentBinMax==3){
					rebin_vector = {0.2, 0.6, 1.2, 1.8, 2.4, 10.0};
				}
			}
			if(Energy==14){
				rebin_vector = {0.2, 0.4, 0.6, 1.0, 1.4, 2.0,  2.6, 10.0};
			}
			if(Energy==19){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.6, 10.0};
			}
			if(Energy==27 || Energy==39 || Energy==62){
				rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.6, 10.0};
			}
		}

		if(strncmp(PID, "Proton",6)==0 ){
			if(Energy==7){
				rebin_vector =  {0.4, 0.8, 1.2, 1.6, 2.2, 2.8, 10.0};
			}
			if(Energy==11){
					rebin_vector = { 0.4, 0.6, 1.0, 1.6, 2.2, 2.8, 10.0};

				//rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.6, 2.2, 2.8, 10.0};
				if(CentBinMin==2 && CentBinMax==3){
					rebin_vector = {0.4, 0.6, 1.0, 1.4, 2.0, 2.6, 10.0};
				}
			}
			if(Energy==14){
				rebin_vector = {0.4, 0.6, 1.0, 1.6, 2.2, 2.8, 10.0};
			}
			if(Energy==19){
				rebin_vector = {0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.4, 10.0};
			}
			if(Energy==27 || Energy==39 || Energy==62){
				rebin_vector = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.2, 2.8 , 10.0};
			}
		}
	}

	return rebin_vector;
}