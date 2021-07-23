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
void DrawEtaPID(TFile *file, Int_t centMin, Int_t centMax, Double_t ptMin, Double_t ptMax, Int_t Energy, TString PID, TString charge, TString EtaGap);
//потоки как функция энергии
void FlowVsEnergyBESupdate(Double_t ptMin, Double_t ptMax, TString EtaGap, TString PID, TString charge, Double_t min_y_axis, Double_t max_y_axis, Double_t max_x_axis,Double_t MinYFlow, Double_t MaxYratio);
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
//разница частиц и античастиц в зависимости от энергии
void ParAntVsEnergyMultyPad( Int_t CentBinMin, Int_t CentBinMax, Double_t ptMin, Double_t ptMax, TString EtaGap,  Double_t min_x_axis, Double_t max_x_axis, Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio, Double_t min_y_axis_dif, Double_t max_y_axis_dif);
void FlowVsPtForBESmultyPad( Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, TString charge, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio, Double_t min_y_axis_dif, Double_t max_y_axis_dif);
void HadronsFlowVsPtVsCentForBES( TString EtaGap, TString PID, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis);
void HadronsFlowVsPtVsCentForBES_ver2( TString EtaGap, TString PID, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis);
void HadronsFlowVsPtVsCentForBES_ver2Int( TString EtaGap, TString PID, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis);
void ParAntFlowVsCentForBES( Int_t ptMin, Int_t ptMax, TString EtaGap, TString PID, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio, Double_t min_y_axis_dif, Double_t max_y_axis_dif);
void ParAntFlowVsPtForBESupdate( Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, Double_t min_x_axis, Double_t max_x_axis);
void DrawFLowVsEnergyBES();

void FlowDifferentMethodEPandSP(Int_t Energy, Int_t harmonic, Int_t CentBinMin, Int_t CentBinMax, TString EtaGap,TString PID, TString charge, Double_t max_y_axis, Double_t max_x_axis, Double_t MinYratio, Double_t MaxYratio);

void DrawSravOldNew27(Int_t Energy, Int_t harmonic, TString EtaGap, TString PID, Double_t min_x_axis, Double_t max_x_axis, Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_dif, Double_t max_y_axis_dif, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio);
void DrawSravOldNew27Hadrons(Int_t Energy, Int_t harmonic,Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, Double_t min_x_axis, Double_t max_x_axis, Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio);
void DrawSravRunIdGroupHadrons(Int_t Energy, Int_t harmonic,Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, Double_t min_x_axis, Double_t max_x_axis, Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio);

int Analysis_BES(Int_t Energy){
	

	makeplotstyle();
	TFile *fstyle = new TFile("style.root");
	TStyle *tsty = (TStyle *)fstyle->Get("PlottingInStyle");
	tsty->cd();

	gROOT->ForceStyle();
	gStyle->SetOptStat(0);
	/*
	DrawSravRunIdGroupHadrons(27, 2, 7,8, "Eta15", -0.05, 3.26, -0.02, 0.24, 0.92, 1.08);
	DrawSravRunIdGroupHadrons(27, 2, 4,6, "Eta15", -0.05, 3.26, -0.02, 0.25, 0.92, 1.08);
	DrawSravRunIdGroupHadrons(27, 2, 0,3, "Eta15", -0.05, 3.26, -0.02, 0.32, 0.92, 1.08);
	DrawSravRunIdGroupHadrons(27, 3, 7,8, "Eta15", -0.05, 3.26, -0.02, 0.15, 0.92, 1.08);
	DrawSravRunIdGroupHadrons(27, 3, 2,6, "Eta15", -0.05, 3.26, -0.02, 0.15, 0.92, 1.08);
	DrawSravRunIdGroupHadrons(27, 3, 2,8, "Eta15", -0.05, 3.26, -0.02, 0.15, 0.92, 1.08);
	
	
	DrawSravOldNew27Hadrons(27, 2, 7,8, "Eta15", -0.05, 3.26, -0.02, 0.24, 0.92, 1.08);
	DrawSravOldNew27Hadrons(27, 2, 4,6, "Eta15", -0.05, 3.26, -0.02, 0.25, 0.92, 1.08);
	DrawSravOldNew27Hadrons(27, 2, 0,3, "Eta15", -0.05, 3.26, -0.02, 0.32, 0.92, 1.08);
	DrawSravOldNew27Hadrons(27, 3, 7,8, "Eta15", -0.05, 3.26, -0.02, 0.15, 0.92, 1.08);
	DrawSravOldNew27Hadrons(27, 3, 2,6, "Eta15", -0.05, 3.26, -0.02, 0.15, 0.92, 1.08);
	DrawSravOldNew27Hadrons(27, 3, 2,8, "Eta15", -0.05, 3.26, -0.02, 0.15, 0.92, 1.08);
	*/

	//FlowDifferentMethodEPandSP(Int_t Energy, Int_t harmonic, Int_t CentBinMin, Int_t CentBinMax, TString EtaGap,TString PID, TString charge, Double_t max_y_axis, Double_t max_x_axis, Double_t MinYratio, Double_t MaxYratio);
	/*
	FlowDifferentMethodEPandSP(Energy, 2, 7, 8, "Eta15", "Hadrons", "", 0.16, 2.97, 0.8, 1.2);
	FlowDifferentMethodEPandSP(Energy, 2, 4, 6, "Eta15", "Hadrons", "", 0.26, 2.97, 0.8, 1.2);
	FlowDifferentMethodEPandSP(Energy, 2, 0, 3, "Eta15", "Hadrons", "", 0.32, 2.97, 0.8, 1.2);

	FlowDifferentMethodEPandSP(Energy, 3, 7, 8, "Eta15", "Hadrons", "", 0.16, 2.97, 0.8, 1.2);
	FlowDifferentMethodEPandSP(Energy, 3, 4, 6, "Eta15", "Hadrons", "", 0.26, 2.97, 0.8, 1.2);
	FlowDifferentMethodEPandSP(Energy, 3, 0, 3, "Eta15", "Hadrons", "", 0.32, 2.97, 0.8, 1.2);
	*/

	HadronsFlowVsPtVsCentForBES("Eta15", "Hadrons",-0.05, 2.86, -0.02, 0.37);
	//HadronsFlowVsPtVsCentForBES_ver2("Eta15", "Hadrons",-0.05, 2.86, -0.02, 0.27);
	//HadronsFlowVsPtVsCentForBES_ver2Int("Eta15", "Hadrons",-0.05, 2.86, -0.3, 5.67);
	
	/*
	DrawSravOldNew27(27, 2, "Eta01", "Pion", -0.05, 3.26, -0.02, 0.24, -0.03, 0.034, 0.82, 1.18);
	DrawSravOldNew27(27, 2, "Eta01", "Kaon", -0.05, 3.26, -0.02, 0.24, -0.03, 0.034, 0.82, 1.18);
	DrawSravOldNew27(27, 2, "Eta01", "Proton", -0.05, 3.26, -0.02, 0.32, -0.03, 0.034, 0.97, 2.18);

	DrawSravOldNew27(27, 3, "Eta01", "Pion", -0.05, 2.96, -0.002, 0.095, -0.03, 0.034, 0.82, 1.18);
	DrawSravOldNew27(27, 3, "Eta01", "Kaon", -0.05, 2.96, -0.002, 0.134, -0.03, 0.034, 0.82, 1.18);
	DrawSravOldNew27(27, 3, "Eta01", "Proton", -0.05, 2.96, -0.002, 0.115, -0.03, 0.034, 0.97, 2.18);
	*/
    //ParAntVsEnergyMultyPad( Int_t CentBinMin, Int_t CentBinMax, Double_t ptMin, Double_t ptMax, TString EtaGap,  Double_t min_x_axis, Double_t max_x_axis, Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio, Double_t min_y_axis_dif, Double_t max_y_axis_dif);
	///ParAntVsEnergyMultyPad(2,8,0.2,2.0,"Eta01", 6.8, 68, -0.005, 0.035, 0.82, 1.18, -0.03, 0.034);
	/*
 //FlowVsPtForBESmultyPad( Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, TString charge, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio, Double_t min_y_axis_dif, Double_t max_y_axis_dif);
	FlowVsPtForBESmultyPad(2, 8, "Eta01", "Pos",-0.04, 2.96, -0.02, 0.22, 0.82, 1.18, -0.03, 0.034);
	FlowVsPtForBESmultyPad(2, 8, "Eta01", "Neg",-0.04, 2.96, -0.02, 0.22, 0.82, 1.18, -0.03, 0.034);
   
    ParAntFlowVsPtForBESupdate(2, 8, "Eta01",-0.04, 2.96);
    ParAntFlowVsPtForBESupdate(2, 6, "Eta01",-0.04, 2.96);
    //HadronsFlowVsPtVsCentForBES(2, 8, "Eta01", "Pion",-0.04, 2.96, -0.02, 0.22, 0.82, 1.18, -0.03, 0.034);
    
    HadronsFlowVsPtVsCentForBES(2, 8, "Eta01", "Pion",-0.05, 2.96, -0.02, 0.22, 0.82, 1.18, -0.03, 0.034);
    HadronsFlowVsPtVsCentForBES(2, 8, "Eta01", "Kaon",-0.05, 2.96, -0.02, 0.22, 0.72, 1.48, -0.03, 0.034);
    HadronsFlowVsPtVsCentForBES(2, 8, "Eta01", "Proton",-0.05, 2.96, -0.02, 0.22, 0.89, 2.07, -0.01, 0.084);

    HadronsFlowVsPtVsCentForBES(2, 6, "Eta01", "Pion",-0.05, 2.96, -0.02, 0.22, 0.82, 1.18, -0.03, 0.034);
    HadronsFlowVsPtVsCentForBES(2, 6, "Eta01", "Kaon",-0.05, 2.96, -0.02, 0.22, 0.72, 1.48, -0.03, 0.034);
    HadronsFlowVsPtVsCentForBES(2, 6, "Eta01", "Proton",-0.05, 2.96, -0.02, 0.22, 0.89, 2.07, -0.01, 0.084);
 	*/
    //FlowVsRapidityForBES(, Double_t ptMin, Double_t ptMax, TString EtaGap, TString PID, TString charge,  Double_t min_y_axis, Double_t max_y_axis){
 	//FlowVsRapidityForBES(2,0.2,2.0,"Eta15","Hadrons", "",0.015,0.065);
 	//FlowVsRapidityForBES(2,0.2,3.2,"Eta15","Hadrons", "",0.015,0.065);
 	//FlowVsRapidityForBES(3,0.2,2.0,"Eta15","Hadrons", "",0.001,0.023);
 	//FlowVsRapidityForBES(3,0.2,3.2,"Eta15","Hadrons", "",0.001,0.023);
	//DrawFlowVsRapidity(Energy);
	//DrawFlowDifferentEtaGap("Hadrons", Energy);
	//DrawFlowVsPtForBES("Hadrons");
	//FlowVsEnergyBESupdate( 0.2, 3.2, "Eta15","Hadrons", "", 0.0, 0.115, 75, 0.5, 1.15);
	//FlowVsEnergyBESupdate( 0.2, 2.0, "Eta15","Hadrons", "", 0.0, 0.115, 75, 0.5, 1.15);


	return 0;

}

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
	TString data[]={"Run10 27 GeV", "Run18 27 GeV"};


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

void DrawSravOldNew27Hadrons(Int_t Energy, Int_t harmonic,Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, Double_t min_x_axis, Double_t max_x_axis, Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio){

	TString data[]={"Run10 27 GeV", "Run18 27 GeV"};
	VEC	rebin_vec = Rebin2( harmonic, Energy, "Hadrons", "", CentBinMin, CentBinMax);
	
	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<Draw_Object *> flow;
	std::vector<Draw_Object *> ratio;

	std::vector<Draw_Object *> flowPer;
	std::vector<Draw_Object *> ratioPer;

	const Int_t color[]={2, 1, 4, 6, 4, 2, 46};
	const Int_t style[]={23, 22, 25, 23, 34, 47, 8, 29};
	const Double_t size[]={1.8,1.8,1.8,1.8,1.8,1.8,1.8};

	TFile *file_read_Run10 = new TFile(Form("%s/file_hadrons/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2.root",outPath,Energy),"READ");	
	TFile *file_read_Run18 = new TFile(Form("%s/file_hadrons/Sys/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2_sys.root",outPath,Energy),"READ");	
	TFile *file_read_Run18per1 = new TFile(Form("%s/file_hadrons/Sys/Flow_%iGeVper1_Hadrons_2binYesTof_dcaEP2_sys.root",outPath,Energy),"READ");	
	TFile *file_read_Run18per2 = new TFile(Form("%s/file_hadrons/Sys/Flow_%iGeVper2_Hadrons_2binYesTof_dcaEP2_sys.root",outPath,Energy),"READ");	


	TFile *file_out = new TFile("OUT.root","RECREATE");

	flow.push_back(new Draw_Object);
    Analysis -> SetParametrs(Energy, file_read_Run10,file_out,harmonic,"Hadrons","","Eta15","");
    flow.back() -> SetParametrsGraph(Form("dataRun10"),Form("#font[42]{Run10 27GeV}"), color[0], style[0], size[0]);
    Analysis -> FlowVsPt_EtaSub(flow.back(), CentBinMin, CentBinMax, rebin_vec);
    //Data line
	flow.push_back(new Draw_Object);
    Analysis -> SetParametrs(Energy, file_read_Run18,file_out,harmonic,"Hadrons","","Eta15","");
    flow.back() -> SetParametrsGraph(Form("dataRun18"),Form("#font[42]{Run18 27GeV}"), color[1], style[1], size[1]);
    Analysis -> FlowVsPt_EtaSub(flow.back(), CentBinMin, CentBinMax, rebin_vec);
    //Data line

    ratio.push_back(new Draw_Object);
    Analysis -> RatioGraphEVAL(ratio.back(),flow[1],flow[0],"");
    //canvas setening
    Draw_Picture_new *canvas = new Draw_Picture_new;
    
    canvas->SetTLine(0., 0., max_x_axis, 0., 1, 2, 2, "Up");
    canvas->SetTLine(0., 0., max_x_axis, 0., 1, 2, 2, "Up");
    canvas->SetTLine(0., 0.95, max_x_axis, 0.95, 1, 2, 2, "Down");
    canvas->SetTLine(0., 1.05, max_x_axis, 1.05, 1, 2, 2, "Down");
    canvas->SetTLine(0., 1., max_x_axis, 1., flow[0] -> GetColor(), 2, 1, "Down");
	canvas->SetAxisToCanvsWithBottomPanel(0., max_x_axis, -0.01, max_y_axis, min_y_axis_ratio, max_y_axis_ratio);
	canvas->SetDrawObject(flow,"Up");
	canvas->SetDrawObject(ratio,"Down");
	canvas->SetLegend(0.1,0.55,0.4,0.88,1);
	TCanvas *result = canvas->CanvasWithBottomPanel(Form("#font[42]{ #scale[1.0]{ #splitline{Au+Au #sqrt{s_{NN}} = %s GeV}{%s, centrality %i-%i %%}}}",Energy_scan_text.at(Energy), particle[GetNumberParticle("Hadrons","")].Data(),cent[CentBinMax+1], cent[CentBinMin]),0.36*max_x_axis,0.8*max_y_axis,
													"#font[42]{p_{T} [GeV/c]}",Form("#font[42]{v_{%i}}",harmonic),Form("#font[42]{#frac{v_{%i}(Run18)}{v_{%i}(Run10)}}", harmonic, harmonic),"");
	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/27New/SravRun10Run18_%iGeV_v%i_pt_Hadrons_%i_%i.png",outPath, Energy,harmonic, cent[CentBinMax+1], cent[CentBinMin]));
	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/27New/SravRun10Run18_%iGeV_v%i_pt_Hadrons_%i_%i.pdf",outPath, Energy,harmonic, cent[CentBinMax+1], cent[CentBinMin]));
	delete result;
	
	///////// разные периоды 

	flowPer.push_back(new Draw_Object);
    Analysis -> SetParametrs(Energy, file_read_Run18,file_out,harmonic,"Hadrons","","Eta15","");
    flowPer.back() -> SetParametrsGraph(Form("dataRun10"),Form("#font[42]{All}"), color[0], style[0], size[0]);
    Analysis -> FlowVsPt_EtaSub(flowPer.back(), CentBinMin, CentBinMax, rebin_vec);
    //Data line
	flowPer.push_back(new Draw_Object);
    Analysis -> SetParametrs(Energy, file_read_Run18per1,file_out,harmonic,"Hadrons","","Eta15","");
    flowPer.back() -> SetParametrsGraph(Form("dataRun18"),Form("#font[42]{period 1}"), color[1], style[1], size[1]);
    Analysis -> FlowVsPt_EtaSub(flowPer.back(), CentBinMin, CentBinMax, rebin_vec);

    //Data line
	flowPer.push_back(new Draw_Object);
    Analysis -> SetParametrs(Energy, file_read_Run18per2,file_out,harmonic,"Hadrons","","Eta15","");
    flowPer.back() -> SetParametrsGraph(Form("dataRun18"),Form("#font[42]{period 2}"), color[2], style[2], size[2]);
    Analysis -> FlowVsPt_EtaSub(flowPer.back(), CentBinMin, CentBinMax, rebin_vec);

    ratioPer.push_back(new Draw_Object);
    Analysis -> RatioGraphEVAL(ratioPer.back(),flowPer[1],flowPer[0],"");
    ratioPer.push_back(new Draw_Object);
    Analysis -> RatioGraphEVAL(ratioPer.back(),flowPer[2],flowPer[0],"");
    //canvas setening
    //canvas setening


	Draw_Picture_new *canvas2 = new Draw_Picture_new;
    
    canvas2->SetTLine(0., 0., max_x_axis, 0., 1, 2, 2, "Up");
    canvas2->SetTLine(0., 0., max_x_axis, 0., 1, 2, 2, "Up");
    canvas2->SetTLine(0., 0.98, max_x_axis, 0.98, 1, 2, 2, "Down");
    canvas2->SetTLine(0., 1.02, max_x_axis, 1.02, 1, 2, 2, "Down");
    canvas2->SetTLine(0., 1., max_x_axis, 1., flow[0] -> GetColor(), 2, 1, "Down");
	canvas2->SetAxisToCanvsWithBottomPanel(0., max_x_axis, -0.01, max_y_axis+0.001, 0.965, 1.035);
	canvas2->SetDrawObject(flowPer,"Up");
	canvas2->SetDrawObject(ratioPer,"Down");
	canvas2->SetLegend(0.1,0.55,0.4,0.88,1);
	TCanvas *result2 = canvas2->CanvasWithBottomPanel(Form("#font[42]{ #scale[1.0]{ #splitline{Au+Au #sqrt{s_{NN}} = %s GeV,Run18}{%s, centrality %i-%i %%}}}",Energy_scan_text.at(Energy), particle[GetNumberParticle("Hadrons","")].Data(),cent[CentBinMax+1], cent[CentBinMin]),0.36*max_x_axis,0.8*max_y_axis,
													"#font[42]{p_{T} [GeV/c]}",Form("#font[42]{v_{%i}}",harmonic),Form("#font[42]{#frac{v_{%i}(Per1,Per2)}{v_{%i}(All)}}", harmonic, harmonic),"");
	result2->SaveAs(Form("%s/picture/2_Charge_Hadrons/27New/SravPer1Per2_%iGeV_v%i_pt_Hadrons_%i_%i.png",outPath, Energy,harmonic, cent[CentBinMax+1], cent[CentBinMin]));
	result2->SaveAs(Form("%s/picture/2_Charge_Hadrons/27New/SravPer1Per2_%iGeV_v%i_pt_Hadrons_%i_%i.pdf",outPath, Energy,harmonic, cent[CentBinMax+1], cent[CentBinMin]));
	delete result2;

	file_read_Run10->Close();
	file_read_Run18->Close();
	file_read_Run18per1->Close();
	file_read_Run18per2->Close();
	file_out->Close();


}

void DrawSravRunIdGroupHadrons(Int_t Energy, Int_t harmonic,Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, Double_t min_x_axis, Double_t max_x_axis, Double_t min_y_axis, Double_t max_y_axis, Double_t min_y_axis_ratio, Double_t max_y_axis_ratio){

	TString data[]={"Run10 27 GeV", "Run18 27 GeV"};
	VEC	rebin_vec = Rebin2( harmonic, Energy, "Hadrons", "", CentBinMin, CentBinMax);
	
	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<Draw_Object *> flow;
	std::vector<Draw_Object *> ratio;

	std::vector<Draw_Object *> flowPer;
	std::vector<Draw_Object *> ratioPer;

	const Int_t color[]={2, 1, 4, 6, 4, 2, 46};
	const Int_t style[]={23, 22, 25, 23, 34, 47, 8, 29};
	const Double_t size[]={1.8,1.8,1.8,1.8,1.8,1.8,1.8};

	TString run[]={"Run gr1","Run gr2","Run gr3","Run gr4","Run gr5"};

	TFile *file_read_Run18 = new TFile(Form("%s/file_hadrons/Sys/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2_sys.root",outPath,Energy),"READ");	


	TFile *file_read[6];
	for(Int_t i=0; i<5; i++){
		file_read[i] = new TFile(Form("%s/file_hadrons/Sys/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2_sys_%i.root",outPath,Energy,i+1),"READ");	
	}	

	TFile *file_out = new TFile("OUT.root","RECREATE");

	flow.push_back(new Draw_Object);
    Analysis -> SetParametrs(Energy, file_read_Run18,file_out,harmonic,"Hadrons","","Eta15","");
    flow.back() -> SetParametrsGraph(Form("dataRun10"),Form("#font[42]{Run18 all}"), color[0], style[0], size[0]);
    Analysis -> FlowVsPt_EtaSub(flow.back(), CentBinMin, CentBinMax, rebin_vec);
    //Data line

	for(Int_t i=0; i<5; i++){
		flow.push_back(new Draw_Object);
    	Analysis -> SetParametrs(Energy, file_read[i],file_out,harmonic,"Hadrons","","Eta15","");
    	flow.back() -> SetParametrsGraph(Form("dataRun10"),Form("#font[42]{%s}", run[i].Data()), color[i+1], style[i+1], size[i+1]);
    	Analysis -> FlowVsPt_EtaSub(flow.back(), CentBinMin, CentBinMax, rebin_vec);

    	ratio.push_back(new Draw_Object);
   		Analysis -> RatioGraphEVAL(ratio.back(),flow[i+1],flow[0],"");

	}
  
    //canvas setening
    Draw_Picture_new *canvas = new Draw_Picture_new;
    
    canvas->SetTLine(0., 0., max_x_axis, 0., 1, 2, 2, "Up");
    canvas->SetTLine(0., 0., max_x_axis, 0., 1, 2, 2, "Up");
    canvas->SetTLine(0., 0.95, max_x_axis, 0.95, 1, 2, 2, "Down");
    canvas->SetTLine(0., 1.05, max_x_axis, 1.05, 1, 2, 2, "Down");
    canvas->SetTLine(0., 1., max_x_axis, 1., flow[0] -> GetColor(), 2, 1, "Down");
	canvas->SetAxisToCanvsWithBottomPanel(0., max_x_axis, -0.01, max_y_axis, min_y_axis_ratio, max_y_axis_ratio);
	canvas->SetDrawObject(flow,"Up");
	canvas->SetDrawObject(ratio,"Down");
	canvas->SetLegend(0.1,0.45,0.4,0.88,1);
	TCanvas *result = canvas->CanvasWithBottomPanel(Form("#font[42]{ #scale[1.0]{ #splitline{Au+Au #sqrt{s_{NN}} = %s GeV}{%s, centrality %i-%i %%}}}",Energy_scan_text.at(Energy), particle[GetNumberParticle("Hadrons","")].Data(),cent[CentBinMax+1], cent[CentBinMin]),0.36*max_x_axis,0.8*max_y_axis,
													"#font[42]{p_{T} [GeV/c]}",Form("#font[42]{v_{%i}}",harmonic),Form("#font[42]{#frac{v_{%i}(Run gr)}{v_{%i}(Run All)}}", harmonic, harmonic),"");
	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/27New/RunSrav_%iGeV_v%i_pt_Hadrons_%i_%i.png",outPath, Energy,harmonic, cent[CentBinMax+1], cent[CentBinMin]));
	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/27New/RunSrav_%iGeV_v%i_pt_Hadrons_%i_%i.pdf",outPath, Energy,harmonic, cent[CentBinMax+1], cent[CentBinMin]));
	delete result;
	
	for(Int_t i=0; i<5; i++){
		file_read[i]->Close();
	}
	file_read_Run18->Close();
	file_out->Close();


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
    		flow.back()->СhangePointGraph(p,(-1)*flow.back()->GetPointXGraph(p), flow.back()->GetPointYGraph(p), 0, flow.back()->GetPointYErrorGraph(p), flow.back()->GetPointYSysErrorGraph(p));
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


void FlowVsEnergyBESupdate(Double_t ptMin, Double_t ptMax, TString EtaGap, TString PID, TString charge, Double_t min_y_axis, Double_t max_y_axis, Double_t max_x_axis,Double_t MinYratio, Double_t MaxYratio){

	VEC rebin_vec = {0.15,ptMin,ptMax,5.2};;
	Int_t point = 1; 
	TString prefix="";
	
	Int_t Arr_energy[6]={11,14,19,27,39,62};
	Double_t x_energy[6]={11.5,14.5,19.6,27,39,62.4};
	Int_t cent_bin_Min[8]={7,6,5,4,3,2,1,0};
	Int_t cent_bin_Max[8]={8,6,5,4,3,2,1,0};

	TFile *file_1[6];
	for(Int_t i=0; i<6; i++){
		if( strncmp(PID, "Hadrons",7)==0){
			file_1[i] = new TFile(Form("%s/file_hadrons/Sys/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2_sys.root",outPath,Arr_energy[i]),"READ");	
	    }else{
			file_1[i] = new TFile(Form("%s/file/Flow_%iGeV_PID_10binYesTofdca11_oldPID_BadRunOld.root",outPath, Arr_energy[i]),"READ");
	    }
	}
	TFile *file_out = new TFile("OUT.root","RECREATE");	

	Analysis_Function *Analysis = new Analysis_Function;
	Draw_Object *object = new Draw_Object;
	std::vector<std::vector<Draw_Object *>> flow ;
	std::vector<std::vector<Draw_Object *>> ratio ;
	Draw_Picture_new *canvasFlow = new Draw_Picture_new;

	Double_t x,y,ex,ey;
	Double_t den,den_err;
	Int_t b=0;

	for(Int_t c=0; c<6; c++){
		flow.push_back({new Draw_Object, new Draw_Object});
		ratio.push_back({new Draw_Object, new Draw_Object});

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
		}
		
		for(Int_t harm=0; harm<2; harm++){
			flow[c][harm]->GetPointDrawObject(6,x,den,ex,den_err);
			Analysis->RatioGraphToOnePoint(ratio[c][harm],den, den_err);			
		}

		flow[c][0] -> SetParametrsGraph("data","v_{2}",4,23,1.4);
		ratio[c][0] -> SetParametrsGraph("data","v_{2}",4,23,1.4);
		
		flow[c][1] -> SetParametrsGraph("data","v_{3}",2,22,1.4);
		ratio[c][1] -> SetParametrsGraph("data","v_{3}",2,22,1.4);

		canvasFlow->SetNumberPad(Form("#font[42]{%s %i-%i %%}", bukva[b].Data(),cent[cent_bin_Max[c]+1], cent[cent_bin_Min[c]]), "Pad", 0.04, 0.85*max_y_axis, 0.09, 0.0);
		b++;

	}
	//
	std::cout<<"good\n";

   	canvasFlow->SetTextInfo(Form("#scale[1.7]{v_{n}}"),"Y",0.5, 0.7,0.3,90);
	canvasFlow->SetTextInfo("#font[42]{#scale[0.65]{#sqrt{s_{NN}[GeV]}}}","X",0.75,0.6,0.5,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.5]{ #splitline{Au+Au, %s}{ %.1f<p_{T}<%.1f GeV/c }}}", particle[GetNumberParticle(PID,charge)].Data(),ptMin, ptMax),"Info",0.3,0.85*max_y_axis,0.09,0.0);
   	canvasFlow->SetTLine(0.03, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasFlow->SetAxisToCanvsNM(0.03, max_x_axis,  min_y_axis, max_y_axis, 0.086);
	canvasFlow->SetAxisToCanvsNM(0.03, max_x_axis,  min_y_axis, max_y_axis, 0.086);
	canvasFlow->SetDrawObject(flow[0],"OnePad");
	canvasFlow->SetLegend(0.1,0.5,0.93,0.7,1);
	TCanvas *result = canvasFlow->CanvasNxM(1200,800, 3, 2, 0.35, 0.25, flow, 2,0,1);

	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/vn_pt_onePad_%.1f_%.1f_%s%s_BES_%s.pdf",outPath, ptMin, ptMax,PID.Data(), charge.Data(), EtaGap.Data()) );
	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/vn_pt_onePad_%.1f_%.1f_%s%s_BES_%s.png",outPath, ptMin, ptMax,PID.Data(), charge.Data(), EtaGap.Data()) );
	delete result;






/*
   	canvas->SetTLine(0., 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvas->SetAxisToCanvsNM(0., max_x_axis,  min_y_axis, max_y_axis, 0.07);
	canvas->SetDrawObject(flow[0],"OnePad");
	canvas->SetLegend(0.1,0.5,0.93,0.7,1);
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
	canvas2->SetLegend(0.1,0.5,0.93,0.7,1);
	TCanvas *result2 = canvas2->CanvasNxM(3, 3, 0.35, 0.35, ratio, Form("#font[42]{ #scale[1.5]{ #splitline{Au+Au, %s}{ %.1f<p_{T}<%.1f GeV/c }}}", particle[GetNumberParticle(PID,charge)].Data(),ptMin, ptMax), 5,0.92*max_ratio_y,
											TitlePad, 60, 0.86*max_ratio_y,
											"#font[42]{#scale[0.65]{#sqrt{s_{NN}[GeV]}}}",Form("#font[42]{#scale[1.0]{#frac{ v_{n}}{ v_{n}(62.4 GeV)}}}"),"",1);

	result2->SaveAs(Form("%s/picture/2_Charge_Hadrons/vn_pt_onePad_ratio_%.1f_%.1f_%s%s_BES_%s.png",outPath, ptMin, ptMax,PID.Data(), charge.Data(), EtaGap.Data()));
	
*/
	//delete result;
	
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
	
	std::vector<TString> NameEtaHadrons = {"Eta15","Eta03","Eta05","Eta07"};
	std::vector<Double_t> EtaVecHadrons = {0.075, 0.15, 0.25, 0.35};
	std::vector<TString> NameEtaPID = {"Eta1","Eta04","Eta07","Eta10"};
	std::vector<Double_t> EtaVecPID = {0.05, 0.2, 0.35, 0.5};
	TString suf = "def";

	const Int_t color[]={2, 1, 4, 6, 4, 2, 46};
	const Int_t style[]={23, 22, 25, 23, 34, 47, 8, 29};
	const Double_t size[]={2.0,2.,2.,2.,2.,2.,2.};

	TFile *file_read;

	if( strncmp(PID, "Hadrons",7)==0){
		file_read = new TFile(Form("%s/file_hadrons/Sys/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2_sys.root",outPath,Energy),"READ");	
    }else{
		file_read = new TFile(Form("%s/file/Flow_%iGeV_PID_10binYesTofdca11_oldPID_BadRunOld.root",outPath, Energy),"READ");
    }

	TFile *file_out = new TFile("OUT.root","RECREATE");

	for(Int_t eta=0; eta<4; eta++){
	    flow.push_back(new Draw_Object);
    	if(strncmp(PID, "Hadrons",7)==0){
    		Analysis -> SetParametrs(Energy, file_read,file_out,harmonic,PID,charge,NameEtaHadrons[eta],prefix);
    		flow.back() -> SetParametrsGraph(Form("data%.3f",EtaVecHadrons[eta]),Form("#font[42]{#eta-gap=%.2f %s}",2.0*EtaVecHadrons[eta], suf.Data()), color[eta], style[eta], size[eta]);
    	}else{
    		Analysis -> SetParametrs(Energy, file_read,file_out,harmonic,PID,charge,NameEtaPID[eta],prefix);
    		flow.back() -> SetParametrsGraph(Form("data%.3f",EtaVecPID[eta]),Form("#font[42]{#eta-gap=%.2f %s}",2.0*EtaVecPID[eta], suf.Data()), color[eta], style[eta], size[eta]);
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
	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/different_eta_gap/DifferentEtaGap_%iGeV_v%i_pt_%s%s_%i_%i.png",outPath, Energy,harmonic,PID.Data(), charge.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/different_eta_gap/DifferentEtaGap_%iGeV_v%i_pt_%s%s_%i_%i.pdf",outPath, Energy,harmonic,PID.Data(), charge.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result;
	file_read->Close();
	file_out->Close();

}
/*
void HadronsFlowVsPtVsCentForBES( TString EtaGap, TString PID, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis){

	VEC rebin_vec;
	TString prefix = "";
	
	Int_t Arr_energy[6]={11,14,19,27,39,62};
	const Int_t color[]={6, 1, 4, 2, 1, 2, 46};
	const Int_t style[]={8, 22, 23, 28, 25, 34, 8, 29};
	const Double_t size[]={1.3,1.4,1.4,1.4,1.4,1.4,1.4};
	
	Int_t centBinMin[6] = {7,6,5,4,2};
	Int_t centBinMax[6] = {8,6,5,4,3};

	Float_t k[2]={1.0,2.0};

	TFile *file_1[6];
	for(Int_t i=0; i<6; i++){
		file_1[i] = new TFile(Form("%s/file_hadrons/Sys/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2_sys.root",outPath, Arr_energy[i]),"READ");
	}
	TFile *file_out = new TFile("OUT.root","RECREATE");

	Draw_Picture_new *canvasFlow = new Draw_Picture_new;
	Draw_Picture_new *canvasDifference = new Draw_Picture_new;
	Draw_Picture_new *canvasRatio = new Draw_Picture_new;
	
	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<std::vector<Draw_Object *>> flow;
	std::vector<Draw_Object *> object;
	
	Int_t b=0;

	for(Int_t j=0; j<6; j++){
		
		canvasFlow->SetTextInfo(Form("#font[42]{%s GeV}", Energy_scan_text.at(Arr_energy[j])), "Pad", 0.35, 0.85*max_y_axis, 0.13, 0.0);
		
		for(Int_t harm=0; harm<2; harm++){
			for(Int_t c=0; c<5; c++){

				rebin_vec = Rebin2( harm+2, Arr_energy[j], PID, "", centBinMin[c], centBinMax[c]);
				
				object.push_back(new Draw_Object);
				Analysis -> SetParametrs(Arr_energy[j], file_1[j], file_out, harm+2, PID, "", EtaGap, prefix);
				Analysis -> FlowVsPt_EtaSub(object.back(), centBinMin[c], centBinMax[c], rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i",harm+2),Form("#font[42]{#scale[1.5]{%i-%i%%}}",cent[centBinMax[c]+1], cent[centBinMin[c]] ), color[c], style[c], size[c]);

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

			}

			flow.push_back(object);
			object.clear();

			if(j==5){
				if(harm==0) canvasFlow->SetTextInfo(Form("#font[42]{(v_{%i})}",harm+2), Form("V%i",harm+2), 0.06, 0.7*max_y_axis, 0.13, 0.0);
				if(harm==1) canvasFlow->SetTextInfo(Form("#font[42]{(%.1f #times v_{%i})}", k[harm],harm+2), Form("V%i",harm+2), 0.06, 0.7*max_y_axis, 0.13, 0.0);
			}

			canvasFlow->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.85*max_y_axis, 0.13, 0.0);
			b++;
		}
	}
	canvasFlow->SetTextInfo("#font[42]{v_{n}}","Y",0.25,0.6,0.5,0.0);
	canvasFlow->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %s}}", particle[GetNumberParticle(PID,"")].Data()),"Info",1.,0.85*max_y_axis,0.13,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{STAR Preliminary}}"),"Prel",0.4,0.85*max_y_axis,0.13,0.0);
	canvasFlow->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasFlow->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.13);
	canvasFlow->SetDrawObject(flow[0],"OnePad");
	canvasFlow->SetLegend(0.1,0.5,0.58,0.98,2);
	TCanvas *result = canvasFlow->CanvasNxM(1000,1280, 6, 2, 0.4, 0.68, flow,1,0,3);

	result->SaveAs(Form("%s/picture/3_PID/Sys/SysBES_vn_%s_%s_pt_cent_0_60.pdf",outPath, PID.Data(), EtaGap.Data()));
	result->SaveAs(Form("%s/picture/3_PID/Sys/SysBES_vn_%s_%s_pt_cent_0_60.png",outPath, PID.Data(), EtaGap.Data()));
	result->SaveAs(Form("%s/picture/3_PID/Sys/SysBES_vn_%s_%s_pt_cent_0_60.C",outPath, PID.Data(), EtaGap.Data()));
	delete result;

	for(Int_t i=0; i<6; i++){
		file_1[i]->Close();
		delete file_1[i];
	}

}
*/

/////картинка для пкс
void HadronsFlowVsPtVsCentForBES( TString EtaGap, TString PID, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis){

	VEC rebin_vec;
	TString prefix = "";
	
	//Int_t Arr_energy[6]={11,14,19,27,39,62};
	Int_t Arr_energy[6]={19,27,39,62};
	const Int_t color[]={6, 1, 4, 2, 1, 2, 46};
	const Int_t style[]={8, 22, 23, 28, 25, 34, 8, 29};
	const Double_t size[]={1.3,1.4,1.4,1.4,1.4,1.4,1.4};
	
	Int_t centBinMin[6] = {7,6,5,4,2};
	Int_t centBinMax[6] = {8,6,5,4,3};

	Float_t k[2]={1.0,2.0};

	max_y_axis=0.25;

	TFile *file_1[6];
	for(Int_t i=0; i<1; i++){
		file_1[i] = new TFile(Form("%s/file_hadrons/Sys/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2_sys.root",outPath, Arr_energy[i]),"READ");
	}
	TFile *file_out = new TFile("OUT.root","RECREATE");

	Draw_Picture_new *canvasFlow = new Draw_Picture_new;
	Draw_Picture_new *canvasDifference = new Draw_Picture_new;
	Draw_Picture_new *canvasRatio = new Draw_Picture_new;
	
	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<std::vector<Draw_Object *>> flow;
	std::vector<Draw_Object *> object;
	
	Int_t b=0;

	//for(Int_t j=0; j<6; j++){
	for(Int_t j=0; j<1; j++){
		
		//canvasFlow->SetTextInfo(Form("#font[42]{%s GeV}", Energy_scan_text.at(Arr_energy[j])), "Pad", 0.35, 0.85*max_y_axis, 0.07, 0.0);
		
		for(Int_t harm=0; harm<2; harm++){
			for(Int_t c=0; c<5; c++){

				rebin_vec = Rebin2( harm+2, Arr_energy[j], PID, "", centBinMin[c], centBinMax[c]);
				
				object.push_back(new Draw_Object);
				Analysis -> SetParametrs(Arr_energy[j], file_1[j], file_out, harm+2, PID, "", EtaGap, prefix);
				Analysis -> FlowVsPt_EtaSub(object.back(), centBinMin[c], centBinMax[c], rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i",harm+2),Form("#font[42]{#scale[1.5]{%i-%i%%}}",cent[centBinMax[c]+1], cent[centBinMin[c]] ), color[c], style[c], size[c]);

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

			}

			flow.push_back(object);
			object.clear();

			//if(j==5){
			if(j==0){
				if(harm==0) canvasFlow->SetTextInfo(Form("#font[42]{(v_{%i})}",harm+2), Form("V%i",harm+2), 0.06, 0.8*max_y_axis, 0.07, 0.0);
				if(harm==1) canvasFlow->SetTextInfo(Form("#font[42]{(%.1f #times v_{%i})}", k[harm],harm+2), Form("V%i",harm+2), 0.06, 0.8*max_y_axis, 0.07, 0.0);
			}

			canvasFlow->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.9*max_y_axis, 0.07, 0.0);
			b++;
		}
	}
	canvasFlow->SetTextInfo("#font[42]{v_{n}}","Y",0.25,0.6,0.5,0.0);
	canvasFlow->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.4,0.7,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{#splitline{Au+Au, #sqrt{s_{NN}}=19.6 GeV}{%s}}}", particle[GetNumberParticle(PID,"")].Data()),"Info",0.35,0.85*max_y_axis,0.07,0.0);
	//canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{STAR Preliminary}}"),"Prel",0.4,0.85*max_y_axis,0.07,0.0);
	canvasFlow->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasFlow->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.07);
	canvasFlow->SetDrawObject(flow[0],"OnePad");
	canvasFlow->SetLegend(0.4,0.7,0.9,0.98,2);
	//TCanvas *result = canvasFlow->CanvasNxM(1000,1280, 6, 2, 0.4, 0.68, flow,1,0,3);
	TCanvas *result = canvasFlow->CanvasNxM(1280,640, 1, 2, 0.4, 0.2, flow,1,0,1);

	result->SaveAs(Form("%s/picture/3_PID/Sys/SysBES_vn_%s_%s_pt_cent_pks19n.pdf",outPath, PID.Data(), EtaGap.Data()));
	result->SaveAs(Form("%s/picture/3_PID/Sys/SysBES_vn_%s_%s_pt_cent_pks19n.png",outPath, PID.Data(), EtaGap.Data()));
	delete result;

	for(Int_t i=0; i<1; i++){
		file_1[i]->Close();
		delete file_1[i];
	}

}
////////////////////////////



void HadronsFlowVsPtVsCentForBES_ver2( TString EtaGap, TString PID, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis){

	VEC rebin_vec;
	VEC rebin_vec_integral={0.2,3.2,5.2};
	TString prefix = "";
	
	Int_t Arr_energy[6]={11,14,19,27,39,62};
	const Int_t color[]={6, 1, 4, 2, 1, 2, 46};
	const Int_t style[]={8, 22, 23, 24, 25, 34, 8, 29};
	const Double_t size[]={1.3,1.4,1.4,1.4,1.4,1.4,1.4};
	
	Int_t centBinMin[6] = {7,6,5,4,2};
	Int_t centBinMax[6] = {8,6,5,4,3};

	Float_t k[2]={1.0,2.0};

	TFile *file_1[6];
	for(Int_t i=0; i<6; i++){
		file_1[i] = new TFile(Form("%s/file_hadrons/Sys/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2_sys.root",outPath, Arr_energy[i]),"READ");
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

				rebin_vec = Rebin2( harm+2, Arr_energy[j], PID, "", centBinMin[c], centBinMax[c]);
				
				object.push_back(new Draw_Object);
				Analysis -> SetParametrs(Arr_energy[j], file_1[j], file_out, harm+2, PID, "", EtaGap, prefix);
				Analysis -> FlowVsPt_EtaSub(object.back(), centBinMin[c], centBinMax[c], rebin_vec);
	    		//Analysis -> FlowVsPtDivideIntegralFLow_EtaSub(object.back(), centBinMin[c], centBinMax[c], rebin_vec, rebin_vec_integral, 0);
	    		object.back() -> SetParametrsGraph(Form("data_v%i_%iGeV_Hadrons",harm+2,Arr_energy[j]),Form("#font[42]{#scale[1.3]{%s GeV}}", Energy_scan_text.at(Arr_energy[j])), color[j], style[j], size[j]);

				if(harm==1){
					for(Int_t n=0; n < object.back()->GetSizeVector(); n++){
						object.back()->СhangePointGraph(n,
														object.back()->GetPointXGraph(n), 
														object.back()->GetPointYGraph(n) * k[1],
														object.back()->GetPointXErrorGraph(n),
														object.back()->GetPointYErrorGraph(n),
														object.back()->GetPointYSysErrorGraph(n) * k[1]);
					}
				}

			}

			flow.push_back(object);
			object.clear();

			if(c==4){
				if(harm==0) canvasFlow->SetTextInfo(Form("#font[42]{(v_{%i})}",harm+2), Form("V%i",harm+2), 0.06, 0.7*max_y_axis, 0.13, 0.0);
				if(harm==1) canvasFlow->SetTextInfo(Form("#font[42]{(%.1f #times v_{%i})}", k[harm],harm+2), Form("V%i",harm+2), 0.06, 0.7*max_y_axis, 0.13, 0.0);
			}

			canvasFlow->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.85*max_y_axis, 0.13, 0.0);
			b++;
		}
	}
	canvasFlow->SetTextInfo("#font[42]{v_{n}}","Y",0.25,0.6,0.5,0.0);
	canvasFlow->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %s}}", particle[GetNumberParticle(PID,"")].Data()),"Info",1.,0.85*max_y_axis,0.13,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{ STAR Preliminary}}"),"Prel",0.4,0.85*max_y_axis,0.13,0.0);
	canvasFlow->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasFlow->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.13);
	canvasFlow->SetDrawObject(flow[0],"OnePad");
	canvasFlow->SetLegend(0.12,0.5,0.98,0.82,3);
	TCanvas *result = canvasFlow->CanvasNxM(1000,1280, 5, 2, 0.4, 0.575, flow,1,0,0);

	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/sys/SysBES_vn_%s_%s_pt_energy_0_60.pdf",outPath, PID.Data(), EtaGap.Data()));
	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/sys/SysBES_vn_%s_%s_pt_energy_0_60.png",outPath, PID.Data(), EtaGap.Data()));
	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/sys/SysBES_vn_%s_%s_pt_energy_0_60.C",outPath, PID.Data(), EtaGap.Data()));
	delete result;

	for(Int_t i=0; i<6; i++){
		file_1[i]->Close();
		delete file_1[i];
	}

}

void HadronsFlowVsPtVsCentForBES_ver2Int( TString EtaGap, TString PID, Double_t min_x_axis, Double_t max_x_axis , Double_t min_y_axis, Double_t max_y_axis){

	VEC rebin_vec;
	VEC rebin_vec_integral={0.2,3.2,5.2};
	TString prefix = "";
	
	Int_t Arr_energy[6]={11,14,19,27,39,62};
	const Int_t color[]={6, 1, 4, 2, 1, 2, 46};
	const Int_t style[]={8, 22, 23, 24, 25, 34, 8, 29};
	const Double_t size[]={1.5,1.5,1.5,1.5,1.5,1.5,1.5};
	
	Int_t centBinMin[6] = {7,6,5,4,2};
	Int_t centBinMax[6] = {8,6,5,4,3};

	Float_t k[2]={1.0,1.0};

	TFile *file_1[6];
	for(Int_t i=0; i<6; i++){
		file_1[i] = new TFile(Form("%s/file_hadrons/Sys/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2_sys.root",outPath, Arr_energy[i]),"READ");
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

				rebin_vec = Rebin2( harm+2, Arr_energy[j], PID, "", centBinMin[c], centBinMax[c]);
				
				object.push_back(new Draw_Object);
				Analysis -> SetParametrs(Arr_energy[j], file_1[j], file_out, harm+2, PID, "", EtaGap, prefix);
				//Analysis -> FlowVsPt_EtaSub(object.back(), centBinMin[c], centBinMax[c], rebin_vec);
	    		Analysis -> FlowVsPtDivideIntegralFLow_EtaSub(object.back(), centBinMin[c], centBinMax[c], rebin_vec, rebin_vec_integral, 0);
	    		object.back() -> SetParametrsGraph(Form("data_v%i_%iGeV_Hadrons",harm+2,Arr_energy[j]),Form("#font[42]{#scale[1.3]{%s GeV}}", Energy_scan_text.at(Arr_energy[j])), color[j], style[j], size[j]);
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
	}
	canvasFlow->SetTextInfo("#font[42]{#frac{v_{n}}{v_{n}^{int}}}","Y",0.25,0.6,0.5,0.0);
	canvasFlow->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %s}}", particle[GetNumberParticle(PID,"")].Data()),"Info",1.,0.85*max_y_axis,0.13,0.0);
	canvasFlow->SetTextInfo(Form("#font[42]{ #scale[1.0]{ STAR Preliminary}}"),"Prel",.4,0.85*max_y_axis,0.13,0.0);
	canvasFlow->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasFlow->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.13);
	canvasFlow->SetDrawObject(flow[0],"OnePad");
	canvasFlow->SetLegend(0.12,0.5,0.98,0.82,3);
	TCanvas *result = canvasFlow->CanvasNxM(1000,1280, 5, 2, 0.4, 0.575, flow,1,0,0);

	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/sys/intSysBES_vn_%s_%s_pt_energy_0_60.pdf",outPath, PID.Data(), EtaGap.Data()));
	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/sys/intSysBES_vn_%s_%s_pt_energy_0_60.png",outPath, PID.Data(), EtaGap.Data()));
	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/sys/intSysBES_vn_%s_%s_pt_energy_0_60.C",outPath, PID.Data(), EtaGap.Data()));
	delete result;

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
		file_1[i] = new TFile(Form("%s/file/Flow_%iGeV_PID_10binYesTofdca11_oldPID_BadRunOld.root",outPath, Arr_energy[i]),"READ");
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
	    		object.back() -> SetParametrsGraph(Form("data_v%i",harm+2),Form("#font[42]{#scale[1.2]{%s GeV}}", Energy_scan_text.at(Arr_energy[j])), color[j], style[j], size[j]);

				if(harm==1){
					for(Int_t n=0; n < object.back()->GetSizeVector(); n++){
						object.back()->СhangePointGraph(n,
														object.back()->GetPointXGraph(n), 
														object.back()->GetPointYGraph(n) * k[1],
														object.back()->GetPointXErrorGraph(n),
														object.back()->GetPointYErrorGraph(n)* k[1],
														object.back()->GetPointYSysErrorGraph(n)* k[1]);
					}
				}

				if(j==5 && par==2){
					if(harm==0) canvasFlow->SetTextInfo(Form("#font[42]{v_{%i}}",harm+2), Form("V%i",harm+2), 0.6*max_x_axis, 0.85*max_y_axis, 0.07, 0.0);
					if(harm==1) canvasFlow->SetTextInfo(Form("#font[42]{%.1f #times v_{%i}}", k[harm],harm+2), Form("V%i",harm+2), 0.6*max_x_axis, 0.85*max_y_axis, 0.07, 0.0);
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
	canvasFlow->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasFlow->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.07);
	canvasFlow->SetDrawObject(flow[0],"OnePad");
	canvasFlow->SetLegend(0.1,0.75,0.98,0.98,3);
	TCanvas *result = canvasFlow->CanvasNxM(1200,1280, 3, 2, 0.4, 0.35, flow,1,0,2);

	result->SaveAs(Form("%s/picture/3_PID/multPad_BES_vn_%s_%s_cent_%i_%i.pdf",outPath, charge.Data(),EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result->SaveAs(Form("%s/picture/3_PID/multPad_BES_vn_%s_%s_cent_%i_%i.png",outPath, charge.Data(),EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result;

	for(Int_t i=0; i<6; i++){
		file_1[i]->Close();
		delete file_1[i];
	}

}

void ParAntFlowVsPtForBESupdate( Int_t CentBinMin, Int_t CentBinMax, TString EtaGap, Double_t min_x_axis, Double_t max_x_axis){

	VEC rebin_vec;
	Int_t point = 1; 
	TString prefix = "TPCandTOF";

	Int_t Arr_energy[6]={11,14,19,27,39,62};
	const Int_t color[]={2, 1, 4, 6, 4, 2, 46};
	const Int_t style[]={24, 22, 23, 8, 25, 34, 8, 29};
	const Double_t size[]={1.3,1.4,1.3,1.4,1.3,1.3,1.3};
	TString arr_par[]={"Pion","Kaon","Proton"};

	TFile *file_1[7];
	for(Int_t i=0; i<6; i++){
		file_1[i] = new TFile(Form("%s/file/Flow_%iGeV_PID_10binYesTofdca11_oldPID_BadRunOld.root",outPath, Arr_energy[i]),"READ");
	}
	TFile *file_out = new TFile("OUT.root","RECREATE");

	Draw_Picture_new *canvasDif = new Draw_Picture_new;
	Draw_Picture_new *canvasRatio = new Draw_Picture_new;
	
	Analysis_Function *Analysis = new Analysis_Function;
	std::vector<std::vector<Draw_Object *>> flow_ratio;
	std::vector<std::vector<Draw_Object *>> flow_difference;
	std::vector<Draw_Object *> object;
	
	Double_t min_y_axis_ratio[3]={0.68,0.68,0.8};
	Double_t max_y_axis_ratio[3]={1.33,1.38,2.23};
	Double_t min_y_axis_dif[3]={-0.022, -0.022, -0.012};
	Double_t max_y_axis_dif[3]={0.017, 0.022, 0.044};

	Double_t coord_y_draw_dif[3];
	Double_t coord_y_draw_ratio[3];

	Int_t b=0;
	Int_t b2=0;

	for(Int_t par=0; par<3; par++){

		coord_y_draw_dif[par]=min_y_axis_dif[par] + (max_y_axis_dif[par] - min_y_axis_dif[par])*0.9;
		coord_y_draw_ratio[par]=min_y_axis_ratio[par] + (max_y_axis_ratio[par] - min_y_axis_ratio[par])*0.95;
		
		for(Int_t harm=0; harm<2; harm++){
		
			for(Int_t j=0; j<6; j++){

				rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], "", CentBinMin, CentBinMax);
	    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"",EtaGap,prefix);
				object.push_back(new Draw_Object);
				Analysis -> DifferentParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i",harm+2),Form("#font[42]{#scale[1.2]{%s GeV}}", Energy_scan_text.at(Arr_energy[j])), color[j], style[j], size[j]);
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
		
			for(Int_t j=0; j<6; j++){

				rebin_vec = Rebin2( harm+2, Arr_energy[j], arr_par[par], "", CentBinMin, CentBinMax);
	    		Analysis -> SetParametrs(Arr_energy[j], file_1[j],file_out,harm+2,arr_par[par],"",EtaGap,prefix);
				object.push_back(new Draw_Object);
				Analysis -> RatioParAnt(object.back(), CentBinMin, CentBinMax, rebin_vec);
	    		object.back() -> SetParametrsGraph(Form("data_v%i",harm+2),Form("#font[42]{#scale[1.2]{%s GeV}}", Energy_scan_text.at(Arr_energy[j])), color[j], style[j], size[j]);
			}
			flow_ratio.push_back(object);
			object.clear();
			if(par==2){
				canvasRatio->SetTextInfo(Form("#font[42]{v_{%i}}", harm+2), Form("V%i",harm+2), 0.85*max_x_axis,0.98*coord_y_draw_ratio[0], 0.08, 0.0);
			}
			canvasRatio->SetNumberPad(Form("#font[42]{%s}", bukva[b2].Data()), "Pad", 0.04, 0.98*coord_y_draw_ratio[par], 0.08, 0.0);
			b2++;
		}
	}

	canvasDif->SetTextInfo(Form("#font[42]{ v_{n}(%s) - v_{n}(%s) }",particle[GetNumberParticle("Pion","Pos")].Data(),particle[GetNumberParticle("Pion","Neg")].Data() ),"Y",0.4, 0.7,0.4,90);
	canvasDif->SetTextInfo(Form("#font[42]{ v_{n}(%s) - v_{n}(%s) }",particle[GetNumberParticle("Kaon","Pos")].Data(),particle[GetNumberParticle("Kaon","Neg")].Data() ),"Y",0.4, 0.4,0.4,90);
	canvasDif->SetTextInfo(Form("#font[42]{ v_{n}(%s) - v_{n}(%s) }",particle[GetNumberParticle("Proton","Pos")].Data(),particle[GetNumberParticle("Proton","Neg")].Data() ),"Y",0.4, 0.07,0.4,90);
	canvasDif->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasDif->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %i-%i %% }}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.6,0.95*coord_y_draw_dif[0],0.08,0.0);
	canvasDif->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasDif->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasDif->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_dif[0], max_y_axis_dif[0], 0.06);
	canvasDif->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_dif[1], max_y_axis_dif[1], 0.06);
	canvasDif->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_dif[2], max_y_axis_dif[2], 0.06);
	canvasDif->SetDrawObject(flow_difference[0],"OnePad");
	canvasDif->SetLegend(0.07,0.08,0.98,0.2,3);
	TCanvas *result = canvasDif->CanvasNxM(1000,1280, 3, 2, 0.4, 0.35, flow_difference,3,0,1);

	result->SaveAs(Form("%s/picture/3_PID/BESupdateDif_vn_%s_cent_%i_%i.pdf",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result->SaveAs(Form("%s/picture/3_PID/BESupdateDif_vn_%s_cent_%i_%i.png",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result;

	canvasRatio->SetTextInfo(Form("#font[42]{#frac{ v_{n}(%s) }{ v_{n}(%s) } }", particle[GetNumberParticle("Pion","Pos")].Data(), particle[GetNumberParticle("Pion","Neg")].Data()),"Y",0.07,0.8,0.3,0);
	canvasRatio->SetTextInfo(Form("#font[42]{#frac{ v_{n}(%s) }{ v_{n}(%s) } }", particle[GetNumberParticle("Kaon","Pos")].Data(), particle[GetNumberParticle("Kaon","Neg")].Data()),"Y",0.07,0.5,0.3,0);
	canvasRatio->SetTextInfo(Form("#font[42]{#frac{ v_{n}(%s) }{ v_{n}(%s) } }", particle[GetNumberParticle("Proton","Pos")].Data(), particle[GetNumberParticle("Proton","Neg")].Data()),"Y",0.07,0.15,0.3,0);
	canvasRatio->SetTextInfo("#font[42]{p_{T} [GeV/c]}","X",0.75,0.6,0.4,0.0);
	canvasRatio->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %i-%i %% }}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.6,0.98*coord_y_draw_ratio[0],0.08,0.0);
	canvasRatio->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasRatio->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_ratio[0], max_y_axis_ratio[0], 0.06);
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_ratio[1], max_y_axis_ratio[1], 0.06);
	canvasRatio->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis_ratio[2], max_y_axis_ratio[2], 0.06);
	canvasRatio->SetDrawObject(flow_ratio[0],"OnePad");
	canvasRatio->SetLegend(0.04,0.7,0.98,0.85,3);
	TCanvas *result2 = canvasRatio->CanvasNxM(1000,1280, 3, 2, 0.4, 0.35, flow_ratio,3,0,0);

	result2->SaveAs(Form("%s/picture/3_PID/BESupdateRatio_vn_%s_cent_%i_%i.pdf",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	result2->SaveAs(Form("%s/picture/3_PID/BESupdateRatio_vn_%s_cent_%i_%i.png",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result2;


	for(Int_t i=0; i<6; i++){
		file_1[i]->Close();
		delete file_1[i];
	}

}


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
	
	Int_t Arr_energy[6]={11,14,19,27,39,62};
	Double_t x_energy[6]={11.5,14.5,19.6,27,39,62.4};
	TString arr_par[]={"Pion","Kaon","Proton"};

	TFile *file_1[6];
	for(Int_t i=0; i<6; i++){
		file_1[i] = new TFile(Form("%s/file/Flow_%iGeV_PID_10binYesTofdca11_oldPID_BadRunOld.root",outPath, Arr_energy[i]),"READ");
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
				object->ClearObject();
			}
	    	res.back() -> SetParametrsGraph(Form("data_v%i",harm+2),Form("#font[42]{#scale[1.1]{%s - %s}}", particle[GetNumberParticle(arr_par[par],"Pos")].Data(), particle[GetNumberParticle(arr_par[par],"Neg")].Data()), color[par], style[par], size[par]);

		}

		flow.push_back(res);
		res.clear();
		canvasDif->SetNumberPad(Form("#font[42]{%s}", bukva[b].Data()), "Pad", 0.04, 0.85*max_y_axis, 0.06, 0.0);
		canvasDif->SetTextInfo(Form("#font[42]{n = %i}", harm+2), Form("V%i",harm+2), 0.85*max_x_axis,0.85*max_y_axis, 0.06, 0.0);
			
		b++;
	}


	canvasDif->SetTextInfo(Form("#font[42]{ v_{n}(X) - v_{n}(#bar{X}) }"),"Y",0.4, 0.4,0.3,90);
	canvasDif->SetTextInfo("#font[42]{#sqrt{s_{NN}} [GeV]}","X",0.5,0.4,0.6,0.0);
	canvasDif->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Au+Au, %i-%i %% }}", cent[CentBinMax+1], cent[CentBinMin]),"Info",0.1*max_x_axis,0.85*max_y_axis,0.06,0.0);
	canvasDif->SetTLine(min_x_axis, 0., max_x_axis, 0., 1, 2, 2, "Up");
	canvasDif->SetTLine(min_x_axis, 1., max_x_axis, 1., 1, 2, 2, "Up");
	canvasDif->SetAxisToCanvsNM(min_x_axis, max_x_axis,  min_y_axis, max_y_axis, 0.06);
	canvasDif->SetDrawObject(flow[0],"OnePad");
	canvasDif->SetLegend(0.1,0.7,0.4,0.98,1);
	TCanvas *result = canvasDif->CanvasNxM(1000,650, 1, 2, 0.4, 0.15, flow,1,0,1);

	result->SaveAs(Form("%s/picture/3_PID/ParAntVsEnergy_vn_%s_cent_%i_%i_pt_%1.f_%1.f.pdf",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin], ptMin, ptMax));
	result->SaveAs(Form("%s/picture/3_PID/ParAntVsEnergy_vn_%s_cent_%i_%i_pt_%1.f_%1.f.png",outPath, EtaGap.Data(), cent[CentBinMax+1], cent[CentBinMin], ptMin, ptMax));
	delete result;

	std::cout<<"good\n";

	for(Int_t i=0; i<6; i++){
		file_1[i]->Close();
		delete file_1[i];
	}
}




void FlowDifferentMethodEPandSP(Int_t Energy, Int_t harmonic, Int_t CentBinMin, Int_t CentBinMax, TString EtaGap,TString PID, TString charge, Double_t max_y_axis, Double_t max_x_axis, Double_t MinYratio, Double_t MaxYratio){

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
	std::vector<TString> NameEtaPID = {"Eta1","Eta04","Eta07","Eta10"};
	std::vector<Double_t> EtaVecPID = {0.05, 0.2, 0.35, 0.5};

	const Int_t color[]={2, 1, 4, 6, 4, 2, 46};
	const Int_t style[]={25, 23, 22, 23, 34, 47, 8, 29};
	const Double_t size[]={2.0,2.0,2.0,2.0,2.0,2.0,1.3};

	TFile *file_read;
	TFile *file_read_SP;

	if( strncmp(PID, "Hadrons",7)==0){
		file_read = new TFile(Form("%s/file_hadrons/Flow_%iGeV_Hadrons_2binYesTof_dcaEP2.root",outPath,Energy),"READ");	
		file_read_SP = new TFile(Form("%s/file_hadrons/Flow_%iGeV_Hadrons_SP.root",outPath,Energy),"READ");	
    }else{
		file_read = new TFile(Form("%s/file/Flow_%iGeV_PID_10binYesTofdca11_oldPID_BadRunOld.root",outPath, Energy),"READ");
    }

	TFile *file_out = new TFile("OUT.root","RECREATE");

	flow.push_back(new Draw_Object);
	Analysis -> SetParametrs(Energy, file_read,file_out,harmonic,PID,charge,NameEtaHadrons[0],prefix);
    flow.back() -> SetParametrsGraph(Form("data%.3f",EtaVecHadrons[0]),Form("#font[42]{#eta-sub EP}"), color[3], style[3], size[3]);
    Analysis -> FlowVsPt_EtaSub(flow.back(), CentBinMin, CentBinMax, rebin_vec);
	
	TString axis[]={"","X","Y"};
	for(Int_t eta=0; eta<3; eta++){
	    flow.push_back(new Draw_Object);
    	if(strncmp(PID, "Hadrons",7)==0){
    		Analysis -> SetParametrs(Energy, file_read_SP,file_out,harmonic,PID,charge,EtaGap,prefix);
    		flow.back() -> SetParametrsGraph(Form("data%.3f",EtaVecHadrons[0]),Form("#font[42]{SP, %s}",axis[eta].Data()), color[eta], style[eta], size[eta]);
    	}else{
    		Analysis -> SetParametrs(Energy, file_read_SP,file_out,harmonic,PID,charge,NameEtaPID[0],prefix);
    		flow.back() -> SetParametrsGraph(Form("data%.3f",EtaVecPID[0]),Form("#font[42]{SP, %s}", axis[eta].Data()), color[eta], style[eta], size[eta]);
    	}
    	Analysis -> FlowVsPt_SP(flow.back(), CentBinMin, CentBinMax, rebin_vec, axis[eta]);
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
													"#font[42]{p_{T} [GeV/c]}",Form("v_{%i}",harmonic),"#font[42]{#frac{SP}{EP}}","");
	result->SaveAs(Form("%s/picture/2_Charge_Hadrons/met/DifferentMetgod_%iGeV_v%i_pt_%s%s_%i_%i.png",outPath, Energy,harmonic,PID.Data(), charge.Data(), cent[CentBinMax+1], cent[CentBinMin]));
	delete result;
	file_read->Close();
	file_out->Close();

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
	
	FlowVsEnergyBESupdate( 0.2, 3.2, "Eta15","Hadrons", "", 0.0, 0.115, 75, 0.5, 1.15);
	FlowVsEnergyBESupdate( 0.2, 2.0, "Eta15","Hadrons", "", 0.0, 0.115, 75, 0.5, 1.15);

}

void DrawFlowDifferentEtaGap(TString mod, Int_t Energy){
	if(strncmp(mod, "Hadrons",7)==0){
		FlowDifferentEtaGap(Energy, 2, 2, 6,"Hadrons","",0.28, 2.93,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 7, 8,"Hadrons","",0.28, 2.93,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 2, 8,"Hadrons","",0.28, 2.93,0.89,1.11);
		FlowDifferentEtaGap(Energy, 3, 2, 6,"Hadrons","",0.15, 2.93,0.8,1.2);
		FlowDifferentEtaGap(Energy, 3, 7, 8,"Hadrons","",0.15, 2.93,0.8,1.2);
		FlowDifferentEtaGap(Energy, 3, 2, 8,"Hadrons","",0.15, 2.93,0.8,1.2);
		FlowDifferentEtaGap(Energy, 3, 3, 8,"Hadrons","",0.15, 2.93,0.8,1.2);
		
		FlowDifferentEtaGap(Energy, 3, 8, 8,"Hadrons","",0.15, 2.93,0.8,1.2);
		FlowDifferentEtaGap(Energy, 3, 7, 7,"Hadrons","",0.15, 2.93,0.8,1.2);
		FlowDifferentEtaGap(Energy, 3, 6, 6,"Hadrons","",0.15, 2.93,0.8,1.2);
		FlowDifferentEtaGap(Energy, 3, 5, 5,"Hadrons","",0.15, 2.93,0.8,1.2);
		FlowDifferentEtaGap(Energy, 3, 4, 4,"Hadrons","",0.15, 2.93,0.8,1.2);
		FlowDifferentEtaGap(Energy, 3, 3, 3,"Hadrons","",0.15, 2.93,0.8,1.2);
		FlowDifferentEtaGap(Energy, 3, 2, 2,"Hadrons","",0.15, 2.93,0.8,1.2);

		FlowDifferentEtaGap(Energy, 2, 8, 8,"Hadrons","",0.287, 2.93,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 7, 7,"Hadrons","",0.287, 2.93,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 6, 6,"Hadrons","",0.287, 2.93,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 5, 5,"Hadrons","",0.287, 2.93,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 4, 4,"Hadrons","",0.287, 2.93,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 3, 3,"Hadrons","",0.287, 2.93,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 2, 2,"Hadrons","",0.287, 2.93,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 1, 1,"Hadrons","",0.287, 2.93,0.89,1.11);
		FlowDifferentEtaGap(Energy, 2, 0, 0,"Hadrons","",0.287, 2.93,0.89,1.11);
	}else{

	}
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

VEC Rebin2(Int_t harmonic, Int_t Energy, TString PID, TString charge, Int_t CentBinMin, Int_t CentBinMax){

	VEC rebin_vector;

	/*
	if(strncmp(PID, "Hadrons",7)==0 ){
		rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 10.0};
		//rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8 , 3.2, 10.0};
	}

	*/
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
	*/

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

	if(harmonic==2){

		if((CentBinMin==2 && CentBinMax==8) || (CentBinMin==2 && CentBinMax==6) ){ //0-60
			
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
	}

	if(harmonic==3){

		if( (CentBinMin==2 && CentBinMax==8) || (CentBinMin==2 && CentBinMax==6)){ //0-60
			
			if(Energy>25){
				if(strncmp(PID, "Pion",4)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.6, 2.2, 2.6, 5.0};
				}
				if(strncmp(PID, "Kaon",4)==0){
					rebin_vector = {0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.8, 5.0};
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
					rebin_vector = {0.2, 0.4, 0.8, 1.2, 1.6, 2.4, 5.0};
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
					rebin_vector = {0.2, 0.6, 1.0, 1.4, 2.0, 5.0};
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
					rebin_vector = {0.2, 0.6, 1.0, 1.6, 2.2, 5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.5, 0.8, 1.4, 2.6, 10.0};
				}
			}




			//////////////////////

			if(Energy==27){
				if(strncmp(PID, "Pion",4)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 5.0};
				}
				if(strncmp(PID, "Kaon",4)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.5, 0.8, 1.2, 1.8, 2.6, 10.0};
				}
			}

			//////////////

		}


		if(CentBinMin==7 && CentBinMax==8){ //0-10
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.6, 1.0, 1.4, 1.8, 2.4, 5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.6, 1.0, 1.4, 2.2, 2.8, 10.0};
			}
		}
		if(CentBinMin==2 && CentBinMax==6){ //10-60
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.8,5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.5, 1.0, 1.4, 2.0, 10.0};
			}
		}
		
		if(CentBinMin==4 && CentBinMax==6){ //10-60
			if(strncmp(PID, "Pion",4)==0 ||  strncmp(PID, "Kaon",4)==0 ){
				rebin_vector = {0.2, 0.4, 0.6, 1.0, 1.4, 1.8, 2.2, 2.8,5.0};
			}
			if(strncmp(PID, "Proton",6)==0 ){
				rebin_vector = {0.2, 0.5, 1.0, 1.4, 2.0, 10.0};
				if(Energy>15){
					rebin_vector = {0.2, 0.5, 1.0, 1.4, 2.0, 2.6, 10.0};
				}
				if(Energy==27){
				if(strncmp(PID, "Pion",4)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 5.0};
				}
				if(strncmp(PID, "Kaon",4)==0){
					rebin_vector = {0.2, 0.4, 0.6, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 5.0};
				}
				if(strncmp(PID, "Proton",6)==0 ){
					rebin_vector = {0.2, 0.5, 0.8, 1.2, 1.8, 2.6, 10.0};
				}
			}
			}
		}
		if(strncmp(PID, "Hadrons",7)==0 ){
			if(Energy==7){
				rebin_vector = {0.2, 0.4, 0.8, 1.2, 1.6, 2.2, 2.8, 10.0};
			}
			if(Energy==11){
					//rebin_vector = {0.2, 0.4, 0.6, 1.0, 1.4, 2.0, 2.6, 10.0};

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